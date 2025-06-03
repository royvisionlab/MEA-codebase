
import sys
sys.path.append('../analysis/')

import numpy as np
import os
import platform
import torch
import pickle
import visionloader as vl
import visionwriter as vw

from vision_utils import ParamsWriter, STAWriter
from cell_typing import cluster_compact, get_cluster_labels, write_cluster_labels

from progressbar import *
import argparse
from symphony_data import Stimulus
import lnp
import gc
import pickle

def crop_frames(frames: np.ndarray, crop_fraction: float=1.0) -> np.ndarray:
    """ Crops the stimulus frames by a fraction of the original size.
    Parameters:
        frames: stimulus tensor of size t, y, x, color (float32 in contrast [-1,1]).
        crop_fraction: fraction by which to crop the stimulus.
    Returns:
        frames: stimulus tensor of size t, y, x, color (float32 in contrast [-1,1]).
    """
    try:
        if crop_fraction < 1.0:
            _, height, width, _ = frames.shape
            x_size = int(width * crop_fraction)
            x_start = int((width - x_size) / 2)
            x_end = x_start + x_size
            y_size = int(height * crop_fraction)
            y_start = int((height - y_size) / 2)
            y_end = y_start + y_size
            frames = frames[:, y_start:y_end, x_start:x_end, :]
        return frames
    except Exception as error:
        print('An error occurred while cropping the frames:', type(error).__name__, '-', error)
        return frames

# Minimal STA analysis.
class STAAnalysisSingle(object):
    def __init__(self, 
            data_path: str,
            file_name: str,
            protocol_id: str='FastNoise',
            sort_algorithm: str='kilosort2',
            save_responses: bool=True):
        self.data_path = data_path
        self.S = Stimulus()
        self.sort_algorithm = sort_algorithm
        self.protocol_id = protocol_id
        self.save_responses = save_responses
        # self.vcd = vl.load_vision_data( os.path.join(data_path,file_name), file_name, include_neurons=True )
        self.vcd = vl.load_vision_data( data_path, file_name, include_neurons=True )
        # Get the cell Ids.
        self.cell_ids = sorted(self.vcd.get_cell_ids())

    def analyze(self, 
                    parameters: dict,
                    start_seed: int,
                    frame_rate: float=60.31807657,
                    stride: int=1,
                    depth: int=31,
                    crop_fraction: float=1.0):
        
        stas = dict()
        spike_counts = dict()
        spike_counts_blue = dict() # Spike counts for blue frames.

        epoch_starts = self.vcd.get_ttl_times()
        print(epoch_starts)

        # Loop through each epoch.
        for count, start in enumerate(epoch_starts):
            
            # Need to get rid of pre/post frames.
            preFrames = round(parameters['preTime'] * 1e-3 * frame_rate)
            stimFrames = round(parameters['stimTime'] * 1e-3 * frame_rate)

            frame_times = np.round(start/20000*1000 + np.arange(stimFrames+1)*1000/frame_rate)

            # Get the frames when the stimulus was on.
            idx = (frame_times-frame_times[0]).searchsorted((parameters['preTime']+parameters['stimTime']),'right')
            frame_times = frame_times[preFrames:idx+1]

            if (self.protocol_id == 'FastNoise'):
                epoch_count = (count % len(parameters['numXStixels']))
                frames = self.S.getFastNoiseFrames(parameters['numXStixels'][epoch_count],
                    parameters['numYStixels'][epoch_count], 
                    parameters['numXChecks'],
                    parameters['numYChecks'],
                    parameters['chromaticClass'],
                    len(frame_times)-1,
                    parameters['stepsPerStixel'][epoch_count], 
                    start_seed + np.floor(count/2).astype(int), #start_seed + count,
                    parameters['frameDwell'][epoch_count] )
            elif (self.protocol_id == 'PinkNoise'):
                frames = self.S.getPinkNoiseFrames( 
                            parameters['numXChecks'], 
                            parameters['numYChecks'],
                            len(frame_times)-1,
                            parameters['noiseContrast'], 
                            parameters['spatialAmplitude'],
                            parameters['temporalAmplitude'],
                            parameters['chromaticClass'],
                            start_seed + count)
            if crop_fraction < 1.0:
                frames = crop_frames(frames=frames, crop_fraction=crop_fraction)
            if stride > 1:
                frames = self.S.upsample_frames(frames, stride)

            # Get the binned spike count.
            binned_spikes = self.get_binned_spikes(frame_times, stride)

            if self.save_responses:
                if count == 0:
                    stimulus = np.zeros((len(epoch_starts),frames.shape[0],frames.shape[1],frames.shape[2],frames.shape[3]))
                    response = np.zeros((len(epoch_starts),binned_spikes.shape[0],binned_spikes.shape[1]))
                stimulus[count,:,:,:,:] = frames
                response[count,:,:] = binned_spikes

            # Compute the STA.
            # stride = 1 # Equivalent to bins/frame
            sta_tensor = self.compute_sta_torch(frames, binned_spikes, 1, depth)

            # Get the total spike count.
            sp_count = np.sum(binned_spikes, axis=0)

            for cell_count, cell in enumerate(self.cell_ids):
                if cell in stas:
                    stas[cell] += sta_tensor[cell_count,:,:,:,:]
                    spike_counts[cell] += sp_count[cell_count]
                    if (parameters['chromaticClass'] != 'achromatic'):
                        spike_counts_blue[cell] += sp_count[cell_count]
                else:
                    stas[cell] = sta_tensor[cell_count,:,:,:,:]
                    spike_counts[cell] = sp_count[cell_count]
                    if (parameters['chromaticClass'] == 'achromatic'):
                        spike_counts_blue[cell] = 0
                    else:
                        spike_counts_blue[cell] = sp_count[cell_count]

        # Sort the STA dictionary by cluster.
        stas = dict(sorted(stas.items()))

        # Get the STA by dividing by the spike count. 
        for cell in stas:
            if spike_counts[cell] > 1:
                stas[cell][:,:,:,:2] /= spike_counts[cell]

            if spike_counts_blue[cell] > 1:
                stas[cell][:,:,:,2] /= spike_counts_blue[cell]
            
            # Normalize to values between 0-1.
            if spike_counts[cell] > 1:
                stas[cell] /= np.amax(np.abs(stas[cell]))

        cluster_id = list(stas.keys())

        sta = np.zeros((len(cluster_id),stas[cluster_id[0]].shape[0],stas[cluster_id[0]].shape[1],stas[cluster_id[0]].shape[2],stas[cluster_id[0]].shape[3]))
        spike_count = np.zeros((len(cluster_id),3)) # Spike counts for the three channels
        for count, cell in enumerate(stas):
            sta[count,:,:,:,:] = stas[cell]
            spike_count[count,0] = spike_counts[cell]
            spike_count[count,1] = spike_counts[cell]
            spike_count[count,2] = spike_counts_blue[cell]
        
        # if self.save_responses:
        #     d_path = os.path.join(self.data_path, 'stimulus_response.p')
        #     with open(d_path, 'wb') as f:
        #         pickle.dump({'stimulus':stimulus,'response':response}, f)
        
        return sta, np.array(cluster_id), spike_count
    
    def get_binned_spikes(self, frame_times, stride=1, sample_rate: float=20000.0):
        '''
        Parameters:
            vcd: Vision data table object
            cells: list of cells to compute STAs for
        Returns:
            A matrix of size frames by cells

        Bins the spike train using the (full resolution) frames of the monitor 
        as bin edges.
        '''

        # If the stride/binsPerFrame is greater than 1, interpolate the frame times.
        if stride > 1:
            frame_times = lnp.get_frame_times(frame_times, stride)

        # Loop through cells and bin spike train.
        binned_spikes = []

        for cell in self.cell_ids:
            spike_times = self.vcd.get_spike_times_for_cell(cell) / sample_rate * 1000 # ms
            binned_spikes.append(np.histogram(spike_times,bins=frame_times)[0])
        
        return np.asarray(binned_spikes).T

    def compute_sta_torch(self, stimulus, binned_spikes, stride, depth=31):
        """
        Computes STA by convolving the spiketrain with the stimulus. This
        implementation uses PyTorch. It doesn't use GPUs because of the memory
        intensiveness, but it is extremely fast relative to numpy even on CPU.

        Parameters:
            stimulus: stimulus tensor of size t, y, x, color (float32 in contrast [-1,1]).
            binned_spikes: spiketrain matrix of size t,cells
            stride: the stride of convolution (interval)
            depth: number of lags for convolution
            norm: boolean whether to normalize STA for nice plotting
        Returns:
            STA tensor, indexed by cell,t,y,x,color channel
        """

        # Intialize STA tensor.
        num_cells = binned_spikes.shape[1]
        n_frames,height,width,n_phosphors = stimulus.shape
        sta_tensor = np.zeros((num_cells,depth,height,width,n_phosphors))

        # Vectorize stimulus.
        stimulus = np.reshape(stimulus,(n_frames,
                                        height * width * n_phosphors)).astype(np.float32).T
        
        # Initialize STA, the offset range, and compute.
        offset_cnt = 0
        sta_tensor = np.zeros((num_cells,depth,height,width,n_phosphors))
        offset_range = np.arange(0,(depth * stride),stride)

        # Initialize progress bar.
        widgets = ['Computing STAs: ', Percentage(), ' ', Bar(marker='*',
                left='[',right=']'),' ', ETA()]
        pbar = ProgressBar(widgets=widgets, maxval=depth)
        pbar.start()

        # Set up torch stuff.
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu') #torch.device('cpu')
        
        binned_spikes = torch.tensor(binned_spikes).to(device).to(torch.float32)
        stimulus = torch.tensor(stimulus).to(device).to(torch.float32)

        for offset in offset_range:

            with torch.no_grad():
                sta = torch.matmul(stimulus[:,:n_frames-offset],
                                binned_spikes[offset:,:])

                # Index appropriately to avoid nans and normalize by number of spikes.
                nonzero_inds = torch.nonzero(torch.sum(binned_spikes[offset:,:],
                                                    dim=0),as_tuple=False)
                sta[:,nonzero_inds] /= torch.sum(binned_spikes[offset:,:],
                                                dim=0)[nonzero_inds]

            sta = sta.detach().cpu().numpy() 
            sta_tensor[:,offset_cnt,...] = np.reshape(sta.T,
                                                    (num_cells,height,
                                                    width,n_phosphors))
            pbar.update(offset_cnt)
            offset_cnt +=1

        pbar.finish()

        if device.type == 'cuda':
            torch.cuda.empty_cache()
        
        # Clear some memory
        del sta
        gc.collect()

        return sta_tensor

    def get_interspike_intervals(self, bin_edges=np.linspace(0,200,401)):
        isi = dict()
        for cell in self.cell_ids:
            spike_times = self.vcd.get_spike_times_for_cell(cell) / 20000 * 1000 # ms
            
            # Compute the interspike interval
            if len(spike_times) > 1:
                isi_tmp = np.diff(spike_times)
                isi[cell] = np.histogram(isi_tmp,bins=bin_edges)[0]
            else:
                isi[cell] = np.zeros((len(bin_edges)-1,)).astype(int)
        
        return isi

    def get_autocorrelation(self, bin_edges=np.linspace(0,200,401)):
        acf_dict = dict() # Autocorrelation dictionary

        # Compute the autocorrelation function.
        acf_tmp = self.get_interspike_intervals(bin_edges=bin_edges)

        for cell in self.cell_ids:
            if cell in acf_dict:
                acf_dict[cell] += acf_tmp[cell].astype(int)
            else:
                acf_dict[cell] = acf_tmp[cell].astype(int)

        # Sort the ACF dictionary by cluster Id.
        acf_dict = dict(sorted(acf_dict.items()))
        cluster_id = list(acf_dict.keys())

        isi = np.zeros((len(cluster_id),acf_dict[cluster_id[0]].shape[0]))
        acf = np.zeros((len(cluster_id),acf_dict[cluster_id[0]].shape[0]))
        for count, cell in enumerate(acf_dict):
            isi[count,:] = acf_dict[cell]
            if np.sum(isi[count,:]) > 1:
                acf[count,:] = isi[count,:] / np.sum(isi[count,:])
            else:
                acf[count,:] = isi[count,:]

        return acf, isi, np.array(cluster_id)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Convert Kilosort output to Vision .neurons format')
    parser.add_argument('data_path', type=str, help='folder containing Kilosort outputs, must already exist')
    parser.add_argument('file_name', type=str, help='name of data file/folder (e.g. data020)')
    parser.add_argument('start_seed', type=int, help='seed of first epoch')
    parser.add_argument('-a','--algorithm', default='kilosort2', type=str, help='Sorting algorithm used (yass or kilosort2)')
    parser.add_argument('-p','--protocol_id', default='FastNoise', type=str, help='Name of Symphony stimulus protocol (default: FastNoise)')
    parser.add_argument('-g','--group', default=None, type=str, help='path of original dataset, needed to determine trigger position')
    parser.add_argument('-b', '--stride', default=1, type=int, help='path to write Vision files, must already exist')
    parser.add_argument('-d','--depth', default=20, type=int, help='depth of STA (default: 20)')
    parser.add_argument('-e','--experiment', type=str, help='Experiment label')
    parser.add_argument('-s','--stixel_size', default=30.0, type=float, help='Size of underlying grid in microns (default: 30.0)')
    parser.add_argument('-i', '--array_id', default=504, type=int, help='Array ID for Vision .globals file')

    args = parser.parse_args()

    if args.experiment is None:
        s = args.data_path.split('/')
        experiment_name = s[-2]
    else:
        experiment_name = args.experiment

    micronsPerStixel = args.stixel_size
    array_id = args.array_id
    data_path = args.data_path
    file_name = args.file_name
    protocol_id = args.protocol_id
    sort_algorithm = args.algorithm
    start_seed = args.start_seed
    stride = args.stride
    depth = args.depth
    stixel_size = args.stixel_size

    # Load the trigger times.
    # triggers = pickle.load( open(os.path.join(args.data_path,'litke_trigger_times.p'), 'rb') )
    # epoch_starts = triggers['trigger_times']
    # a = STAAnalysisSingle(data_path,file_name,protocol_id=protocol_id,sort_algorithm=sort_algorithm)
    a = STAAnalysisSingle(data_path, file_name, protocol_id=protocol_id, sort_algorithm=sort_algorithm)

    if args.protocol_id == 'FastNoise':
        if array_id > 3500:
            crop_fraction = 0.8
            parameters = {'numXChecks':203, 
                          'numYChecks':127, 
                          'numXStixels':[69,69], 
                          'numYStixels':[44,44],
                          'frameDwell':[1,1],
                          'stepsPerStixel':[3,3],
                          'preTime':250,'stimTime':160000,'tailTime':250,'chromaticClass':'BY'}
        elif array_id >= 1501 and array_id <= 3500:
            crop_fraction = 0.5
            parameters = {'numXChecks':304, 
                        'numYChecks':190, 
                        'numXStixels':[103,103], 
                        'numYStixels':[65,65],
                        'frameDwell':[1,2],
                        'stepsPerStixel':[3,3],
                        'preTime':250,'stimTime':160000,'tailTime':250,'chromaticClass':'BY'}
        else:
            crop_fraction = 0.7
            parameters = {'numXChecks':203, 
                          'numYChecks':127, 
                          'numXStixels':[69,69], 
                          'numYStixels':[44,44],
                          'frameDwell':[1,1],
                          'stepsPerStixel':[3,3],
                          'preTime':250,'stimTime':160000,'tailTime':250,'chromaticClass':'BY'}
        # parameters = {'numXChecks':152, 
        #               'numYChecks':95, 
        #               'numXStixels':[52,52], 
        #               'numYStixels':[33,33],
        #               'frameDwell':[1,1],
        #               'stepsPerStixel':[3,3],
        #               'preTime':250,'stimTime':280000,'tailTime':250,'chromaticClass':'BY'}
        # parameters = {'numXChecks':100, 
        #               'numYChecks':75, 
        #               'numXStixels':[51,35,51,35], 
        #               'numYStixels':[39,26,39,26],
        #               'frameDwell':[2,2,1,1],
        #               'stepsPerStixel':[2,3,2,3],
        #               'preTime':250,'stimTime':180000,'tailTime':250,'chromaticClass':'BY'}
        # parameters = {'numXChecks':200, 
        #               'numYChecks':150, 
        #               'numXStixels':[51,101,68,51,51,35], 
        #               'numYStixels':[39,76,51,39,39,26],
        #               'frameDwell':[1,2,2,6,1,1],
        #               'stepsPerStixel':[4,2,3,4,4,6],
        #               'preTime':250,'stimTime':180000,'tailTime':250,'chromaticClass':'BY'}
        # parameters = {'numXChecks':100, 'numYChecks':75, 'numXStixels':35, 'numYStixels':26,'frameDwell':1,'stepsPerStixel':3,'preTime':250,'stimTime':180000,'tailTime':250,'chromaticClass':'BY'}
        # parameters = {'numXChecks':100, 'numYChecks':76, 'numXStixels':51, 'numYStixels':39,'frameDwell':1,'stepsPerStixel':2,'preTime':250,'stimTime':180000,'tailTime':250,'chromaticClass':'BY'}
    elif args.protocol_id == 'PinkNoise':
        parameters = {'numXChecks':100, 'numYChecks':75, 'preTime':250,'stimTime':60000,'tailTime':250, 'noiseContrast':0.35, 'spatialAmplitude':0.5, 'temporalAmplitude': 0.1, 'chromaticClass':'BY'}
    else:
        parameters = {'numXChecks':152, 
                      'numYChecks':95, 
                      'numXStixels':[52,52], 
                      'numYStixels':[33,33],
                      'frameDwell':[1,1],
                      'stepsPerStixel':[1,1],
                      'preTime':250,'stimTime':180000,'tailTime':250,'chromaticClass':'B'}

    try:
        sta, cluster_id, spike_count = a.analyze(parameters=parameters, frame_rate=59.941548817817917, start_seed=start_seed, stride=stride, depth=depth, crop_fraction=crop_fraction)
    except Exception as error:
        print('An error occurred while computing the STA:', type(error).__name__, '-', error)
        sys.exit(1)
    
    print("Computed STAs")

    # Compute the ACF
    acf,isi,_ = a.get_autocorrelation()

    # Compute the spatiotemporal RF parameters.
    # timecourse_matrix, spatial_maps, gauss_params = lnp.compute_sta_params_fast(sta, snr_threshold=2.0)
    timecourse_matrix, spatial_maps, _, gauss_params, _, _, _ = lnp.compute_spatiotemporal_maps(sta)

    x0 = gauss_params[:, 0]
    y0 = sta.shape[2] - gauss_params[:, 1]
    sigma_x = gauss_params[:, 2]
    sigma_y = gauss_params[:, 3]
    theta = gauss_params[:, 2] #np.zeros(gauss_params.shape[0]) #gauss_params[:, 1]
    # x0 = gauss_params[:, 2]
    # y0 = sta.shape[2] - gauss_params[:, 3]
    # sigma_x = gauss_params[:, 4]
    # sigma_y = gauss_params[:, 5]
    # theta = gauss_params[:, 1]

    # Estimate the RF area from the Gaussian parameters.
    rf_area = np.pi * sigma_x * sigma_y
    
    # Get the total spike count.
    total_spikes = np.sum(spike_count[:,(0,2)], axis=1)

    # filepath = os.path.join(args.data_path, args.file_name, args.file_name + '.mat')
    filepath = os.path.join(args.data_path, args.file_name + '.mat')
    
    wr = ParamsWriter(filepath=filepath.replace('.mat','.params'), cluster_id=cluster_id)
    wr.write(timecourse_matrix[:,::-1,:], acf, total_spikes, x0, y0, sigma_x, sigma_y, theta, isi_binning=0.5)
    print("Wrote parameters file.")

    # Test the STAWriter class.
    sw = STAWriter(filepath=filepath.replace('.mat','.sta'))
    sw.write(sta=sta, cluster_id=cluster_id, stixel_size=micronsPerStixel, frame_refresh=1000/60)
    print("Wrote STA file.")

    width = float(sta.shape[3])
    height = float(sta.shape[2])
    if int(experiment_name[:8]) < 20230926:
        pixelsPerStixel = int(800.0 / width)
        mu_per_pixel = 3.8
    else:
        pixelsPerStixel = int(912.0 * 2.0 / width)
        if experiment_name[-1] == 'H':
            mu_per_pixel = 3.37
        else:
            mu_per_pixel = 3.24
    
    if array_id == 1501:
        micronsPerStixel = 20.0
    else:
        micronsPerStixel = 30.0
    refreshPeriod = 1000.0/60.0 

    runtime_movie_params = vl.RunTimeMovieParamsReader(pixelsPerStixelX = pixelsPerStixel,
                pixelsPerStixelY = pixelsPerStixel,
                width = width,
                height = height, #float(sta.shape[2]),
                micronsPerStixelX = micronsPerStixel,
                micronsPerStixelY = micronsPerStixel,
                xOffset = 0.0,
                yOffset = 0.0,
                interval = int(1), # same as stride, I think 
                monitorFrequency = 60.0,
                framesPerTTL = 1,
                refreshPeriod = refreshPeriod, 
                nFramesRequired = -1, # No idea what this means..
                droppedFrames = []) 
    
    with vw.GlobalsFileWriter(args.data_path, args.file_name) as gfw:
        gfw.write_simplified_litke_array_globals_file(array_id & 0xFFF, # FIXME get rid of this after we figure out what happened with 120um
                                                        0,
                                                        0,
                                                        'Kilosort converted',
                                                        '',
                                                        0,
                                                        a.vcd.n_samples)
        gfw.write_run_time_movie_params(runtime_movie_params)

    # Get the clusters.
    cell_clusters = cluster_compact(timecourse_matrix, acf, rf_area, total_spikes, variance_or_components=0.85,
            refractory_threshold=0.2, snr_threshold=2.0, isi_binning=0.5, count_threshold=300)
    category_list = get_cluster_labels(cell_clusters, cluster_id)
    # write_cluster_labels(os.path.join(args.data_path, args.file_name, 'online_clusters.txt'), category_list)
    write_cluster_labels(os.path.join(args.data_path, 'online_clusters.txt'), category_list)
    print("Wrote automated clustering")

    # Save a .mat file.
    # mdic = {'sta': sta, 'timecourse_matrix':timecourse_matrix, 'acf':acf, 'rf_area': rf_area, 'total_spikes':total_spikes, 'spatial_maps':spatial_maps, 'gauss_params':gauss_params}
    # hdf5storage.savemat( filepath, mdic, format=7.3, matlab_compatible=True, compress=False )
    