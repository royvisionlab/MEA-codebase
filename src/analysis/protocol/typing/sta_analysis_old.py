
import numpy as np
# import matplotlib.pyplot as plt
import os, sys
import platform
import cv2
import torch
from scipy.linalg import blas
from progressbar import *
from tqdm import tqdm
from scipy.io import savemat
import pickle
import hdf5storage
import argparse
from symphony_data import Metadata, Dataset, Stimulus
import lnp
import gc
from vision_utils import ParamsWriter, STAWriter, find_refractory_violations, find_bad_cells
import visionloader as vl
import visionwriter as vw
import config as cfg
import cv2
from typing import Tuple

sys.path.append(os.path.join(os.path.dirname(__file__), '../../../utilities'))
from yass_output_to_neurons import get_litke_triggers

def compute_rf_contours(significance_map: np.ndarray, max_contour_values: int=30) -> Tuple[np.ndarray, np.ndarray]:
    # Get the contour data.
    simple_contour_xy = np.zeros((max_contour_values,2))
    try:
        if np.sum(significance_map) > 2:
            _, thresh = cv2.threshold((255*significance_map).astype(np.uint8),0,255,0)
            cntrs = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
            cntrs = cntrs[0] if len(cntrs) == 2 else cntrs[1]
            c_xy = cntrs[0][:,0,:] # contour x,y coordinates
            simple_contour_area = cv2.contourArea(cntrs[0])
            num_values = c_xy.shape[0]
            for ii in range(max_contour_values):
                if ii < num_values:
                    simple_contour_xy[ii,:] = c_xy[ii,:]
                else:
                    simple_contour_xy[ii,:] = c_xy[-1,:]
            simple_contour_xy[-1,:] = c_xy[0,:]
    except:
        simple_contour_area = 0.0
    return simple_contour_xy, simple_contour_area

def fit_rf_hull_and_contour(space_map: np.ndarray, max_contour_values: int=30):
    rf_map = space_map.copy()
    for ii in range(rf_map.shape[2]):
        if np.nanmax(-rf_map[:,:,ii]) > np.nanmax(rf_map[:,:,ii]):
            rf_map[:,:,ii] = -rf_map[:,:,ii]
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        if np.ndim(rf_map) == 3:
            rf_map = np.nanmean(rf_map[:,:,(0,2)], axis=2)
        # Remove values less than zero.
        rf_map[rf_map < 0.0] = 0.0
        rf_map = rf_map / np.nanmax(rf_map)
        rf_map = np.nan_to_num(rf_map, nan=0.0, posinf=0.0, neginf=0.0)
    if np.all(rf_map == 0.0):
        return np.zeros((max_contour_values,2)), 0.0
    # Find the optimal threshold for the RF map.
    threshold = lnp.find_rf_threshold(rf_map)
    # Apply the threshold to the RF map.
    rf_map[rf_map < threshold] = 0.0
    rf_map[rf_map >= threshold] = 1.0
    try:
        simple_contour_xy, simple_contour_area = compute_rf_contours(rf_map, max_contour_values=max_contour_values)
    except Exception as error:
        print('An error occurred while computing the RF contour:', type(error).__name__, '-', error)
        simple_contour_xy = np.zeros((max_contour_values,2))
        simple_contour_area = 0.0
    return simple_contour_xy, simple_contour_area


class STAAnalysis(object):
    def __init__(self, 
            experimentName: str,
            sort_algorithm: str='yass', 
            protocol_id: str='FastNoise',
            frame_rate: float=60.31807657):
        self.experimentName = experimentName
        self.M = Metadata(experimentName)
        self.D = Dataset(experimentName)
        self.S = Stimulus()
        self.sort_algorithm = sort_algorithm
        self.protocol_id = protocol_id
        self.protocol = None
        self.frame_rate = frame_rate

    def get_protocol_from_file_name(self, file_name):
        self.protocol = self.M.search_data_file(self.protocol_id, file_name=file_name)
        return None
    
    def get_protocol_from_group(self, groupStr):
        self.protocol = self.M.searchProtocol(self.protocol_id, groupStr)
        return None

    def get_stimulus_frames(self, parameters, frame_times):
        if (self.protocol_id == 'FastNoise') or (self.protocol_id == 'AdaptNoiseColorSteps'):
            if 'gaussianFilter' in parameters:
                gaussianFilter = bool(parameters['gaussianFilter'])
                filterSdStixels = parameters['filterSdStixels']
            else:
                gaussianFilter = False
                filterSdStixels = 1.0
            frames = self.S.getFastNoiseFrames(
                int(parameters['numXStixels']),
                int(parameters['numYStixels']), 
                int(parameters['numXChecks']),
                int(parameters['numYChecks']),
                parameters['chromaticClass'],
                len(frame_times)-1,
                int(parameters['stepsPerStixel']), 
                int(parameters['seed']),
                int(parameters['frameDwell']),
                gaussianFilter=gaussianFilter,
                filterSdStixels=filterSdStixels)
        elif (self.protocol_id == 'SpatialNoise') or (self.protocol_id == 'AdaptNoiseSpatial'):
            if 'gaussianFilter' in parameters:
                gaussianFilter = bool(parameters['gaussianFilter'])
                filterSdStixels = parameters['filterSdStixels']
            else:
                gaussianFilter = False
                filterSdStixels = 1.0
            f_times = frame_times - frame_times[0]
            total_frames = len(frame_times)-1 #len(np.where(np.logical_and((f_times > parameters['preTime']),(f_times <= parameters['preTime']+parameters['stimTime'])))[0])
            unique_frames = len(np.where(np.logical_and((f_times > parameters['preTime']),(f_times <= parameters['preTime']+parameters['uniqueTime'])))[0])
            frames = self.S.get_spatial_noise_frames(
                int(parameters['numXStixels']),
                int(parameters['numYStixels']), 
                int(parameters['numXChecks']),
                int(parameters['numYChecks']),
                parameters['chromaticClass'],
                unique_frames,
                total_frames - unique_frames,
                int(parameters['stepsPerStixel']), 
                int(parameters['seed']),
                int(parameters['frameDwell']),
                gaussianFilter=gaussianFilter,
                filterSdStixels=filterSdStixels)

        elif (self.protocol_id == 'PinkNoise'):
            frames = self.S.getPinkNoiseFrames( 
                parameters['numXChecks'], 
                parameters['numYChecks'],
                len(frame_times)-1,
                parameters['rmsContrast'], 
                parameters['spatialAmplitude'],
                parameters['temporalAmplitude'],
                parameters['chromaticClass'],
                int(parameters['seed']))

        elif (self.protocol_id == 'SparseNoise'):
            frames = self.S.get_sparse_noise_frames(
                int(parameters['numXStixels']),
                int(parameters['numYStixels']), 
                int(parameters['numXChecks']),
                int(parameters['numYChecks']),
                parameters['chromaticClass'],
                len(frame_times)-1,
                int(parameters['stepsPerStixel']), 
                int(parameters['seed']),
                int(parameters['frameDwell']),
                parameters['pixelDensity'])

        return frames

    def downsample_data(self, frames, binned_spikes, down_factor):
        n_frames = np.floor(binned_spikes.shape[0] / down_factor).astype(int)
        tmp_frames = np.zeros((n_frames,frames.shape[1],frames.shape[2],frames.shape[3]))
        tmp_spikes = np.zeros((n_frames, binned_spikes.shape[1]))
        for frame in range(n_frames):
            fr_idx = frame*down_factor+ np.arange(down_factor)
            tmp_frames[frame,...] = np.mean(frames[fr_idx,...],axis=0)
            tmp_spikes[frame,:] = np.sum(binned_spikes[fr_idx,:],axis=0)
        binned_spikes = tmp_spikes
        frames = tmp_frames
        return frames, binned_spikes

    def normalize_sta(self, sta_tensor):
        """
        Normalizes STA according to some constants. Puts it in 0-1, with mean .5.
        Makes a copy, does not mutate.
        Parameters:
            sta_tensor: STA tensor of size cells, t, y, x, color
        Returns:
            normalized STA tensor of same size as input
        """
        sta_tensor_norm = sta_tensor.copy()
        sta_tensor_norm /= np.amax(np.abs(sta_tensor_norm),
                            axis=(1,2,3,4)).reshape((-1,1,1,1,1))
        sta_tensor_norm *= 0.5
        sta_tensor_norm += 0.5

        return sta_tensor_norm

    def compute_sta_blas(stimulus,binned_spikes,stride,depth=30,norm=True):
        """
        Computes STA by convolving the spiketrain with the stimulus. This 
        implementation relies on the BLAS interfaces that numpy wraps around.

        Parameters:
            stimulus: stimulus tensor of size t, y, x, color (8 bit ints).
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

        # Check that the FORTRAN type is reasonable.
        blas_obj = blas.get_blas_funcs('gemm',(stimulus,binned_spikes))

        if blas_obj.typecode != 'd':
            assert False, "The stimulus or spikes array were not of expected type."

        # Initialize STA, the offset range, and compute.
        offset_cnt = 0
        sta_tensor = np.zeros((num_cells,depth,height,width,n_phosphors))
        offset_range = np.arange(0,(depth * stride),stride)

        for offset in offset_range:

            # Ensure the spike and stimulus vector keep their byte order.
            spikes = np.array(binned_spikes[offset:,:],order="F")

            assert spikes.flags.fortran,\
                "Could not get spikes vector in FORTRAN order for BLAS routines."
            assert stimulus[:,:n_frames-offset].flags.fortran,\
                "Could not get stimulus in FORTRAN order for BLAS routines."

            # Compute the STA and write to the tensor.
            sta = blas.dgemm(alpha=1.0, a=stimulus[:,:n_frames-offset],
                            b=spikes)

            # Index appropriately to avoid nans and normalize by number of spikes.
            nonzero_inds = np.nonzero(np.sum(spikes,0))[0]
            sta[:,nonzero_inds] /= np.sum(spikes,0)[nonzero_inds]

            del spikes

            sta_tensor[:,offset_cnt,:,:,:] = np.reshape(sta.T,
                                                    (num_cells,height,
                                                    width,n_phosphors))

            del sta
            offset_cnt +=1

        return sta_tensor

    def compute_sta_torch(self, stimulus, binned_spikes, stride, depth=30, norm=False):
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

        # Set up torch stuff.
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu') #torch.device('cpu')

        # Vectorize stimulus.
        stimulus = np.reshape(stimulus,(n_frames,
                                        height * width * n_phosphors)).astype(np.float32).T
        
        # Initialize STA, the offset range, and compute.
        offset_cnt = 0
        sta_tensor = np.zeros((num_cells,depth,height,width,n_phosphors))
        offset_range = np.arange(0,(depth * stride),stride)

        # Initialize progress bar.
        # widgets = ['Computing STAs: ', Percentage(), ' ', Bar(marker='*',
        #         left='[',right=']'),' ', ETA()]
        # pbar = ProgressBar(widgets=widgets, maxval=depth)
        # pbar.start()

        
        
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
            # pbar.update(offset_cnt)
            offset_cnt +=1

        # pbar.finish()

        if device.type == 'cuda':
            torch.cuda.empty_cache()
        
        # Clear some memory
        del sta
        gc.collect()

        if not norm:
            return sta_tensor

        return self.normalize_sta(sta_tensor)

    def save_sta_mat(self, filepath, mdic):
        # savemat(filepath, mdic)
        hdf5storage.savemat( filepath, mdic, format=7.3, matlab_compatible=True, compress=False )
        return None

    def load_sta_mat(self, filepath):
        mdic = hdf5storage.loadmat(filepath)
        return mdic['sta'], mdic['cluster_id'], mdic['spike_count']

    def save_sta_pickle(self, filepath, sta, cluster_id, spike_count):
        d = {'sta':sta, 'cluster_id':cluster_id, 'spike_count': spike_count}
        with open(filepath, 'wb') as file:
            pickle.dump(d, file)
        return None

    def load_sta_pickle(self, filepath):
        d = pickle.load( open(filepath, 'rb') )
        return d['sta'], d['cluster_id']

    def get_autocorrelation(self, bin_edges=np.linspace(0,300,601)):
        acf_dict = dict() # Autocorrelation dictionary

        # protocol = self.M.searchProtocol(searchStr, groupStr)
        protocol = self.protocol

        # Get the number of groups.
        nGroups = len(protocol['group'])
        # Loop through the EpochGroups
        for groupCount, group in enumerate(protocol['group']):

            nBlocks = len(group['block'])
            # Grab the epoch blocks
            for blockCount, block in enumerate(group['block']):
                print('Processing group ' + str(groupCount+1) + ' of ' + str(nGroups) + ' and block ' + str(blockCount+1) + ' of ' + str(nBlocks))

                dataFile = self.M.get_data_file_from_block(block)

                # Try to load a data file.
                filepath = os.path.join(self.M.sortpath, self.experimentName, dataFile, self.sort_algorithm)

                self.D.load_vision_table(filepath, dataFile)

                # Compute the autocorrelation function.
                acf_tmp = self.D.get_interspike_intervals(bin_edges=bin_edges)

                for cell in self.D.cell_ids:
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

    def compute_nonlinearity(self,
                    sta,
                    cluster_id,
                    idx_good: np.ndarray=None,
                    stride: int=2,
                    num_bins: int=25,
                    device_type: str='gpu'):

        # protocol = self.M.searchProtocol(protocol_id, group_str)
        protocol = self.protocol
        generator = np.zeros((len(cluster_id),1))
        response = np.zeros((len(cluster_id),1))

        # Get the number of groups.
        nGroups = len(protocol['group'])
        # Loop through the EpochGroups
        for groupCount, group in enumerate(protocol['group']):

            nBlocks = len(group['block'])
            # Grab the epoch blocks
            for blockCount, block in enumerate(group['block']):
                print('Processing group ' + str(groupCount+1) + ' of ' + str(nGroups) + ' and block ' + str(blockCount+1) + ' of ' + str(nBlocks))

                dataFile = self.M.get_data_file_from_block(block)

                # Try to load a data file.
                filepath = os.path.join(self.M.sortpath, self.experimentName, dataFile, self.sort_algorithm)

                self.D.load_vision_table(filepath, dataFile)

                epoch_frames = self.M.get_frame_times_ms( block )
                # Loop through each epoch.
                for count, epoch in enumerate(tqdm(block['epoch'], desc='Computing generator for block ' + str(blockCount+1) + ' of ' + str(nBlocks))):
                    # print('Analyzing epoch ' + str(count))
                    if count < len(epoch_frames):
                        parameters = epoch['parameters']

                        if self.protocol_id == 'FastNoise':
                            if 'stepsPerStixel' not in parameters:
                                parameters['stepsPerStixel'] = np.round(parameters['stixelSize'] / parameters['gridSize']).astype(int)
                        
                        # Need to get rid of pre/post frames.
                        preFrames = round(parameters['preTime']*1e-3 * 60)
                        stimFrames = round(parameters['stimTime']*1e-3 * 60)

                        frame_times = epoch_frames[count]

                        # Get the frames when the stimulus was on.
                        idx = (frame_times-frame_times[0]).searchsorted((parameters['preTime']+parameters['stimTime']),'right')
                        frame_times = frame_times[preFrames : idx+1]

                        frames = self.get_stimulus_frames(parameters, frame_times)

                        # Check for frame drops.
                        frame_times, transition_frames = self.S.check_frame_times(frame_times)

                        # Fix any frame drops.
                        if np.any(transition_frames > 1):
                            frames = self.S.fix_frames(frames, transition_frames)
                        
                        if stride > 1:
                            frames = self.S.upsample_frames(frames, stride)

                        # Get the binned spike count.
                        gen_idx = np.isin(cluster_id, self.D.cell_ids)
                        binned_spikes = self.D.get_binned_spikes(frame_times, stride)

                        if (self.protocol_id == 'SparseNoise'):
                            frame_dwell = parameters['frameDwell']
                            frames, binned_spikes = self.downsample_data(frames=frames, binned_spikes=binned_spikes, down_factor=frame_dwell)

                        binned_spikes = binned_spikes.T
                        response_tmp = np.zeros((len(cluster_id),binned_spikes.shape[1]))
                        response_tmp[gen_idx,:] = binned_spikes
                        response = np.append(response, response_tmp, axis = 1)

                        # Compute the generator signal. Stride of convolution must be 1 as stride accounted for in frames and binned_spikes, just like in compute_sta.
                        if device_type == 'gpu':
                        # Figure out how many chunks to use for computing the generator.
                            num_chunks = np.ceil(np.size(sta)*4.0 / (5000000000.0-np.size(frames)*4)).astype(int) #np.ceil(np.size(sta)*4.0 / (10000000000.0-np.size(frames)*4)).astype(int)
                            cells_per_chunk = np.floor(sta.shape[0]/num_chunks).astype(int)
                            gen_tmp = np.zeros((sta.shape[0],binned_spikes.shape[1]))
                            for chunk in range(num_chunks):
                                cell_idx = range(chunk*cells_per_chunk, np.min([sta.shape[0],(chunk+1)*cells_per_chunk]).astype(int))
                                gen_tmp[cell_idx,:] = lnp.compute_generator_signal_torch(frames, sta[cell_idx,...], stride=1, alpha=1/10000, device_type=device_type)
                        else:
                            gen_tmp = lnp.compute_generator_signal_torch(frames, sta, stride=1, alpha=1/10000, device_type=device_type)
                        generator = np.append(generator, gen_tmp, axis = 1)
                        del gen_tmp, binned_spikes, response_tmp

        # Bin to get the nonlinearity
        if idx_good is None:
            x_bin, y_bin = lnp.bin_nonlinearities(generator=generator[:,1:], response=response[:,1:], num_bins=num_bins)
        else:
            x_bin, y_bin = lnp.bin_nonlinearities(generator=generator[:,1:], response=response[idx_good,1:], num_bins=num_bins)
        
        return x_bin, y_bin

    def compute_sta_old(self, 
                    stride: int=2,
                    depth: int=61):

        # protocol = self.M.searchProtocol(protocol_id, groupStr)
        protocol = self.protocol
        stas = dict()
        spike_counts = dict()
        spike_counts_blue = dict() # Spike counts for blue frames.

        total_count = 0

        # Get the number of groups.
        nGroups = len(protocol['group'])
        # Loop through the EpochGroups
        for groupCount, group in enumerate(protocol['group']):

            nBlocks = len(group['block'])
            # Grab the epoch blocks
            for blockCount, block in enumerate(group['block']):
                print('Processing group ' + str(groupCount+1) + ' of ' + str(nGroups) + ' and block ' + str(blockCount+1) + ' of ' + str(nBlocks))

                dataFile = self.M.get_data_file_from_block(block)
                print(dataFile)

                # if (dataFile == 'data027') or (dataFile == 'data028') or (dataFile == 'data032') or (dataFile == 'data033'): #(dataFile == 'data005') or (dataFile == 'data006') or (dataFile == 'data010') or (dataFile == 'data011'):
                #     continue

                # Try to load a data file.
                filepath = os.path.join(self.M.sortpath, self.experimentName, dataFile, self.sort_algorithm)

                self.D.load_vision_table(filepath, dataFile)

                epoch_frames = self.M.get_frame_times_ms( block )
                # Loop through each epoch.
                for count, epoch in enumerate(tqdm(block['epoch'], desc='Computing STA for block ' + str(blockCount+1) + ' of ' + str(nBlocks))):
                # for count, epoch in enumerate(block['epoch']):
                    if count < len(epoch_frames):
                        parameters = epoch['parameters']

                        if self.protocol_id == 'FastNoise':
                            if 'stepsPerStixel' not in parameters:
                                parameters['stepsPerStixel'] = np.round(parameters['stixelSize'] / parameters['gridSize']).astype(int)
                        
                        # Need to get rid of pre/post frames.
                        pre_frames = round(parameters['preTime']*1e-3 * 60)
                        stimFrames = round(parameters['stimTime']*1e-3 * 60)

                        frame_times = epoch_frames[count]

                        # Get the frames when the stimulus was on.
                        idx = (frame_times-frame_times[0]).searchsorted((parameters['preTime'] + parameters['stimTime']), 'right')
                        frame_times = frame_times[pre_frames : idx+1]

                        frames = self.get_stimulus_frames(parameters, frame_times)

                        # Downsample the frames if they are too big.
                        if frames.shape[2] > 200:
                            f_copy = frames.copy()
                            frames = np.zeros((f_copy.shape[0], 100, 75, f_copy.shape[3]))
                            for t_val in range(f_copy.shape[0]):
                                for c_val in range(f_copy.shape[3]):
                                    frames[t_val,:,:,c_val] = cv2.resize(f_copy[t_val,:,:,c_val], dsize=(75, 100), interpolation=cv2.INTER_AREA)

                        # Check for frame drops.
                        frame_times, transition_frames = self.S.check_frame_times(frame_times)

                        # print(frame_times.shape)
                        # print(transition_frames.shape)
                        # print(frames.shape)
                        # print(np.argwhere(transition_frames > 1))

                        # Fix any frame drops.
                        if np.any(transition_frames > 1):
                            frames = self.S.fix_frames(frames, transition_frames)
                        
                        if stride > 1:
                            frames = self.S.upsample_frames(frames, stride)

                        # Get the binned spike count.
                        binned_spikes = self.D.get_binned_spikes(frame_times, stride)
                        binned_spikes = binned_spikes.astype(float)

                        if frames.shape[2] > 200:
                            num_chunks = np.ceil((float(frames.shape[2])/100.0)**2 + 2).astype(int)
                            cells_per_chunk = np.ceil(float(binned_spikes.shape[1]) / float(num_chunks)).astype(int)
                            pix_per_chunk = np.ceil(float(frames.shape[2]) / float(num_chunks)).astype(int)
                        else:
                            num_chunks = 1
                            cells_per_chunk = binned_spikes.shape[0]

                        if (self.protocol_id == 'SparseNoise'):
                            frame_dwell = parameters['frameDwell']
                            frames, binned_spikes = self.downsample_data(frames=frames, binned_spikes=binned_spikes, down_factor=frame_dwell)

                        print('num chunks:' + str(num_chunks))
                        # Compute the STA. Use BLAS for really large datasets.
                        if num_chunks == 1:
                            sta_tensor = self.compute_sta_torch(frames, binned_spikes, 1, depth)
                        else:
                            sta_tensor = np.zeros((binned_spikes.shape[1],depth,frames.shape[1],frames.shape[2],3), dtype=np.float32)
                            for chunk in range(num_chunks):
                                cell_idx = np.arange(chunk*cells_per_chunk, np.min([binned_spikes.shape[1],(chunk+1)*cells_per_chunk]).astype(int))
                                pix_idx = np.arange(chunk*pix_per_chunk, np.min([frames.shape[2],(chunk+1)*pix_per_chunk]).astype(int))
                                sta_tensor[:,:,:,pix_idx,:] = self.compute_sta_torch(frames[:,:,pix_idx,:], binned_spikes, 1, depth)
                                # sta_tensor[cell_idx,:,:,:,:] = self.compute_sta_torch(frames, binned_spikes[:,cell_idx], 1, depth)
                        # if (frames.shape[3] > 150) and (depth > 5):
                        #     sta_tensor = lnp.compute_sta_blas(frames,binned_spikes,stride=1,depth=depth,norm=False)
                        # else:
                        #     sta_tensor = self.compute_sta_torch(frames, binned_spikes, 1, depth)
                        
                        # Get the total spike count.
                        sp_count = np.sum(binned_spikes, axis=0)

                        for cell_count, cell in enumerate(self.D.cell_ids):
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
                        
                        total_count += 1

                        del sta_tensor, binned_spikes

        # Sort the STA dictionary by cluster.
        stas = dict(sorted(stas.items()))

        # Get the STA by dividing by the spike count. 
        for cell in stas:
            if (parameters['chromaticClass'] != 'achromatic'):
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
        
        return sta, np.array(cluster_id), spike_count

    def compute_sta(self, 
                    stride: int=2,
                    depth: int=61):

        # protocol = self.M.searchProtocol(protocol_id, groupStr)
        protocol = self.protocol
        stas = dict()
        spike_counts = dict()
        spike_counts_blue = dict() # Spike counts for blue frames.

        total_count = 0

        # Get the number of groups.
        nGroups = len(protocol['group'])
        # Loop through the EpochGroups
        for groupCount, group in enumerate(protocol['group']):

            nBlocks = len(group['block'])
            # Grab the epoch blocks
            for blockCount, block in enumerate(group['block']):
                print('Processing group ' + str(groupCount+1) + ' of ' + str(nGroups) + ' and block ' + str(blockCount+1) + ' of ' + str(nBlocks))

                dataFile = self.M.get_data_file_from_block(block)
                print(dataFile)

                # if (dataFile == 'data027') or (dataFile == 'data028') or (dataFile == 'data032') or (dataFile == 'data033'): #(dataFile == 'data005') or (dataFile == 'data006') or (dataFile == 'data010') or (dataFile == 'data011'):
                #     continue

                # Try to load a data file.
                filepath = os.path.join(self.M.sortpath, self.experimentName, dataFile, self.sort_algorithm)

                self.D.load_vision_table(filepath, dataFile)

                epoch_frames = self.M.get_frame_times_ms( block )
                # Loop through each epoch.
                for count, epoch in enumerate(tqdm(block['epoch'], desc='Computing STA for block ' + str(blockCount+1) + ' of ' + str(nBlocks))):
                # for count, epoch in enumerate(block['epoch']):
                    if count < len(epoch_frames):
                        parameters = epoch['parameters']

                        if self.protocol_id == 'FastNoise':
                            if 'stepsPerStixel' not in parameters:
                                parameters['stepsPerStixel'] = np.round(parameters['stixelSize'] / parameters['gridSize']).astype(int)
                        
                        # Need to get rid of pre/post frames.
                        pre_frames = round(parameters['preTime']*1e-3 * 60)

                        frame_times = epoch_frames[count]

                        # Get the frames when the stimulus was on.
                        idx = (frame_times-frame_times[0]).searchsorted((parameters['preTime'] + parameters['stimTime']), 'right')
                        frame_times = frame_times[pre_frames : idx+1]

                        frames = self.get_stimulus_frames(parameters, frame_times)

                        # Check for frame drops.
                        frame_times, transition_frames = self.S.check_frame_times(frame_times)

                        # Fix any frame drops.
                        if np.any(transition_frames > 1):
                            frames = self.S.fix_frames(frames, transition_frames)
                        
                        if stride > 1:
                            frames = self.S.upsample_frames(frames, stride)

                        # Get the binned spike count.
                        binned_spikes = self.D.get_binned_spikes(frame_times, stride)
                        binned_spikes = binned_spikes.astype(float)

                        if (self.protocol_id == 'SparseNoise'):
                            frame_dwell = parameters['frameDwell']
                            frames, binned_spikes = self.downsample_data(frames=frames, binned_spikes=binned_spikes, down_factor=frame_dwell)
                        # Compute the STA. Use BLAS for really large datasets.
                        sta_tensor = self.compute_sta_torch(frames, binned_spikes, 1, depth)
                        
                        # Get the total spike count.
                        sp_count = np.sum(binned_spikes, axis=0)

                        for cell_count, cell in enumerate(self.D.cell_ids):
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
                        
                        total_count += 1

                        del sta_tensor, binned_spikes

        # Sort the STA dictionary by cluster.
        stas = dict(sorted(stas.items()))
        
        return stas

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Compute the spike-triggered average from noise data.')
    parser.add_argument('experimentName', type=str, help='folder containing Kilosort outputs, must already exist')
    parser.add_argument('-a','--algorithm', default='kilosort2', type=str, help='Sorting algorithm used (yass or kilosort2)')
    parser.add_argument('-p','--protocol', default='FastNoise', type=str, help='Name of Symphony stimulus protocol (default: FastNoise)')
    parser.add_argument('-g','--group', default=None, type=str, help='search string within epoch group label (default: None)')
    parser.add_argument('-f','--file', nargs='+', default=None, type=str, help='name of single data file to analyze i.e. data003 (default: None)')
    parser.add_argument('-s', '--stride', default=2, type=int, help='path to write Vision files, must already exist')
    parser.add_argument('-d','--depth', default=61, type=int, help='name of dataset that we want to write, i.e. data000')
    parser.add_argument('-n','--nonlinearity', default=False, action="store_true", help='compute the nonlinearity?')
    parser.add_argument('-c','--chunk', default=None, type=str, help='chunk name')

    args = parser.parse_args()

    experiment_name = args.experimentName
    sort_algorithm = args.algorithm
    protocol_id = args.protocol

    a = STAAnalysis(experiment_name, sort_algorithm=sort_algorithm, protocol_id=protocol_id)
    if args.file != None:
        file_names = args.file
        if type(args.file) is str:
            file_names = [args.file]
        a.get_protocol_from_file_name(file_names)
        if type(args.file) is list:
            label = '_'.join(args.file)
        else:
            label = args.file
    else:
        a.get_protocol_from_group(args.group)
        label = ''

    if args.protocol == 'SparseNoise':
        stride = 1
        depth = 3
    else:
        stride = args.stride
        depth = args.depth

    # if (platform.node() == 'mike'):
    #     file_dir = '/home/mike/ftp/files/'  
    # else:
    #     file_dir = '/gscratch/retina/data/sorted/' + args.experimentName + '/'

    SORT_PATH, _ = cfg.get_data_paths()
    file_dir = SORT_PATH + args.experimentName + '/'
    file_str = args.experimentName + '_' + args.algorithm + '_' + label #SORT_PATH + args.experimentName + '/' + args.experimentName + '_' + args.algorithm + '_' + label
    # file_str = args.experimentName + '_' + args.algorithm + '_' + label
    
    # Compute the ACF
    acf, isi, _ = a.get_autocorrelation()

    # Find the cells that obey refractoriness.
    _, idx_good = find_refractory_violations(acf, refractory_threshold = 0.2)

    # Compute the STA.
    sta, cluster_id, spike_count = a.compute_sta(stride, depth)

    nonlinearity_bins = 25
    width = float(sta.shape[3])
    height = float(sta.shape[2])
    if int(experiment_name[:8]) < 20230926:
        pixelsPerStixel = int(800.0 / width)
        mu_per_pixel = 3.8
    else:
        pixelsPerStixel = int(912.0 * 2.0 / width)
        mu_per_pixel = 2.43
    micronsPerStixel = np.round(pixelsPerStixel * mu_per_pixel).astype(float)
    refreshPeriod = 1000.0/a.frame_rate/float(stride) # STA refresh period in msec.

    # Compute the spatiotemporal RF parameters.
    if args.protocol == 'SparseNoise':
        # Compute the nonlinearity.
        x_bin, y_bin = a.compute_nonlinearity(sta, cluster_id=cluster_id, stride=stride, num_bins=nonlinearity_bins)
        mdic = {'sta': sta, 'cluster_id': cluster_id, 'spike_count': spike_count, 'acf':acf, 'isi':isi, 'x_bin': x_bin, 'y_bin': y_bin }
        savemat(file_dir + file_str + '_SparseNoise.mat', mdic)
    else:

        if args.chunk is None:
            out_file_path = file_dir + file_str
            globals_path = file_dir
            globals_name = file_str
        else:
            out_file_path = file_dir + args.chunk + '/' + args.algorithm + '/' + args.algorithm
            globals_path = file_dir + args.chunk + '/' + args.algorithm + '/'
            globals_name = args.algorithm

        # Compute the spatiotemporal RF parameters.
        timecourse_matrix, spatial_maps, significance_maps, hull_parameters, hull_area, hull_vertices, rank1_r2 = lnp.compute_spatiotemporal_maps(sta)

        # Compute the RF contour.
        max_contour_values=30
        simple_contour_xy = np.zeros((spatial_maps.shape[0], max_contour_values, 2))
        simple_contour_area = np.zeros((spatial_maps.shape[0]))
        for c_idx in range(spatial_maps.shape[0]):
            try:
                simple_contour_xy[c_idx,...], simple_contour_area[c_idx] = fit_rf_hull_and_contour(space_map=spatial_maps[c_idx,...], max_contour_values=max_contour_values)
            except:
                pass
        # Flip the y coordinates.
        simple_contour_xy[:, :, 1] = sta.shape[2] - simple_contour_xy[:, :, 1]
        # hull_vertices[:, :, 1] = sta.shape[2] - hull_vertices[:, :, 1]

        # Fix the zero values in the hull vertices.
        for i in range(hull_vertices.shape[0]):
            zero_idx = np.where(hull_vertices[i,:,0] == 0)[0]
            if len(zero_idx) > 0:
                hull_vertices[i,zero_idx,0] = hull_vertices[i,0,0]
                hull_vertices[i,zero_idx,1] = hull_vertices[i,0,1]

        hull_vertices[:, :, 1] = sta.shape[2] - hull_vertices[:, :, 1]
        
        # Normalize the spatial maps.
        spatial_maps = spatial_maps / np.max(np.abs(spatial_maps), axis=(1,2,3))[:,None,None,None]
        spatial_maps = np.nan_to_num(spatial_maps, nan=0.0, posinf=0.0, neginf=0.0) # Replace NaNs with zeros.

        # Set the last frame of the STA to the spatial map.
        last_frame = sta[:,-1,:,:,:]
        sta[:,-1,:,:,:] = spatial_maps

        # Unpack the hull parameters.
        x0 = hull_parameters[:, 0]
        y0 = sta.shape[2] - hull_parameters[:, 1]
        sigma_x = hull_parameters[:, 2]
        sigma_y = hull_parameters[:, 3]
        theta = hull_parameters[:, 4]

        # For binary stimuli.
        timecourse_variance = -np.sqrt(timecourse_matrix*timecourse_matrix)
        stv = 1 - (sta*sta)

        #timecourse_matrix, isi, spike_count, x0, y0, sigma_x, sigma_y, 
            #   theta: np.ndarray=None, isi_binning: float=0.5, 
            #   contour_xy: np.ndarray=None, contour_area: np.ndarray=None, simple_contour_xy: np.ndarray=None, 
            #   simple_contour_area: np.ndarray=None, timecourse_variance: np.ndarray=None, ei_parameters: np.ndarray=None
        with ParamsWriter(filepath=out_file_path + '.params', cluster_id=cluster_id) as wr:
            wr.write(timecourse_matrix=timecourse_matrix[:,::-1,:], isi=acf, spike_count=np.sum(spike_count[:,(0,2)],axis=1), x0=x0, y0=y0, 
                        sigma_x=sigma_x, sigma_y=sigma_y, theta=theta, isi_binning=0.5)
            # wr.write(timecourse_matrix=timecourse_matrix[:,::-1,:], isi=acf, spike_count=np.sum(spike_count[:,(0,2)],axis=1), x0=x0, y0=y0, 
            #             sigma_x=sigma_x, sigma_y=sigma_y, theta=theta, isi_binning=0.5, contour_xy=hull_vertices, contour_area=hull_area,
            #             simple_contour_xy=simple_contour_xy, simple_contour_area=simple_contour_area, timecourse_variance=timecourse_variance[:,::-1,:])

        runtime_movie_params = vl.RunTimeMovieParamsReader(pixelsPerStixelX = pixelsPerStixel,
                    pixelsPerStixelY = pixelsPerStixel,
                    width = width,
                    height = height, #float(sta.shape[2]),
                    micronsPerStixelX = micronsPerStixel,
                    micronsPerStixelY = micronsPerStixel,
                    xOffset = 0.0,
                    yOffset = 0.0,
                    interval = int(stride), # same as stride, I think 
                    monitorFrequency = a.frame_rate,
                    framesPerTTL = 1,
                    refreshPeriod = refreshPeriod, 
                    nFramesRequired = -1, # No idea what this means..
                    droppedFrames = []) 
        
        # Get the array id.
        data_path_tmp = SORT_PATH.replace('sorted/','raw/') + args.experimentName + '/' + file_names[0]
        if os.path.isdir(data_path_tmp):
            _, _, array_id, _ = get_litke_triggers(data_path_tmp, RW_BLOCKSIZE=2000000, TTL_THRESHOLD=-1000)
        else:
            array_id = 504
        
        with vw.GlobalsFileWriter(globals_path, globals_name) as gfw:
            gfw.write_simplified_litke_array_globals_file(array_id & 0xFFF, # FIXME get rid of this after we figure out what happened with 120um
                                                            0,
                                                            0,
                                                            'Kilosort converted',
                                                            '',
                                                            0,
                                                            a.D.vcd.n_samples)
            gfw.write_run_time_movie_params(runtime_movie_params) 
        
        if args.nonlinearity: # Compute the nonlinearity.
            x_bin, y_bin = a.compute_nonlinearity(sta, cluster_id=cluster_id, stride=stride, num_bins=nonlinearity_bins, device_type='gpu')
        else:
            x_bin = []
            y_bin = []

        # mdic = {'sta': sta, 'cluster_id': cluster_id, 'spike_count': spike_count, 'acf':acf, 'isi':isi,
        #         'timecourse_matrix': timecourse_matrix, 'spatial_maps': spatial_maps, 
        #         'significance_maps': significance_maps, 'hull_area': hull_area, 'hull_vertices': hull_vertices, 
        #         'hull_parameters': hull_parameters, 'rank1_r2': rank1_r2, 'x_bin': x_bin, 'y_bin': y_bin }

        mdic2 = {'cluster_id': cluster_id, 'spike_count': spike_count, 'acf':acf, 'isi':isi,
                'timecourse_matrix': timecourse_matrix, 'spatial_maps': spatial_maps, 
                'significance_maps': significance_maps, 'hull_area': hull_area, 'hull_vertices': hull_vertices, 
                'hull_parameters': hull_parameters, 'rank1_r2': rank1_r2, 'x_bin': x_bin, 'y_bin': y_bin,
                'simple_contour_xy': simple_contour_xy, 'simple_contour_area': simple_contour_area } 

        # Write the .STA file.
        #sta: np.ndarray=None, ste: np.ndarray=None, cluster_id: np.ndarray=None, stixel_size: float=30.0, frame_refresh
        # start_time = time.time()
        with STAWriter(filepath=out_file_path + '.sta') as sw:
            sw.write(sta=sta, ste=stv, cluster_id=cluster_id, stixel_size=micronsPerStixel, frame_refresh=refreshPeriod)
        # print('Writing STA took ' + str(time.time() - start_time) + ' seconds.')

        
#         import visionloader as vl
#         import time
# sta_by_cell_id = dict()
# for cell_count, cell in enumerate(cluster_id):
#     sta_by_cell_id[cell] = vl.STAContainer(stixel_size=30.0,refresh_time=refreshPeriod,sta_offset=0,
#         red=np.transpose(sta[cell_count,::-1,:,:,0],(2,1,0)),
#         red_error=np.transpose(stv[cell_count,::-1,:,:,0],(2,1,0)),
#         green=np.transpose(sta[cell_count,::-1,:,:,1],(2,1,0)),
#         green_error=np.transpose(stv[cell_count,::-1,:,:,1],(2,1,0)),
#         blue=np.transpose(sta[cell_count,::-1,:,:,2],(2,1,0)),
#         blue_error=np.transpose(stv[cell_count,::-1,:,:,2],(2,1,0)))
#         analysis_folder_path='/data/data/sorted/20231003C/chunk1/kilosort2'
#         dataset_name = 'kilosort2'
# start_time = time.time()
# with vw.STAWriter(analysis_folder_path=analysis_folder_path,
#             dataset_name=dataset_name,
#             width=int(height),
#             height=int(width),
#             depth=depth,
#             refresh_time=refreshPeriod,
#             sta_offset=0,
#             stixel_size=30.0) as f:
#     f.write_sta_by_cell_id(sta_by_cell_id)
# print('Writing STA took ' + str(time.time() - start_time) + ' seconds.')
        # sw = STAWriter(filepath=out_file_path + '.sta')
        # sw.write(sta, cluster_id, stixel_size=micronsPerStixel, frame_refresh=refreshPeriod)

        # Replace the last frame of the STA with the last frame of the original STA.
        # sta[:,-1,:,:,:] = last_frame
        
        # Save MAT files.
        # Save the shortened params file first.
        savemat(out_file_path + '_params.mat', mdic2)

        # a.save_sta_mat(out_file_path + '_sta.mat', mdic)
        