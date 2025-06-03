import os, sys
import numpy as np
import gc
import torch
import config as cfg
import argparse
from symphony_data import Dataset, Stimulus
from tqdm import tqdm
import hdf5storage
from lnp import find_rf_threshold, calculate_time_course_svd, compute_spatial_map, compute_spatiotemporal_maps, compute_generator_signal_torch, bin_nonlinearities
import warnings
import cv2
from typing import Tuple
from vision_utils import ParamsWriter, STAWriter
import visionloader as vl
import visionwriter as vw

sys.path.append(os.path.join(os.path.dirname(__file__), '../../../utilities'))
from yass_output_to_neurons import get_litke_triggers

def get_array_properties(SORT_PATH, experiment_name, file_name):
    num_samples = 0
    # Get the array id.
    for ii in range(len(file_name)):
        data_path_tmp = SORT_PATH.replace('sorted/','raw/') + experiment_name + '/' + file_name[ii]
        if os.path.isdir(data_path_tmp):
            _, _, array_id, n_samples = get_litke_triggers(data_path_tmp, RW_BLOCKSIZE=2000000, TTL_THRESHOLD=-1000)
            num_samples += n_samples
        else:
            array_id = 504
    return int(array_id), int(num_samples)

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
    threshold = find_rf_threshold(rf_map)
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

def write_globals_file(globals_path: str, globals_name: str, array_id: int, num_samples: int, sta_width: int, sta_height: int, micronsPerStixel: float, pixelsPerStixel: int, mean_frame_rate: float, stride: int):
    refreshPeriod = 1000.0/mean_frame_rate/float(stride) # STA refresh period in msec.
    runtime_movie_params = vl.RunTimeMovieParamsReader(pixelsPerStixelX = pixelsPerStixel,
                    pixelsPerStixelY = pixelsPerStixel,
                    width = sta_width,
                    height = sta_height, 
                    micronsPerStixelX = micronsPerStixel,
                    micronsPerStixelY = micronsPerStixel,
                    xOffset = 0.0,
                    yOffset = 0.0,
                    interval = int(stride), # same as stride, I think 
                    monitorFrequency = mean_frame_rate,
                    framesPerTTL = 1,
                    refreshPeriod = refreshPeriod, 
                    nFramesRequired = -1, # No idea what this means..
                    droppedFrames = []) 
    
    with vw.GlobalsFileWriter(globals_path, globals_name) as gfw:
        gfw.write_simplified_litke_array_globals_file(array_id & 0xFFF, # FIXME get rid of this after we figure out what happened with 120um
                                                        0,
                                                        0,
                                                        'Kilosort converted',
                                                        '',
                                                        0,
                                                        num_samples)
        gfw.write_run_time_movie_params(runtime_movie_params) 

def write_params_file(out_path: str, cluster_id: np.ndarray, timecourse_matrix: np.ndarray, isi: np.ndarray, spike_count: np.ndarray, 
            hull_parameters: np.ndarray, isi_binning: float=0.5, contour_xy: np.ndarray=None, contour_area: np.ndarray=None, simple_contour_xy: np.ndarray=None, 
            simple_contour_area: np.ndarray=None, timecourse_variance: np.ndarray=None, ei_parameters: np.ndarray=None, sta_size: np.ndarray=None):
    # Unpack the hull parameters.
    x0 = hull_parameters[:, 0]
    y0 = sta_size[1] - hull_parameters[:, 1]
    sigma_x = hull_parameters[:, 2]
    sigma_y = hull_parameters[:, 3]
    theta = hull_parameters[:, 4]
    if timecourse_variance is not None:
        timecourse_variance=timecourse_variance[:,::-1,:] # Reverse the timecourse variance.
    with ParamsWriter(filepath=out_path + '.params', cluster_id=cluster_id) as wr:
        wr.write(timecourse_matrix=timecourse_matrix[:,::-1,:], isi=isi, spike_count=spike_count, x0=x0, y0=y0, sigma_x=sigma_x, sigma_y=sigma_y, theta=theta, isi_binning=isi_binning,
             contour_xy=contour_xy, contour_area=contour_area, simple_contour_xy=simple_contour_xy, simple_contour_area=simple_contour_area, timecourse_variance=timecourse_variance,ei_parameters=ei_parameters)

def write_sta_file(out_path: str, sta: np.ndarray, ste: np.ndarray, cluster_id: np.ndarray, stixel_size: float=30.0):
    with STAWriter(filepath=out_path + '.sta') as wr:
        wr.write(sta=sta, ste=ste, cluster_id=cluster_id, stixel_size=stixel_size)

def compute_interspike_interval_distribution(spike_times: np.ndarray, bin_edges: np.ndarray=np.linspace(0,300,601)) -> np.ndarray:
    """ Computes the interspike intervals from a spike train. 
    Parameters:
        spike_times: object array of spike times.
    Returns:
        isi: vector of interspike intervals.
    """
    isi = np.zeros((spike_times.shape[0],len(bin_edges)-1))
    for ii in range(spike_times.shape[0]):
        isi_tmp = list()
        for jj in range(spike_times.shape[1]):
            sp = spike_times[ii,jj]
            # Check if it's  numpy array
            if isinstance(sp, np.ndarray):
                if sp.shape[0] > 1:
                    isi_tmp.extend(np.diff(sp))
        # Compute the interspike interval
        isi_tmp = np.array(isi_tmp).astype(float)
        if len(isi_tmp) > 0:
            isi_tmp = np.histogram(isi_tmp,bins=bin_edges)[0]
            if np.sum(isi_tmp) > 0:
                isi[ii,:] = isi_tmp / np.sum(isi_tmp)
    return isi

def compute_sta_parameters(sta: np.ndarray) -> Tuple[np.ndarray,np.ndarray,np.ndarray,np.ndarray,np.ndarray,np.ndarray,np.ndarray]:
    # Compute the spatiotemporal RF parameters.
    timecourse_matrix, spatial_maps, significance_maps, hull_parameters, hull_area, hull_vertices, rank1_r2 = compute_spatiotemporal_maps(sta)

    # Normalize the spatial maps.
    spatial_maps = spatial_maps / np.max(np.abs(spatial_maps), axis=(1,2,3))[:,None,None,None]
    spatial_maps = np.nan_to_num(spatial_maps, nan=0.0, posinf=0.0, neginf=0.0) # Replace NaNs with zeros.
    return timecourse_matrix, spatial_maps, significance_maps, hull_parameters, hull_area, hull_vertices, rank1_r2

def upsample_frames(
        frames: np.ndarray, 
        upsample_factor: int) -> np.ndarray:
    """ Upsamples the stimulus frames by an integer factor along the time dimension. 
    Parameters:
        frames: stimulus tensor of size t, y, x, color (float32 in contrast [-1,1]).
        upsample_factor: integer factor by which to upsample the stimulus.
    Returns:
        stim: stimulus tensor of size t*upsample_factor, y, x, color (float32 in contrast [-1,1]).
    """
    nt = frames.shape[0]*upsample_factor
    stimulus_dimensions = np.ndim(frames)
    if stimulus_dimensions == 1:
        frames = frames[:,np.newaxis]
        n_frames,n_phosphors = frames.shape
        stim = np.zeros((nt,n_phosphors), dtype=np.float32)
    elif stimulus_dimensions == 2:
        n_frames,n_phosphors = frames.shape
        stim = np.zeros((nt,n_phosphors), dtype=np.float32)
    else:
        n_frames,height,width,n_phosphors = frames.shape
        stim = np.zeros((nt,height,width,n_phosphors), dtype=np.float32)
    for k in range(nt):
        idx = np.floor(k / upsample_factor).astype(int)
        stim[k,...] = frames[idx,...]
    if stimulus_dimensions == 1:
        stim = stim[:,0]
    return stim

def get_stimulus_frames(protocol_id, parameters, frame_times):
    S = Stimulus()
    if (protocol_id == 'FastNoise') or (protocol_id == 'AdaptNoiseColorSteps'):
        if 'gaussianFilter' in parameters:
            gaussianFilter = bool(parameters['gaussianFilter'])
            filterSdStixels = parameters['filterSdStixels']
        else:
            gaussianFilter = False
            filterSdStixels = 1.0
        frames = S.getFastNoiseFrames(
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
    elif (protocol_id == 'SpatialNoise') or (protocol_id == 'AdaptNoiseSpatial'):
        if 'gaussianFilter' in parameters:
            gaussianFilter = bool(parameters['gaussianFilter'])
            filterSdStixels = parameters['filterSdStixels']
        else:
            gaussianFilter = False
            filterSdStixels = 1.0
        f_times = frame_times - frame_times[0]
        total_frames = len(frame_times)-1 #len(np.where(np.logical_and((f_times > parameters['preTime']),(f_times <= parameters['preTime']+parameters['stimTime'])))[0])
        unique_frames = len(np.where(np.logical_and((f_times > parameters['preTime']),(f_times <= parameters['preTime']+parameters['uniqueTime'])))[0])
        frames = S.get_spatial_noise_frames(
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
    elif (protocol_id == 'PinkNoise'):
        frames = S.getPinkNoiseFrames( 
            parameters['numXChecks'], 
            parameters['numYChecks'],
            len(frame_times)-1,
            parameters['rmsContrast'], 
            parameters['spatialAmplitude'],
            parameters['temporalAmplitude'],
            parameters['chromaticClass'],
            int(parameters['seed']))
    elif (protocol_id == 'SparseNoise'):
        frames = S.get_sparse_noise_frames(
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

def get_spatial_noise_frames(
        protocol_id,
        numXStixels, 
        numYStixels, 
        numXChecks, 
        numYChecks, 
        chromaticClass, 
        numFrames, 
        stepsPerStixel, 
        seed, 
        frameDwell,
        pixelDensity=0.01) -> np.ndarray:
    S = Stimulus()
    if protocol_id == 'SparseNoise':
        frames = S.get_sparse_noise_frames(int(numXStixels), int(numYStixels), int(numXChecks), int(numYChecks), chromaticClass,
                numFrames, int(stepsPerStixel), int(seed), int(frameDwell), pixelDensity)
    else:
        frames = S.getFastNoiseFrames(int(numXStixels), int(numYStixels), int(numXChecks), int(numYChecks), chromaticClass,
            numFrames, int(stepsPerStixel), int(seed), int(frameDwell), gaussianFilter=False, filterSdStixels=1.0)
    # frames = S.get_spatial_noise_frames(int(numXStixels), int(numYStixels), int(numXChecks), int(numYChecks), chromaticClass, 
    #     9590, 1200, int(stepsPerStixel), int(seed), int(frameDwell), gaussianFilter=False, filterSdStixels=1.0)
    return frames

def downsample_data(frames: np.ndarray, binned_spikes: np.ndarray, down_factor: int) -> Tuple[np.ndarray, np.ndarray]:
    """ Downsamples the stimulus frames and binned spikes by an integer factor along the time dimension.
    Parameters:
        frames: stimulus tensor of size t, y, x, color (float32 in contrast [-1,1]).
        binned_spikes: spiketrain matrix of size t,cells
        down_factor: integer factor by which to downsample the stimulus.
    Returns:
        frames: stimulus tensor of size t/down_factor, y, x, color (float32 in contrast [-1,1]).
        binned_spikes: spiketrain matrix of size t/down_factor,cells
    """
    n_frames = np.floor(binned_spikes.shape[0] / down_factor).astype(int)
    tmp_frames = np.zeros((n_frames,frames.shape[1],frames.shape[2],frames.shape[3]))
    tmp_spikes = np.zeros((n_frames, binned_spikes.shape[1]))
    for frame in range(n_frames):
        fr_idx = frame*down_factor + np.arange(down_factor)
        tmp_frames[frame, ...] = np.mean(frames[fr_idx, ...],axis=0)
        tmp_spikes[frame, :] = np.sum(binned_spikes[fr_idx, :],axis=0)
    binned_spikes = tmp_spikes
    frames = tmp_frames
    return frames, binned_spikes

def compute_sta(
        stimulus: np.ndarray, 
        binned_spikes: np.ndarray, 
        depth: int=61) -> np.ndarray:
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
    stimulus_dimensions = np.ndim(stimulus) 
    # Intialize STA tensor.
    num_cells = binned_spikes.shape[1]
    if stimulus_dimensions == 1:
        stimulus = stimulus[:,np.newaxis]
        n_frames,n_phosphors = stimulus.shape
        sta_tensor = np.zeros((num_cells,depth,n_phosphors))
        stimulus = stimulus.astype(np.float32).T
    elif stimulus_dimensions == 2:
        n_frames,n_phosphors = stimulus.shape
        sta_tensor = np.zeros((num_cells,depth,n_phosphors))
        stimulus = stimulus.astype(np.float32).T
    else:
        n_frames,height,width,n_phosphors = stimulus.shape
        sta_tensor = np.zeros((num_cells,depth,height,width,n_phosphors))
        # Vectorize stimulus.
        stimulus = np.reshape(stimulus,(n_frames,
            height * width * n_phosphors)).astype(np.float32).T
    # Set up torch stuff.
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu') 
    # Initialize STA, the lag range, and compute.
    lags = np.arange(0,depth,1)
    binned_spikes = torch.tensor(binned_spikes).to(device).to(torch.float32)
    stimulus = torch.tensor(stimulus).to(device).to(torch.float32)
    for lag_count, lag in enumerate(lags):
        with torch.no_grad():
            sta = torch.matmul(stimulus[:,:n_frames-lag],
                            binned_spikes[lag:,:])
            # Index appropriately to avoid nans and normalize by number of spikes.
            nonzero_inds = torch.nonzero(torch.sum(binned_spikes[lag:,:],
                                                dim=0),as_tuple=False)
            sta[:,nonzero_inds] /= torch.sum(binned_spikes[lag:,:],
                                            dim=0)[nonzero_inds]
        sta = sta.detach().cpu().numpy()
        if stimulus_dimensions > 2:
            sta_tensor[:,lag_count,...] = np.reshape(sta.T,
                                                (num_cells,height,
                                                width,n_phosphors))
        else:
            sta_tensor[:,lag_count] = sta.T
    binned_spikes = binned_spikes.detach().cpu().numpy()
    # Multiply by the spike count for each cell; average later.
    if stimulus_dimensions <= 2:
        sta_tensor *= np.sum(binned_spikes,axis=0)[:,None,None]
    else:
        sta_tensor *= np.sum(binned_spikes,axis=0)[:,None,None,None,None]
    if device.type == 'cuda':
        torch.cuda.empty_cache()
    # Clear some memory
    del sta
    gc.collect()
    return sta_tensor

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

def compute_nonlinearity(
        protocol_id: str,
        params: dict, 
        sta: np.ndarray,
        binned_spikes: np.ndarray, 
        unique_frames: int, 
        stride: int=2, 
        num_bins: int=25,
        crop_fraction: float=1.0) -> Tuple[np.ndarray, np.ndarray]:
    # Allocate memory to the generator and response.
    generator = np.zeros((len(cluster_id),1))
    response = np.zeros((len(cluster_id),1))
    device_type='gpu'
    
    # Determine the number of processing chunks.
    if params['numXChecks'][0] > (200*crop_fraction):
        num_chunks = np.ceil((float(params['numXChecks'][0])/100.0)**2).astype(int)
        cells_per_chunk = np.ceil(float(binned_spikes.shape[0]) / float(num_chunks)).astype(int)
    else:
        num_chunks = 1
        cells_per_chunk = binned_spikes.shape[0]
    for idx in tqdm(range(len(params['numXStixels'])), desc='Computing Nonlinearity'):
        if protocol_id == 'SparseNoise':
            pixelDensity = float(params['pixelDensity'][idx])
        else:
            pixelDensity = 0.0
        frames = get_spatial_noise_frames(
            protocol_id,
            int(params['numXStixels'][idx]),
            int(params['numYStixels'][idx]), 
            int(params['numXChecks'][idx]),
            int(params['numYChecks'][idx]),
            params['chromaticClass'][idx],
            unique_frames,
            int(params['stepsPerStixel'][idx]), 
            int(params['seed'][idx]),
            int(params['frameDwell'][idx]),
            pixelDensity=pixelDensity)
        if crop_fraction < 1.0:
            frames = crop_frames(frames=frames, crop_fraction=crop_fraction)
        
        spikes_tmp = binned_spikes[:,idx,:].T
        if (protocol_id == 'SparseNoise'):
            frame_dwell = int(params['frameDwell'][idx])
            frames, spikes_tmp = downsample_data(frames=frames, binned_spikes=spikes_tmp, down_factor=frame_dwell)
        # Append the spike responses.
        response = np.append(response, spikes_tmp.T, axis = 1)
        if stride > 1:
            frames = upsample_frames(frames, stride)
        if num_chunks == 1:
            gen_tmp = compute_generator_signal_torch(frames, sta, stride=1, alpha=1/10000, device_type=device_type)
        else:
            gen_tmp = np.zeros((sta.shape[0],frames.shape[0]))
            for ii in range(num_chunks):
                cell_idx = np.arange(ii*cells_per_chunk, np.min(((ii+1)*cells_per_chunk, binned_spikes.shape[0])))
                gen_tmp[cell_idx,:] = compute_generator_signal_torch(frames, sta[cell_idx,...], stride=1, alpha=1/10000, device_type=device_type)
        generator = np.append(generator, gen_tmp, axis = 1)
    x_bin, y_bin = bin_nonlinearities(generator=generator[:,1:], response=response[:,1:], num_bins=num_bins)
    # hdf5storage.savemat('/data/data/sorted/20240717C/sparse/kilosort2.5/generator.mat', {'generator':generator,'response':response,'sta':sta}, format=7.3, matlab_compatible=True, compress=False )
    # Get the nonzero indices.
    del gen_tmp, spikes_tmp, generator, response
    gc.collect()
    return x_bin, y_bin

def compute_spatial_sta(
        protocol_id: str,
        params: dict, 
        binned_spikes: np.ndarray, 
        unique_frames: int, 
        stride: int=2, 
        depth: int=61,
        crop_fraction: float=1.0) -> Tuple[np.ndarray, np.ndarray]:
    # Determine the number of processing chunks.
    if params['numXChecks'][0] > (200*crop_fraction):
        num_chunks = np.ceil((float(params['numXChecks'][0])/100.0)**2).astype(int)
        cells_per_chunk = np.ceil(float(binned_spikes.shape[0]) / float(num_chunks)).astype(int)
    else:
        num_chunks = 1
        cells_per_chunk = binned_spikes.shape[0]
    for idx in tqdm(range(len(params['numXStixels'])), desc='Computing STA'):
        if protocol_id == 'SparseNoise':
            pixelDensity = float(params['pixelDensity'][idx])
        else:
            pixelDensity = 0.0
        frames = get_spatial_noise_frames(
            protocol_id,
            int(params['numXStixels'][idx]),
            int(params['numYStixels'][idx]), 
            int(params['numXChecks'][idx]),
            int(params['numYChecks'][idx]),
            params['chromaticClass'][idx],
            unique_frames,
            int(params['stepsPerStixel'][idx]), 
            int(params['seed'][idx]),
            int(params['frameDwell'][idx]),
            pixelDensity=pixelDensity)
        if crop_fraction < 1.0:
            frames = crop_frames(frames=frames, crop_fraction=crop_fraction)
        if idx == 0:
            sta = np.zeros((binned_spikes.shape[0],depth,frames.shape[1],frames.shape[2],3))
        
        spikes_tmp = binned_spikes[:,idx,:].T
        if (protocol_id == 'SparseNoise'):
            frame_dwell = int(params['frameDwell'][idx])
            frames, spikes_tmp = downsample_data(frames=frames, binned_spikes=spikes_tmp, down_factor=frame_dwell)
        if stride > 1:
            frames = upsample_frames(frames, stride)
        if num_chunks == 1:
            sta_tmp = compute_sta(stimulus=frames.astype(np.float32), binned_spikes=spikes_tmp, depth=depth)
            sta += sta_tmp
        else:
            for ii in range(num_chunks):
                cell_idx = np.arange(ii*cells_per_chunk, np.min(((ii+1)*cells_per_chunk, binned_spikes.shape[0])))
                sta_tmp = compute_sta(stimulus=frames.astype(np.float32), binned_spikes=spikes_tmp[:,cell_idx], depth=depth)
                sta[cell_idx,...] += sta_tmp
    spike_count = np.sum(binned_spikes,axis=(1,2))
    # Get the nonzero indices.
    nonzero_inds = np.nonzero(spike_count)[0]
    # Divide by spike_count-1 to get the variance.
    sta[nonzero_inds,...] /= spike_count[nonzero_inds][:,None,None,None,None]
    del sta_tmp, spikes_tmp
    gc.collect()
    sta[np.isnan(sta)]=0.0 # Set NaNs to zero
    return sta, spike_count

def compute_stv_old(sta: np.ndarray, stimulus: np.ndarray, binned_spikes: np.ndarray, depth: int=61):
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
        STV tensor, indexed by cell,t,y,x,color channel
    """
    stimulus_dimensions = np.ndim(stimulus) 
    # Intialize STA tensor.
    num_cells = binned_spikes.shape[1]
    if stimulus_dimensions == 1:
        stimulus = stimulus[:,np.newaxis]
        n_frames,n_phosphors = stimulus.shape
        stv_tensor = np.zeros((num_cells,depth,n_phosphors))
        stimulus = stimulus.astype(np.float32).T
    elif stimulus_dimensions == 2:
        n_frames,n_phosphors = stimulus.shape
        stv_tensor = np.zeros((num_cells,depth,n_phosphors))
        stimulus = stimulus.astype(np.float32).T
    else:
        n_frames,height,width,n_phosphors = stimulus.shape
        stv_tensor = np.zeros((num_cells,depth,height,width,n_phosphors))
        sta_tensor = np.reshape(sta,(num_cells,depth,height*width*n_phosphors))
        # Vectorize stimulus.
        stimulus = np.reshape(stimulus,(n_frames,
            height * width * n_phosphors)).astype(np.float32).T
    # Set up torch stuff.
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu') 
    # Initialize STA, the lag range, and compute.
    lag_count = 0
    lags = np.arange(0,depth,1)
    binned_spikes = torch.tensor(binned_spikes).to(device).to(torch.float32)
    stimulus = torch.tensor(stimulus).to(device).to(torch.float32)
    sta_tensor = torch.tensor(sta_tensor).to(device).to(torch.float32)
    for lag in lags:
        stv = torch.zeros((height*width*n_phosphors, num_cells)).to(device).to(torch.float32)
        with torch.no_grad():
            for cell_count in range(num_cells):
                stv[:,cell_count] = torch.matmul(
                    torch.pow( stimulus[:,:n_frames-lag] - sta_tensor[cell_count,lag,:][:,None], 2),
                                binned_spikes[lag:,cell_count])
            # Index appropriately to avoid nans and normalize by number of spikes.
            nonzero_inds = torch.nonzero(torch.sum(binned_spikes[lag:,:],
                                                   dim=0),as_tuple=False)
            stv[:,nonzero_inds] /= torch.sum(binned_spikes[lag:,:],
                                            dim=0)[nonzero_inds]
        stv = stv.detach().cpu().numpy()
        if stimulus_dimensions > 2:
            stv_tensor[:,lag_count,...] = np.reshape(stv.T,
                                                (num_cells,height,
                                                width,n_phosphors))
        else:
            stv_tensor[:,lag_count] = stv.T
        lag_count +=1
    if device.type == 'cuda':
        torch.cuda.empty_cache()
    # Clear some memory
    del stv
    gc.collect()
    return stv_tensor

def compute_stv_test(sta: np.ndarray, stimulus: np.ndarray, binned_spikes: np.ndarray, depth: int=61):
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
        STV tensor, indexed by cell,t,y,x,color channel
    """
    stimulus_dimensions = np.ndim(stimulus) 
    # Intialize STA tensor.
    num_cells = binned_spikes.shape[1]
    if stimulus_dimensions == 1:
        stimulus = stimulus[:,np.newaxis]
        n_frames,n_phosphors = stimulus.shape
        stv_tensor = np.zeros((num_cells,depth,n_phosphors))
        stimulus = stimulus.astype(np.float32).T
    elif stimulus_dimensions == 2:
        n_frames,n_phosphors = stimulus.shape
        stv_tensor = np.zeros((num_cells,depth,n_phosphors))
        stimulus = stimulus.astype(np.float32).T
    else:
        n_frames,height,width,n_phosphors = stimulus.shape
        stv_tensor = np.zeros((num_cells,depth,height,width,n_phosphors))
        # Vectorize stimulus.
        stimulus = np.reshape(stimulus,(n_frames,
            height * width * n_phosphors)).astype(np.float32).T
    # Set up torch stuff.
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu') 
    # Initialize STA, the lag range, and compute.
    lags = np.arange(0,depth)
    binned_spikes = torch.tensor(binned_spikes).to(device).to(torch.float32)
    stimulus = torch.tensor(stimulus).to(device).to(torch.float32)
    sta_tensor = np.reshape(sta,(num_cells,depth,height*width*n_phosphors))
    sta_tensor = torch.tensor(sta_tensor).to(device).to(torch.float32)
    for ii in range(sta_tensor.shape[0]):
        if torch.sum(binned_spikes[lags[-1]:,ii]) > 0:
            for lag_count, lag in enumerate(lags):
                with torch.no_grad():
                    sp_count = torch.sum(binned_spikes[lag:,ii])
                    if sp_count > 0:
                        s_tmp = stimulus[:,:n_frames-lag] - sta_tensor[ii,lag,:][:,None]
                        stv = torch.matmul(torch.mul(s_tmp, s_tmp),
                                binned_spikes[lag:,ii])
                        stv /= sp_count
                # stv = stv.detach().cpu().numpy()
                if stimulus_dimensions > 2:
                    stv_tensor[ii,lag_count,...] = np.reshape(stv.detach().cpu().numpy().T, (height,width,n_phosphors))
                else:
                    stv_tensor[:,lag_count] = stv.detach().cpu().numpy().T
    binned_spikes = binned_spikes.detach().cpu().numpy()
    # Multiply by the spike count for each cell; average later.
    if stimulus_dimensions <= 2:
        stv_tensor *= np.sum(binned_spikes,axis=0)[:,None,None]
    else:
        stv_tensor *= np.sum(binned_spikes,axis=0)[:,None,None,None,None]
    if device.type == 'cuda':
        torch.cuda.empty_cache()
    # Clear some memory
    del stv
    gc.collect()
    return stv_tensor

def compute_spatial_stv_test(params: dict, sta: np.ndarray, binned_spikes: np.ndarray, stride: int=2, depth: int=61) -> np.ndarray:
    stv = np.zeros((binned_spikes.shape[0],depth,int(params['numYChecks'][0]),int(params['numXChecks'][0]),3))
    for idx in tqdm(range(len(params['numXStixels'])), desc='Computing STV'):
        frames = get_spatial_noise_frames(
            int(params['numXStixels'][idx]),
            int(params['numYStixels'][idx]), 
            int(params['numXChecks'][idx]),
            int(params['numYChecks'][idx]),
            params['chromaticClass'][idx],
            int(params['numFrames'][idx]),
            int(params['stepsPerStixel'][idx]), 
            int(params['seed'][idx]),
            int(params['frameDwell'][idx]))
        if stride > 1:
            frames = upsample_frames(frames, stride)
        stv_tmp = compute_stv_test(sta=sta, stimulus=frames.astype(np.float32), binned_spikes=binned_spikes[:,idx,:].T, depth=depth)
        stv += stv_tmp
    spike_count = np.sum(binned_spikes,axis=(1,2))
    # Get the nonzero indices.
    nonzero_inds = np.nonzero(spike_count-1)[0]
    # Divide by spike_count-1 to get the variance.
    stv[nonzero_inds,...] /= (spike_count[nonzero_inds]-1)[:,None,None,None,None]
    # Subtract the STA from the STV (after multiplying by the spike count).
    # stv[nonzero_inds,...] -= (sta[nonzero_inds,...]*sta[nonzero_inds,...]) * (spike_count[nonzero_inds]/(spike_count[nonzero_inds]-1))[:,None,None,None,None]
    del stv_tmp
    gc.collect()
    stv[np.isnan(stv)]=0.0 # Set NaNs to zero
    return stv

def compute_stv(stimulus: np.ndarray, binned_spikes: np.ndarray, depth: int=61):
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
        STV tensor, indexed by cell,t,y,x,color channel
    """
    stimulus_dimensions = np.ndim(stimulus) 
    # Intialize STA tensor.
    num_cells = binned_spikes.shape[1]
    if stimulus_dimensions == 1:
        stimulus = stimulus[:,np.newaxis]
        n_frames,n_phosphors = stimulus.shape
        stv_tensor = np.zeros((num_cells,depth,n_phosphors))
        stimulus = stimulus.astype(np.float32).T
    elif stimulus_dimensions == 2:
        n_frames,n_phosphors = stimulus.shape
        stv_tensor = np.zeros((num_cells,depth,n_phosphors))
        stimulus = stimulus.astype(np.float32).T
    else:
        n_frames,height,width,n_phosphors = stimulus.shape
        stv_tensor = np.zeros((num_cells,depth,height,width,n_phosphors))
        # Vectorize stimulus.
        stimulus = np.reshape(stimulus,(n_frames,
            height * width * n_phosphors)).astype(np.float32).T
    # Set up torch stuff.
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu') 
    # Initialize STA, the lag range, and compute.
    lags = np.arange(0,depth)
    binned_spikes = torch.tensor(binned_spikes).to(device).to(torch.float32)
    stimulus = torch.tensor(stimulus*stimulus).to(device).to(torch.float32)
    # Square the stimulus.
    # stimulus = torch.pow(stimulus,2)
    # stimulus = torch.mul(stimulus, stimulus)
    for lag_count, lag in enumerate(lags):
        with torch.no_grad():
            stv = torch.matmul(stimulus[:,:n_frames-lag],
                            binned_spikes[lag:,:])
            # Index appropriately to avoid nans and normalize by number of spikes.
            nonzero_inds = torch.nonzero(torch.sum(binned_spikes[lag:,:],
                                                dim=0),as_tuple=False)
            stv[:,nonzero_inds] /= torch.sum(binned_spikes[lag:,:],
                                            dim=0)[nonzero_inds]
        stv = stv.detach().cpu().numpy()
        if stimulus_dimensions > 2:
            stv_tensor[:,lag_count,...] = np.reshape(stv.T,
                                                (num_cells,height,
                                                width,n_phosphors))
        else:
            stv_tensor[:,lag_count] = stv.T
    binned_spikes = binned_spikes.detach().cpu().numpy()
    # Multiply by the spike count for each cell; average later.
    if stimulus_dimensions <= 2:
        stv_tensor *= np.sum(binned_spikes,axis=0)[:,None,None]
    else:
        stv_tensor *= np.sum(binned_spikes,axis=0)[:,None,None,None,None]
    if device.type == 'cuda':
        torch.cuda.empty_cache()
    # Clear some memory
    del stv
    gc.collect()
    return stv_tensor

def compute_spatial_stv(params: dict, sta: np.ndarray, binned_spikes: np.ndarray, stride: int=2, depth: int=61) -> np.ndarray:
    stv = np.zeros((binned_spikes.shape[0],depth,int(params['numYChecks'][0]),int(params['numXChecks'][0]),3))
    for idx in tqdm(range(len(params['numXStixels'])), desc='Computing STV'):
        frames = get_spatial_noise_frames(
            int(params['numXStixels'][idx]),
            int(params['numYStixels'][idx]), 
            int(params['numXChecks'][idx]),
            int(params['numYChecks'][idx]),
            params['chromaticClass'][idx],
            int(params['numFrames'][idx]),
            int(params['stepsPerStixel'][idx]), 
            int(params['seed'][idx]),
            int(params['frameDwell'][idx]))
        if stride > 1:
            frames = upsample_frames(frames, stride)
        stv_tmp = compute_stv(stimulus=frames.astype(np.float32), binned_spikes=binned_spikes[:,idx,:].T, depth=depth)
        stv += stv_tmp
    spike_count = np.sum(binned_spikes,axis=(1,2))
    # Get the nonzero indices.
    nonzero_inds = np.nonzero(spike_count-1)[0]
    # Divide by spike_count-1 to get the variance.
    stv[nonzero_inds,...] /= (spike_count[nonzero_inds]-1)[:,None,None,None,None] #(spike_count[nonzero_inds]-1)[:,None,None,None,None]
    # Subtract the STA from the STV (after multiplying by the spike count).
    # stv[nonzero_inds,...] -= (sta[nonzero_inds,...]*sta[nonzero_inds,...]) * (spike_count[nonzero_inds]/(spike_count[nonzero_inds]-1))[:,None,None,None,None]
    del stv_tmp
    gc.collect()
    stv[np.isnan(stv)]=0.0 # Set NaNs to zero
    return stv

def compute_spatial_stv_old(params: dict, sta: np.ndarray, binned_spikes: np.ndarray, stride: int=2, depth: int=61) -> np.ndarray:
    stv = np.zeros((binned_spikes.shape[0],depth,int(params['numYChecks'][0]),int(params['numXChecks'][0]),3))
    for idx in tqdm(range(len(params['numXStixels'])), desc='Computing STV'):
        frames = get_spatial_noise_frames(
            int(params['numXStixels'][idx]),
            int(params['numYStixels'][idx]), 
            int(params['numXChecks'][idx]),
            int(params['numYChecks'][idx]),
            params['chromaticClass'][idx],
            int(params['numFrames'][idx]),
            int(params['stepsPerStixel'][idx]), 
            int(params['seed'][idx]),
            int(params['frameDwell'][idx]))
        if stride > 1:
            frames = upsample_frames(frames, stride)
        stv_tmp = compute_stv(sta=sta, stimulus=frames.astype(np.float32), binned_spikes=binned_spikes[:,idx,:].T, depth=depth)
        stv += stv_tmp
    spike_count = np.sum(binned_spikes,axis=(1,2))
    for ii in range(stv.shape[0]):
        if spike_count[ii] > 0:
            stv[ii,...] /= spike_count[ii]
    del stv_tmp
    gc.collect()
    stv[np.isnan(stv)]=0.0 # Set NaNs to zero
    return stv

def normalize_sta(sta: np.ndarray, norm_type: str='max') -> np.ndarray:
    for i in range(sta.shape[0]):
        if norm_type == 'max':
            v = np.max(np.abs(sta[i,...]))
        if v > 0:
            sta[i] /= v
    sta[np.isnan(sta)]=0.0 # Set NaNs to zero
    return sta

def compute_stv_parameters(stv: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    n_cells,stv_depth,x,y,n_phosphors = stv.shape
    timecourse_matrix = np.zeros((n_cells,stv_depth,n_phosphors))
    spatial_maps = np.zeros((n_cells,x,y,n_phosphors))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        for i in tqdm(range(n_cells), desc ="Computing STV parameters"):
            cell_stv = stv[i,...].copy()
            cell_stv -= np.nanmean(cell_stv)
            cell_stv /= np.max(np.abs(cell_stv))
            # Compute the time course.
            time_course = calculate_time_course_svd(cell_stv, threshold=2.0)
            space_map = compute_spatial_map(cell_stv, time_course)
            timecourse_matrix[i,...] = time_course
            spatial_maps[i,...] = space_map
    return timecourse_matrix, spatial_maps



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Spatiotemporal noise analysis')
    parser.add_argument('experimentName', type=str, help='Name of experiment')
    parser.add_argument('-a','--algorithm', default='kilosort2', type=str, help='Sorting algorithm used (yass or kilosort2)')
    parser.add_argument('-p','--protocol', default='SpatialNoise', type=str, help='Name of stimulus protocol to analyze')
    parser.add_argument('-g','--group', default=None, type=str, help='Search string for EpochGroup (optional)')
    parser.add_argument('-f','--file', nargs='+', default=None, type=str, help='Data file identifier')
    parser.add_argument('-c','--chunk_name', default=None, type=str, help='Name of chunk containing the file(s)')
    parser.add_argument('-x','--crop_fraction', default=1.0, type=float, help='Fraction of the stixels to use for the calculations')
    parser.add_argument('-s','--stride', default=2, type=int, help='Bins per frame')
    parser.add_argument('-d','--depth', default=61, type=int, help='Number of frames to include in the STA')
    parser.add_argument('-n','--nonlinearity', default=False, action="store_true", help='compute the nonlinearity?')
    parser.add_argument('-v','--compute_stv', default=False, action="store_true", help='compute the spike-triggered variance?')
    parser.add_argument('-r','--sample_rate', default=20000.0, type=float, help='Data sample rate')

    args = parser.parse_args()

    if type(args.file) is str:
        args.file = [args.file]

    stride = args.stride
    depth = args.depth
    experiment_name = args.experimentName
    protocol_id = args.protocol
    group_id = args.group
    sort_algorithm=args.algorithm
    file_name=args.file
    sample_rate = args.sample_rate
    chunk_name = args.chunk_name
    crop_fraction = args.crop_fraction
    nonlinearity = args.nonlinearity

    if protocol_id == 'SparseNoise':
        stride = 1
        depth = 3
        nonlinearity = True

    if int(experiment_name[:8]) < 20230926:
        marginal_frame_rate = 60.31807657 # Upper bound on the frame rate to make sure that we don't miss any frames.
        frame_offset = stride
    else:
        marginal_frame_rate = 59.941548817817917 # Upper bound on the frame rate to make sure that we don't miss any frames.
        frame_offset = 0

    if file_name != None:
        if type(file_name) is str:
            file_name = [file_name]
        if type(file_name) is list:
            label = '_'.join(file_name)
        else:
            label = file_name
    else:
        label = ''

    # Get the output file path.
    SORT_PATH, _ = cfg.get_data_paths()
    SAVE_PATH = cfg.get_save_path()
    file_dir = SORT_PATH + experiment_name + '/'
    file_str = experiment_name + '_' + sort_algorithm + '_' + label

    if chunk_name is None:
        out_file_path = file_dir + file_str
        globals_path = file_dir
        globals_name = file_str
    else:
        out_file_path = file_dir + chunk_name + '/' + sort_algorithm + '/' + sort_algorithm
        globals_path = file_dir + chunk_name + '/' + sort_algorithm + '/'
        globals_name = sort_algorithm
    
    # Load the dataset.
    d = Dataset(experiment_name)

    # Pull the frame times.
    frame_times = d.get_frame_times(protocolStr=protocol_id, file_name=file_name)

    if protocol_id == 'FastNoise':
        param_names = ['stimulusClass', 'chromaticClass', 'noiseClass', 'frameDwell', 'stepDuration', 'gridSize',
            'numXChecks', 'numYChecks', 'numFrames', 'numXStixels', 'numYStixels', 'stixelSize', 'stepsPerStixel', 'seed', 'stimTime']
    else:
        param_names = ['stimulusClass', 'chromaticClass', 'noiseClass', 'frameDwell', 'stepDuration', 'gridSize', 'uniqueTime', 'repeatTime',
            'numXChecks', 'numYChecks', 'numFrames', 'numXStixels', 'numYStixels', 'stixelSize', 'stepsPerStixel', 'seed', 'stimTime','pixelDensity']

    print('Computing spike counts...')
    spike_dict, cluster_id, params, unique_params, pre_pts, stim_pts, tail_pts, mean_frame_rate = d.get_count_at_frame_multiple(
        protocolStr=protocol_id, groupStr=group_id, param_names=param_names, sort_algorithm=sort_algorithm, file_name=file_name, frame_rate=marginal_frame_rate, stride=stride)
    
    # Get the array properties.
    array_id, num_samples = get_array_properties(SORT_PATH, experiment_name, file_name)
    # array_id = 3501
    if num_samples == 0:
        num_samples = d.vcd.n_samples
    
    pre_pts = np.array(pre_pts).astype(int)
    stim_pts = np.array(stim_pts).astype(int)
    tail_pts = np.array(tail_pts).astype(int)

    if len(params['gridSize']) > 0:
        grid_size = float(params['gridSize'][0]) # Get the size of the underlying grid in microns.
    else:
        grid_size = 30.0

    # Determine the number of bins to analyze.
    if protocol_id == 'FastNoise':
        unique_frames = int(params['numFrames'][0])
    else:
        if len(params['uniqueTime']) > 0:
            unique_frames = int(np.ceil(float(params['uniqueTime'][0]) / 1000.0 * marginal_frame_rate))
        else:
            unique_frames = int(params['numFrames'][0])
    num_bins = int(unique_frames*stride)

    # Get the spike counts for each cluster.
    psth = np.zeros((len(cluster_id), spike_dict[cluster_id[0]].shape[0], spike_dict[cluster_id[0]].shape[1]))
    for ii in range(len(cluster_id)):
        psth[ii,:,:] = spike_dict[cluster_id[ii]]

    ## STA FILE
    # Compute the spike-triggered average (STA).
    sta, spike_count = compute_spatial_sta(
        protocol_id=protocol_id,
        params=params, 
        binned_spikes=psth[:,:,pre_pts[0]+frame_offset:pre_pts[0]+frame_offset+num_bins], 
        unique_frames=unique_frames, 
        stride=stride, 
        depth=depth,
        crop_fraction=crop_fraction)
    sta = normalize_sta(sta, norm_type='max') # Normalize the STA.
    sta_size=sta.shape[1:] # Get the size of the STA (t,y,x,c).
    # Write the STA file.
    if chunk_name is not None:
        write_sta_file(out_path=out_file_path, sta=sta, ste=None, cluster_id=cluster_id, stixel_size=grid_size)

    ## PARAMETERS FILE
    try:
        # Compute the ISI distribution.
        spike_times, cluster_id_times, _, _, _, _, _ = d.get_spike_times_and_parameters(
            protocolStr=protocol_id, groupStr=group_id, param_names=param_names, sort_algorithm=sort_algorithm, file_name=file_name, bin_rate=1000.0, sample_rate=sample_rate)
        isi = compute_interspike_interval_distribution(spike_times=spike_times, bin_edges=np.linspace(0,300,601))
    
        # Compute the STA parameters.
        timecourse_matrix, spatial_maps, significance_maps, hull_parameters, hull_area, hull_vertices, rank1_r2 = compute_sta_parameters(sta)
        
        # Compute the Gaussian fit to the STA.
        from fit_sta import compute_rf_parameters
        gauss_params = np.zeros((sta.shape[0], 6))
        for ii in range(sta.shape[0]):
            try:
                gauss_params[ii,:] = compute_rf_parameters(sta[ii,...])
            except Exception as error:
                gauss_params[ii,1:] = hull_parameters[ii,:] # Use the hull parameters if the Gaussian fit fails.
        

        # Compute the spike-triggered variance (STV).
        if args.compute_stv:
            # stv = compute_spatial_stv(params=params, sta=sta, binned_spikes=psth[:,:,pre_pts[0]:pre_pts[0]+num_bins], unique_frames=unique_frames, stride=stride, depth=depth)
            # timecourse_variance, spatial_variance = compute_stv_parameters(stv.copy())
            # stv = normalize_sta(stv, norm_type='max') # Normalize the STV.
            # For binary stimuli.
            timecourse_variance = -np.sqrt(timecourse_matrix*timecourse_matrix)
            stv = 1 - (sta*sta)
        
        # Compute the RF contour.
        max_contour_values=30
        simple_contour_xy = np.zeros((spatial_maps.shape[0], max_contour_values, 2))
        simple_contour_area = np.zeros((spatial_maps.shape[0]))
        for c_idx in range(spatial_maps.shape[0]):
            try:
                simple_contour_xy[c_idx,...], simple_contour_area[c_idx] = fit_rf_hull_and_contour(space_map=spatial_maps[c_idx,...], max_contour_values=max_contour_values)
            except Exception as error:
                print('An error occurred while computing the RF contour:', type(error).__name__, '-', error)
        # Flip the y coordinates.
        simple_contour_xy[:, :, 1] = sta_size[1] - simple_contour_xy[:, :, 1]

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

        # Write the params file.
        if chunk_name is not None:
            write_params_file(out_path=out_file_path, cluster_id=cluster_id, timecourse_matrix=timecourse_matrix, isi=isi, spike_count=spike_count, 
                hull_parameters=hull_parameters, isi_binning=0.5, contour_xy=hull_vertices, contour_area=hull_area, simple_contour_xy=simple_contour_xy, 
                simple_contour_area=simple_contour_area, timecourse_variance=None, ei_parameters=None, sta_size=sta_size)
    except Exception as error:
        print('An error occurred while computing the STA parameters:', type(error).__name__, '-', error)

    ## GLOBALS FILE
    try:
        sta_width = int(sta.shape[3])
        sta_height = int(sta.shape[2])
        if int(experiment_name[:8]) < 20230926:
            pixelsPerStixel = int(800.0 / sta_width)
            mu_per_pixel = 3.8
        else:
            pixelsPerStixel = int(912.0 * 2.0 / sta_width)
            mu_per_pixel = 3.34
        micronsPerStixel = grid_size #np.round(pixelsPerStixel * mu_per_pixel).astype(float)
        
        if chunk_name is not None:
            write_globals_file(globals_path, globals_name, array_id, num_samples, sta_width, sta_height, micronsPerStixel, pixelsPerStixel, mean_frame_rate, stride)
    except Exception as error:
        print('An error occurred while writing the globals file:', type(error).__name__, '-', error)
    
    ## Nonlinearity
    if nonlinearity: # Compute the nonlinearity.
        nonlinearity_bins = 25
        x_bin, y_bin = compute_nonlinearity(
            protocol_id=protocol_id,
            params=params,
            sta=sta,
            binned_spikes=psth[:,:,pre_pts[0]+frame_offset:pre_pts[0]+frame_offset+num_bins],  
            unique_frames=unique_frames, 
            stride=stride, 
            num_bins=nonlinearity_bins,
            crop_fraction=crop_fraction)
    else:
        x_bin = []
        y_bin = []
    
    # Save the data.
    mdic = {'cluster_id': cluster_id, 'spike_count': spike_count, 'isi':isi,
        'timecourse_matrix': timecourse_matrix, 'spatial_maps': spatial_maps, 
        'significance_maps': significance_maps, 'hull_area': hull_area, 'hull_vertices': hull_vertices, 
        'hull_parameters': hull_parameters, 'rank1_r2': rank1_r2, 'x_bin': x_bin, 'y_bin': y_bin,
        'simple_contour_xy': simple_contour_xy, 'simple_contour_area': simple_contour_area } 
    
    if args.compute_stv:
        mdic.update({'stv': stv, 'timecourse_variance': timecourse_variance})
    
    if os.path.exists(out_file_path + '_params.mat'):
        os.remove(out_file_path + '_params.mat')
    hdf5storage.savemat(out_file_path + '_params.mat', mdic, format=7.3, matlab_compatible=True, compress=False )