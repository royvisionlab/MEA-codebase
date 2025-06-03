# Code for computing cross-correlograms and autocorrelograms.
# Written by: Mike Manookin, University of Washington, 2023.

import os
import numpy as np
from typing import Union, Tuple
import visionloader as vl
from tqdm import tqdm
import correlations.correlation_cpp_extensions as corrcpp

def get_spike_times(
        sort_path: str, 
        sorter_name: str, 
        file_names: Union[str,list], 
        sample_rate: float=20000.0) -> Tuple[np.ndarray, np.ndarray]:
    """ Get the spike times for a given sorter and file name. 
    Arguments:
        sort_path: Path to the sorted data.
        sorter_name: Name of the sorter.
        file_names: Name of the file(s) to load. Can be a string or a list of strings.
        sample_rate: Sampling rate of the data in Hz (default: 20000.0).
    Returns:
        spike_times: Spike times for the cells. Object array: [num cells] with each element being an array of spike times for that cell. [num spikes]
        cell_ids: Cell IDs for the cells. [num cells]
    """
    if isinstance(file_names, str):
        file_names = [file_names]
    spike_dict = dict()
    spike_offset = 0.0
    for file_name in file_names:
        file_path = os.path.join(sort_path, file_name, sorter_name)
        v_table = vl.load_vision_data(file_path, file_name, include_neurons=True)
        cell_ids = v_table.get_cell_ids()
        max_spike = 0.0
        for cell_id in cell_ids:
            sp = v_table.get_spike_times_for_cell(cell_id)/sample_rate*1000.0 # Spike times in ms. 
            if len(sp) > 0:
                max_spike = np.max([max_spike, np.max(sp)])
                if cell_id in spike_dict.keys():
                    spike_dict[cell_id] = np.concatenate((spike_dict[cell_id], sp + spike_offset))
                else:
                    spike_dict[cell_id] = sp + spike_offset
        spike_offset += (max_spike + 200.0)
    # Convert the spike dictionary to numpy object array.
    spike_dict = {k: v for k, v in sorted(spike_dict.items(), key=lambda item: item[0])}
    cell_ids = np.array(list(spike_dict.keys()))
    spike_times = np.array(list(spike_dict.values()), dtype=object)
    return spike_times, cell_ids

def compute_correlation_matrix(x: np.ndarray, y: np.ndarray):
    """ Compute the correlation matrix between matrices x and y. """
    n_dim = np.ndim(x)
    if x.shape[n_dim-1] != y.shape[n_dim-1]:
        raise ValueError('x and y must have the same number of columns.')
    if n_dim == 1:
        corr = corrcpp.correlation_coefficient(x, y)
    elif n_dim == 2:
        corr = np.zeros((x.shape[0], y.shape[0]))
        for i in tqdm(range(x.shape[0]), desc='Computing correlation matrix'):
            for j in range(y.shape[0]):
                corr[i,j] = corrcpp.correlation_coefficient(x[i,:], y[j,:])
    return corr

def compute_correlogram(
        spike_times_1: np.ndarray, 
        spike_times_2: np.ndarray, 
        bin_width: float=1.0, 
        max_lag: float=50.0,
        normalization: str='probability') -> Tuple[np.ndarray, np.ndarray]:
    """ Compute the cross-correlogram between two spike trains. 
    Arguments:
        spike_times_1: Spike times for the first cell. [num spikes]
        spike_times_2: Spike times for the second cell. object array: [num cells] with each element being an array of spike times for that cell. [num spikes]
        bin_width: Bin width for the CCG in msec (default: 1.0).
        max_lag: Maximum lag for the CCG in msec (default: 50.0).
    Returns:
        ccg: Cross-correlogram between the two spike trains.
        lags: Lags for the cross-correlogram (in msec).
    """
    num_bins = np.round(max_lag/bin_width).astype(int)
    # Compute the lags in msec.
    lags = np.linspace(-max_lag,max_lag, 2*num_bins+1)
    st1 = spike_times_1.astype(int)
    st1 = np.sort(st1)
    if spike_times_2.dtype == np.dtype('O'): # Object array.
        ccg = np.zeros((len(spike_times_2), 2*num_bins+1))
        for cell_count in range(len(spike_times_2)): #cell_count in tqdm(range(len(spike_times_2)), desc='Computing CCG'):
            st2 = spike_times_2[cell_count].astype(int)
            st2 = np.sort(st2)
            y = np.array(corrcpp.crosscorrelogram(st1, st2, bin_width, max_lag, normalization))
            ccg[cell_count,:] = y
    else:
        ccg = np.array(corrcpp.crosscorrelogram(st1, spike_times_2.astype(int), bin_width, max_lag, normalization))
    ccg = np.nan_to_num(ccg, nan=0.0, posinf=0.0, neginf=0.0) # Convert NaNs to zeros.
    return ccg, lags

def compute_autocorrelations(
        spike_times: np.ndarray, 
        bin_width: float=1.0, 
        max_lag: float=50.0,
        positive_lags_only: bool=True) -> Tuple[np.ndarray, np.ndarray]:
    """ Compute the autocorrelogram for a single spike train. 
    Arguments:
        spike_times: Spike times for the cell. Can be either an int array for a single cell ([num spikes]) or an object array for multiple cells ([num cells x num spikes]).
        bin_width: Bin width for the CCG in msec (default: 1.0).
        max_lag: Maximum lag for the CCG in msec (default: 50.0).
    Returns:
        ccg: Autocorrelogram for the spike train.
        lags: Lags for the autocorrelogram (in msec).
    """
    num_bins = np.round(max_lag/bin_width).astype(int)
    # Compute the lags in msec.
    if positive_lags_only:
        lags = np.linspace(0,max_lag, num_bins+1)
    else:
        lags = np.linspace(-max_lag,max_lag, 2*num_bins+1)
    if spike_times.dtype == np.dtype('O'): # Object array.
        if positive_lags_only:
            ccg = np.zeros((len(spike_times), num_bins+1))
        else:
            ccg = np.zeros((len(spike_times), 2*num_bins+1))
        for cell_count in range(len(spike_times)): 
            s_times = spike_times[cell_count].astype(int)
            s_times = np.sort(s_times)
            if positive_lags_only:
                y = np.array(corrcpp.autocorrelogram(s_times, bin_width, max_lag))
            else:
                y = np.array(corrcpp.crosscorrelogram(s_times, s_times, bin_width, max_lag))
            ccg[cell_count,:] = y
    else:
        s_times = spike_times.astype(int)
        s_times = np.sort(s_times)
        if positive_lags_only:
            ccg = np.array(corrcpp.autocorrelogram(s_times, bin_width, max_lag))
        else:
            ccg = np.array(corrcpp.crosscorrelogram(s_times, s_times, bin_width, max_lag))
    ccg = np.nan_to_num(ccg, nan=0.0, posinf=0.0, neginf=0.0) # Convert NaNs to zeros.
    return ccg, lags

def compute_crosscorrelations(spike_times: np.ndarray, cell_ids: np.ndarray, target_ids: Union[int,list]=None, bin_width: float=1.0, max_lag: float=50.0, normalization: str='probability') -> Tuple[np.ndarray, np.ndarray]:
    """ Compute the cross-correlogram between two spike trains. 
    Arguments:
        spike_times: Spike times for the cells. Object array: [num cells] with each element being an array of spike times for that cell. [num spikes]
        cell_ids: Cell IDs for the cells. [num cells]
        target_ids: Cell IDs for the target cells. If None, then compute the CCG for all cells.
        bin_width: Bin width for the CCG in msec (default: 1.0).
        max_lag: Maximum lag for the CCG in msec (default: 50.0).
    Returns:
        ccg: Cross-correlogram between the two spike trains.
        lags: Lags for the cross-correlogram (in msec).
    """
    if target_ids is None:
        target_ids = cell_ids
    if isinstance(target_ids, int):
        target_ids = [target_ids]
    num_bins = np.round(max_lag/bin_width).astype(int)
    ccg = np.zeros((len(target_ids), len(cell_ids), 2*num_bins+1))
    for cell_count, target_id in enumerate(tqdm(target_ids, desc='Computing CCG')):
        if target_id in cell_ids:
            src_idx = np.where(cell_ids == target_id)[0][0]
            source_spikes = spike_times[src_idx]
            if len(source_spikes) > 0:
                ccg[cell_count,:,:], lags = compute_correlogram(spike_times_1=source_spikes, spike_times_2=spike_times, bin_width=bin_width, max_lag=max_lag, normalization=normalization)
    return ccg, lags
