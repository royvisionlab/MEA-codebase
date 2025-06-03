
import sys, os, platform
import visionloader as vl
import visionwriter as vw
from vision_utils import STAWriter, ParamsWriter
import lnp
import argparse
import numpy as np
import config as cfg

import visionloader as vl
import numpy as np
import numpy.ma as ma
import os
from typing import Tuple, Union
from tqdm import tqdm
import re
import pickle

import os
import numpy as np
import visionloader as vl
from tqdm import tqdm
from typing import Tuple
import numpy.ma as ma
import sta_utils as su
from lnp import zero_sta_low_freq, compute_spatial_map, calculate_time_course_count
from analysis_utils import compute_pairwise_correlation


def get_cell_ids_from_chunk(data_dir: str, file_name: str) -> np.ndarray:
    """ Get the cell ids from a chunk.
    Parameters:
        data_dir: The path to the data directory containing the Vision files.
        file_name: The name of the vision file (excluding the extension).
    Returns:
        cell_ids: The cell ids.
    """
    # Load the Vision table.
    v_table = vl.load_vision_data(data_dir, file_name, include_params=True)
    # Get the cell ids.
    cell_ids = v_table.get_cell_ids()
    return np.asarray(cell_ids)

def compute_wave_correlation(wave_matrix1: np.ndarray, wave_matrix2: np.ndarray) -> np.ndarray:
    """ Compute the correlation between the waveforms.
    Parameters:
        wave_matrix1: The matrix of the waveforms for the first chunk.
        wave_matrix2: The matrix of the waveforms for the second chunk.
    Returns:
        wave_corr_matrix: The correlation matrix of the waveforms.
    """
    wave_corr_matrix = compute_pairwise_correlation(wave_matrix1, wave_matrix2)
    # wave_corr_matrix = np.zeros((wave_matrix1.shape[0], wave_matrix2.shape[0])).astype(np.float32)
    # for i in tqdm(range(wave_matrix1.shape[0]), desc='Computing waveform correlation matrix'):
    #     for j in range(wave_matrix2.shape[0]):
    #         corr = np.corrcoef(wave_matrix1[i,:], wave_matrix2[j,:])[0,1]
    #         wave_corr_matrix[i,j] = corr
    return wave_corr_matrix

def compute_electrode_shifts(electrode_locations: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """ Compute the electrode shifts. 
    Parameters:
        electrode_locations: The locations of the electrodes.
    Returns:
        electrodes_original: The indices of the original electrodes.
        electrode_shifts: The indices of the shifted electrodes.
        xy_shifts: The x-y shifts of the electrodes. """
    # Calculate the array spacing manually.
    dxy = np.linalg.norm(electrode_locations[0,:] - electrode_locations, axis=1)
    dxy = np.sort(dxy)
    array_spacing = dxy[1]
    half_spacing = array_spacing / 2
    y_spacing = np.sqrt(array_spacing**2 - half_spacing**2)
    # Shift the electrode locations.
    xy_shifts = np.array([[0,0], [-array_spacing,0], [array_spacing,0], 
            [half_spacing,y_spacing], [half_spacing,-y_spacing], [-half_spacing,y_spacing], [-half_spacing,-y_spacing]])
    electrodes_original = np.zeros(xy_shifts.shape[0], dtype=object)
    electrode_shifts = np.zeros(xy_shifts.shape[0], dtype=object)
    for ii in range(xy_shifts.shape[0]):
        xy_shift = xy_shifts[ii,:]
        electrode_locations_shift = electrode_locations + xy_shift
        # Now find the indices of the shifted locations.
        idx_1 = list()
        idx_2 = list()
        for i in range(electrode_locations.shape[0]):
            d = np.linalg.norm(electrode_locations[i,:] - electrode_locations_shift, axis=1)
            if np.min(d) < array_spacing/6:
                idx_1.append(i)
                idx_2.append(np.argmin(d))
        idx_1 = np.array(idx_1)
        idx_2 = np.array(idx_2)
        electrodes_original[ii] = idx_1
        electrode_shifts[ii] = idx_2
    return electrodes_original, electrode_shifts, xy_shifts

def compute_ei_corr(ei_1: np.ndarray, ei_2: np.ndarray, mask_invalid: bool=False) -> np.ndarray:
    """ Compute the correlation between the EI maps of the two chunks at all possible shifts by one electrode.
    Parameters:
        ei_1: The EI maps for the first chunk.
        ei_2: The EI maps for the second chunk.
        electrodes_original: The indices of the original electrodes.
        electrode_shifts: The indices of the shifted electrodes.
    Returns:
        corr: The correlation between the EI maps of the two chunks at all possible shifts by one electrode.
    """
    corr = compute_pairwise_correlation(m_1=ei_1, m_2=ei_2, mask_invalid=mask_invalid)
    # corr = np.zeros((len(electrode_shifts),ei_1.shape[0],ei_2.shape[0])).astype(np.float32)
    # for ii in tqdm(range(len(electrode_shifts)), desc='Computing EI correlation (unshifted)'):
    #     corr[ii,...] = compute_pairwise_correlation(m_1=ei_1[:, electrodes_original[ii].astype(int)], m_2=ei_2[:, electrode_shifts[ii].astype(int)], mask_invalid=mask_invalid)
    return corr

def compute_shifted_ei_corr(ei_1: np.ndarray, ei_2: np.ndarray, electrodes_original: np.ndarray, electrode_shifts: np.ndarray, mask_invalid: bool=False) -> np.ndarray:
    """ Compute the correlation between the EI maps of the two chunks at all possible shifts by one electrode.
    Parameters:
        ei_1: The EI maps for the first chunk.
        ei_2: The EI maps for the second chunk.
        electrodes_original: The indices of the original electrodes.
        electrode_shifts: The indices of the shifted electrodes.
    Returns:
        corr: The correlation between the EI maps of the two chunks at all possible shifts by one electrode.
    """
    corr = np.zeros((len(electrode_shifts),ei_1.shape[0],ei_2.shape[0])).astype(np.float32)
    for ii in tqdm(range(len(electrode_shifts)), desc='Computing shifted EI correlation'):
        corr[ii,...] = compute_pairwise_correlation(m_1=ei_1[:, electrodes_original[ii].astype(int)], m_2=ei_2[:, electrode_shifts[ii].astype(int)], mask_invalid=mask_invalid)
    return corr

def compute_ei_correlation(cell_ids: np.ndarray, ei_raw_matrix: np.ndarray) -> np.ndarray:
    # Compute the raw EI correlation matrix.
    ei_corr_matrix = np.zeros((len(cell_ids), len(cell_ids))).astype(np.float32)
    for i in tqdm(range(len(cell_ids)), desc='Computing EI correlation matrix'):
        for j in range(i, len(cell_ids)):
            if i != j:
                corr = np.corrcoef(ei_raw_matrix[i,:], ei_raw_matrix[j,:])[0,1]
                ei_corr_matrix[i,j] = corr
                ei_corr_matrix[j,i] = corr

def center_wave(wave: np.ndarray, wave_len: int=201, wave_center: int=100, offset: int=None) -> np.ndarray:
    """ Center a waveform by padding the front and back with zeros. Places the minimum at the wave_center of a total wave_len length
    Parameters:
        wave: The waveform to be centered.
        wave_len: The length of the centered waveform.
        wave_center: The index of the center of the waveform.
        offset: The index of the minimum of the waveform.
    Returns:
        wave: The centered waveform.
    """
    if offset is None:
        offset = np.argmin(wave)
    pad_front = wave_center - offset
    pad_back = wave_len - len(wave) - pad_front
    if pad_front >= 0 and pad_back >= 0:
        wave = np.pad(wave, (pad_front, pad_back), mode='constant')
    elif pad_front >=0:
        wave = np.pad(wave, (pad_front, 0), mode='constant')
        wave = wave[:pad_back]
    elif pad_back >=0:
        wave = np.pad(wave, (0, pad_back), mode='constant')
        wave = wave[-1*pad_front:]
    else:
        wave = wave[-1*pad_front:pad_back]
    return wave

def compute_chunk_properties(data_dir: str, file_name: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """ Compute the correlations within a chunk.
    Parameters:
        data_dir: The path to the data directory containing the Vision files.
        file_name: The name of the vision file (excluding the extension).
    Returns:
        cell_ids: The cell ids.
        ei_corr_matrix: The correlation matrix of the EI maps.
        ei_sub_corr_matrix: The correlation matrix of the EI subunits.
        wave_corr_matrix: The correlation matrix of the waveforms. 
        electrode_locations: The locations of the electrodes.
    """
    # Load the Vision table.
    v_table = vl.load_vision_data(data_dir, file_name, include_ei=True)
    # Get the electrode map.
    electrode_locations = v_table.get_electrode_map()
    # Get the cell ids.
    cell_ids = v_table.get_cell_ids()
    ei_raw_matrix = list()
    ei_sub_matrix = list() #np.zeros((len(cell_ids), electrode_locations.shape[0]))
    wave_matrix = list()
    valid_idx = list()
    for idx, cell_id in enumerate(tqdm(cell_ids, desc='Computing chunk properties')):
        try:
            ei = v_table.get_ei_for_cell(cell_id).ei
            ei_energy = np.mean(np.power(ei, 2), axis=1)
            ei_raw_matrix.append(ei_energy)
            electrode = np.argmax(ei_energy)
            max_wave = ei[electrode, :] / np.max(np.abs(ei[electrode, :]))
            wave_matrix.append( center_wave(max_wave, wave_center=20, wave_len=ei.shape[1]) )
            highest_two_electrodes = np.argsort(ei_energy)[-2:]
            ei_subset = ei_energy.copy()
            ei_subset[highest_two_electrodes] = np.nan
            # ei_sub_matrix[idx,:] = ei_subset.astype(np.float32)
            ei_sub_matrix.append(ei_subset.astype(np.float32))
            valid_idx.append(idx)
        except Exception as e:
            print('Error: ', e)
    # Keep only the valid cells.
    cell_ids = np.array(cell_ids)
    cell_ids = cell_ids[valid_idx]
    ei_sub_matrix = np.array(ei_sub_matrix)
    # ei_sub_matrix = ei_sub_matrix[valid_idx,:]
    wave_matrix = np.array(wave_matrix)
    ei_raw_matrix = np.array(ei_raw_matrix)
    # Make sure there are no NaNs.
    ei_raw_matrix[np.isnan(ei_raw_matrix)] = 0
    # ei_sub_matrix[np.isnan(ei_sub_matrix)] = 0
    wave_matrix[np.isnan(wave_matrix)] = 0
    return cell_ids, ei_raw_matrix, ei_sub_matrix, wave_matrix, electrode_locations

def get_interspike_intervals(spike_times: np.ndarray, bin_width: float=0.5, max_time: float=300.0) -> np.ndarray:
    bin_edges = np.arange(0,max_time+bin_width,bin_width)
    # Compute the interspike interval
    if len(spike_times) > 1:
        isi_tmp = np.diff(spike_times)
        isi = np.histogram(isi_tmp,bins=bin_edges)[0]/float(len(isi_tmp))
    else:
        isi = np.zeros((len(bin_edges)-1,)).astype(int)
    return isi

def compute_chunk_isi(
        data_dir: str, 
        file_name: str, 
        bin_width: float=0.5, 
        max_time: float=300.0,
        sample_rate: float=20000.0) -> Tuple[np.ndarray, np.ndarray]:
    """ Extract the interspike intervals for the chunk.
    Parameters:
        data_dir: The path to the data directory containing the Vision files.
        file_name: The name of the vision file (excluding the extension).
        bin_width: The width of the bins for the ISI.
        max_time: The maximum time for the ISI.
        sample_rate: The sample rate of the data.
    Returns:
        cell_ids: The cell ids.
        isi: The interspike intervals.
    """
    # Load the Vision table.
    v_table = vl.load_vision_data(data_dir, file_name, include_neurons=True)
    # Get the cell ids.
    cell_ids = v_table.get_cell_ids()
    # spike_matrix = np.zeros(len(cell_ids), dtype=object)
    isi = list()
    for ii, cell_id in enumerate(cell_ids):
        spike_times = v_table.get_spike_times_for_cell(cell_id)/sample_rate*1000.0
        # spike_matrix[ii] = spike_times
        isi.append(get_interspike_intervals(spike_times=spike_times, bin_width=bin_width, max_time=max_time))
    return cell_ids, np.asarray(isi)

def get_stas_from_ids(data_dir: str, file_name: str, cell_ids: np.ndarray) -> np.ndarray:
    """ Extract the STAs for the chunk. 
    Parameters:
        data_dir: The path to the data directory containing the Vision files.
        file_name: The name of the vision file (excluding the extension).
        cell_ids: The cell ids.
    Returns:
        sta_matrix: The matrix of the STAs.
    """
    sta_matrix = list()
    valid_idx = list()
    for ii, cell_id in enumerate(cell_ids):
        try: 
            with vl.STAReader(data_dir, file_name) as sta_reader:
                sta_container = sta_reader.get_sta_for_cell_id(cell_id)
            cell_sta = su.get_sta_tensor(sta_container) # [width, height, RGB, time]
            # Move the last axis to the front.
            cell_sta = np.moveaxis(cell_sta, -1, 0) # [time, width, height, RGB]
            sta_matrix.append(cell_sta)
            valid_idx.append(ii)
        except Exception as e:
            print('Error: ', e)
    return np.asarray(sta_matrix), np.asarray(valid_idx)

def get_chunk_stas(data_dir: str, file_name: str) -> Tuple[np.ndarray, np.ndarray]:
    """ Extract the STAs for the chunk.
    Parameters:
        data_dir: The path to the data directory containing the Vision files.
        file_name: The name of the vision file (excluding the extension).
    Returns:
        cell_ids: The cell ids.
        sta_matrix: The matrix of the STAs.
    """
    try:
        # Load the Vision table.
        v_table = vl.load_vision_data(data_dir, file_name, include_params=True)
        # Get the cell ids.
        cell_ids = v_table.get_cell_ids()
        sta_matrix = list()
        for cell_id in cell_ids:
            with vl.STAReader(data_dir, file_name) as sta_reader:
                sta_container = sta_reader.get_sta_for_cell_id(cell_id)
            cell_sta = su.get_sta_tensor(sta_container) # [width, height, RGB, time]
            # Move the last axis to the front.
            cell_sta = np.moveaxis(cell_sta, -1, 0)
            sta_matrix.append(cell_sta)
        return cell_ids, np.asarray(sta_matrix)
    except Exception as e:
        print('Error: ', e)
        return None, None

def find_best_match(
        data_path: str, 
        experiment_name: str, 
        chunk_name1: str, 
        chunk_name2: str, 
        sorter_name1: str, 
        sorter_name2: str, 
        use_isi: bool=False,
        collapse_time: bool=True,
        compute_ei_shifts: bool=False,
        corr_weights: list=[1.0, 1.0, 1.0, 1.0, 3.0]) -> dict:
    """ Find the best match between the clusters.
    Parameters:
        ei_corr_matrix: The correlation matrix of the EI maps.
        ei_sub_corr_matrix: The correlation matrix of the EI subunits.
        wave_corr_matrix: The correlation matrix of the waveforms.
    Returns:
        best_match: The best match between the clusters.
    """
    data_dir1 = os.path.join(data_path, experiment_name, chunk_name1, sorter_name1)
    data_dir2 = os.path.join(data_path, experiment_name, chunk_name2, sorter_name2)
    cell_ids1, ei_raw_matrix1, ei_sub_matrix1, wave_matrix1, electrode_locations1 = compute_chunk_properties(data_dir=data_dir1, file_name=sorter_name1)
    cell_ids2, ei_raw_matrix2, ei_sub_matrix2, wave_matrix2, _ = compute_chunk_properties(data_dir=data_dir2, file_name=sorter_name2)
    if compute_ei_shifts:
        electrodes_original, electrode_shifts, _ = compute_electrode_shifts(electrode_locations1)
        ei_corr_raw = compute_shifted_ei_corr(ei_raw_matrix1, ei_raw_matrix2, electrodes_original, electrode_shifts)
        ei_sub_raw = compute_shifted_ei_corr(ei_sub_matrix1, ei_sub_matrix2, electrodes_original, electrode_shifts, mask_invalid=True)
    else:
        ei_corr_raw = compute_ei_corr(ei_raw_matrix1, ei_raw_matrix2, mask_invalid=False)
        ei_sub_raw = compute_ei_corr(ei_sub_matrix1, ei_sub_matrix2, mask_invalid=True)
    wave_corr_matrix = compute_wave_correlation(wave_matrix1, wave_matrix2)
    # Get the correlation between the ISI.
    if use_isi:
        _, isi_1 = compute_chunk_isi(data_dir=data_dir1, file_name=sorter_name1, bin_width=2.0)
        _, isi_2 = compute_chunk_isi(data_dir=data_dir2, file_name=sorter_name2, bin_width=2.0)
        isi_corr_matrix = compute_pairwise_correlation(isi_1, isi_2, mask_invalid=True)
    else:
        isi_corr_matrix = np.zeros((len(cell_ids1), len(cell_ids2))).astype(np.float32)
    if os.path.exists(os.path.join(data_dir1, sorter_name1 + '.sta')) and os.path.exists(os.path.join(data_dir2, sorter_name2 + '.sta')):
        sta_corr_matrix = compute_sta_correlation_from_ids(data_dir1=data_dir1, data_dir2=data_dir2, file_name1=sorter_name1, file_name2=sorter_name2, cell_ids1=cell_ids1, cell_ids2=cell_ids2, collapse_time=collapse_time)
    else:
        sta_corr_matrix = np.zeros((len(cell_ids1), len(cell_ids2))).astype(np.float32)
    best_match = dict()
    for ii in range(ei_corr_raw.shape[0]):
        if np.ndim(ei_corr_raw) == 3:
            ei_corr = np.max(ei_corr_raw[:,ii,:], axis=0)
            ei_sub_corr = np.max(ei_sub_raw[:,ii,:], axis=0)
        else:
            ei_corr = ei_corr_raw[ii,:]
            ei_sub_corr = ei_sub_raw[ii,:]
        wave_corr = wave_corr_matrix[ii,:]
        isi_corr = isi_corr_matrix[ii,:]
        sta_corr = sta_corr_matrix[ii,:]
        corr = np.mean([corr_weights[0]*ei_corr, corr_weights[1]*ei_sub_corr, corr_weights[2]*wave_corr, corr_weights[3]*isi_corr, corr_weights[4]*sta_corr], axis=0)
        idx = np.argmax(corr)
        best_match[cell_ids1[ii]] = [cell_ids2[idx]]
    return best_match

def find_merge_candidates(cell_ids_1: np.ndarray, 
        cell_ids_2: np.ndarray,
        ei_corr_matrix: np.ndarray, 
        ei_sub_corr_matrix: np.ndarray, 
        wave_corr_matrix: np.ndarray,
        isi_corr_matrix: np.ndarray = None,
        merge_thresholds: list=[[0.93,0.93,0.9], [0.99,0.91,0.96], [0.92,0.99,0.5], [0.99,0.5,0.99]],
        isi_corr_threshold: float=0.6) -> Tuple[dict, np.ndarray]:
    """ Find the clusters that should be merged. 
    Parameters:
        ei_corr_matrix: The correlation matrix of the EI maps.
        ei_sub_corr_matrix: The correlation matrix of the EI subunits.
        wave_corr_matrix: The correlation matrix of the waveforms.
    Returns:
        merge_candidates: The indices of the clusters that should be merged.
        already_merged: The indices of the clusters that have already been merged.
    """
    # Compute the evaluation matrix.
    evaluations = list()
    for threshs in merge_thresholds:
        tests = list()
        tests.append(ei_corr_matrix > threshs[0])
        tests.append(ei_sub_corr_matrix > threshs[1])
        tests.append(wave_corr_matrix > threshs[2])
        if isi_corr_matrix is not None:
            tests.append(isi_corr_matrix > isi_corr_threshold)
        evaluation = np.array(np.logical_and.reduce(tests, axis=0))
        evaluations.append(evaluation)
    evaluations_t = np.array(np.logical_or.reduce(evaluations, axis=0))
    # Find the merge candidates.
    merge_c = np.nonzero(evaluations_t)[0]
    merge_candidates = dict()
    already_merged = list()
    if len(merge_c) > 0:
        for id in merge_c:
            if id in already_merged:
                continue
            merge_idx = np.nonzero(evaluations_t[id])[0]
            merge_list = list()
            for ii in merge_idx:
                if ii not in already_merged:
                    already_merged.append(ii)
                    merge_list.append(cell_ids_2[ii])
            if len(merge_list) > 0:
                merge_candidates[cell_ids_1[id]] = merge_list
    already_merged = sorted(cell_ids_2[already_merged])
    return merge_candidates, already_merged

def compute_sta_correlation_for_merge(cell_sta: np.ndarray, merge_sta: np.ndarray) -> np.ndarray:
    """ Compute the correlation between the STAs.
    Parameters:
        cell_sta: The STAs for the first chunk.
        merge_sta: The STAs for the second chunk.
    Returns:
        sta_corr: The correlation between the STAs.
    """
    cell_sta = cell_sta[...,(0,2)]
    merge_sta = merge_sta[...,(0,2)]
    sta_corr = compute_pairwise_correlation(cell_sta.reshape((cell_sta.shape[0],-1)), merge_sta.reshape((merge_sta.shape[0],-1)))
    return sta_corr.ravel()

def compute_sta_correlation_from_ids(
    data_dir1: str, 
    data_dir2: str, 
    file_name1: str, 
    file_name2: str, 
    cell_ids1: np.ndarray, 
    cell_ids2: np.ndarray,
    time_pts: int=31,
    collapse_time: bool=True) -> np.ndarray:
    if os.path.exists(os.path.join(data_dir1, file_name1 + '.sta')) and os.path.exists(os.path.join(data_dir2, file_name2 + '.sta')):
        sta_1, valid_idx1 = get_stas_from_ids(data_dir=data_dir1, file_name=file_name1, cell_ids=cell_ids1)
        sta_2, valid_idx2 = get_stas_from_ids(data_dir=data_dir2, file_name=file_name2, cell_ids=cell_ids2)
        if len(valid_idx1) == 0 or len(valid_idx2) == 0:
            return np.zeros((len(cell_ids1), len(cell_ids2))).astype(np.float32)
        assert sta_1.shape[2] == sta_2.shape[2] and sta_1.shape[3] == sta_2.shape[3], 'Spatial STA dimensions do not match!'
        sta_1 = sta_1[...,(0,2)]
        sta_2 = sta_2[...,(0,2)]
        # Reverse the time axis.
        sta_1 = sta_1[:,::-1,:,:,:] - np.median(sta_1)
        sta_2 = sta_2[:,::-1,:,:,:] - np.median(sta_2)
        if collapse_time:
            space_map_1 = np.zeros((sta_1.shape[0], sta_1.shape[2], sta_1.shape[3], sta_1.shape[4]))
            for ii in range(sta_1.shape[0]):
                cell_sta = sta_1[ii,...].copy()
                cell_sta = zero_sta_low_freq(cell_sta, threshold_frequency=2.0)
                time_course = calculate_time_course_count(cell_sta, max_stixels=3)
                space_map_1[ii,...] = compute_spatial_map(cell_sta, time_course)
                if np.mean(time_course[3:7]) < 0:
                    space_map_1[ii,...] *= -1
            space_map_2 = np.zeros((sta_2.shape[0], sta_2.shape[2], sta_2.shape[3], sta_2.shape[4]))
            for ii in range(sta_2.shape[0]):
                cell_sta = sta_2[ii,...].copy()
                cell_sta = zero_sta_low_freq(cell_sta, threshold_frequency=1.0)
                time_course = calculate_time_course_count(cell_sta, max_stixels=3)
                space_map_2[ii,...] = compute_spatial_map(cell_sta, time_course)
                if np.mean(time_course[3:7]) < 0:
                    space_map_2[ii,...] *= -1
            sta_corr_matrix = compute_pairwise_correlation(space_map_1.reshape((space_map_1.shape[0],-1)), space_map_2.reshape((space_map_2.shape[0],-1)))
        else:
            nt_1 = sta_1.shape[1]
            nt_2 = sta_2.shape[1]
            time_pts = np.min([time_pts, nt_1, nt_2])
            sta_1 = sta_1[:,:time_pts,:,:,:]
            sta_2 = sta_2[:,:time_pts,:,:,:]
            sta_corr_matrix = compute_pairwise_correlation(sta_1.reshape((sta_1.shape[0],-1)), sta_2.reshape((sta_2.shape[0],-1)))
    else:
        sta_corr_matrix = np.zeros((len(cell_ids1), len(cell_ids2))).astype(np.float32)
    return sta_corr_matrix

def evaluate_merge_sta(data_dir1: str, data_dir2: str, file_name1: str, file_name2: str, merge_candidates: dict, already_merged: np.ndarray, sta_threshold: float=0.1) -> Tuple[dict, list]:
    """ Evaluate the merge candidates based on the STAs. 
    Parameters:
        data_dir1: The path to the first data directory.
        data_dir2: The path to the second data directory.
        file_name1: The name of the first file.
        file_name2: The name of the second file.
        merge_candidates: The merge candidates.
        already_merged: The indices of the clusters that have already been merged.
        sta_threshold: The threshold for the STA.
    Returns:
        merge_candidates: The indices of the clusters that should be merged.
        already_merged: The indices of the clusters that have already been merged.
    """
    if os.path.exists(os.path.join(data_dir1, file_name1 + '.sta')) and os.path.exists(os.path.join(data_dir2, file_name2 + '.sta')):
        invalid_keys = list()
        for key, merge_ids in tqdm(merge_candidates.items(), 'Evaluating merge candidate STAs'):
            sta_corr = compute_sta_correlation_from_ids(data_dir1, data_dir2, file_name1, file_name2, [key], merge_ids, collapse_time=True)
            # cell_sta = get_stas_from_ids(data_dir=data_dir1, file_name=file_name1, cell_ids=[key])
            # merge_sta = get_stas_from_ids(data_dir=data_dir2, file_name=file_name2, cell_ids=merge_ids)
            # sta_corr = compute_sta_correlation_for_merge(cell_sta, merge_sta)
            # Valid merges are only those above a certain threshold.
            valid_merge = np.where(sta_corr > sta_threshold)[0]
            invalid_merge = np.where(sta_corr <= sta_threshold)[0]
            if len(invalid_merge) > 0:
                invalid_idx = merge_ids[invalid_merge]
                for idx in invalid_idx:
                    already_merged = np.delete(already_merged, np.argwhere(already_merged==idx))
            merge_ids = np.array(merge_ids).astype(int)
            if len(valid_merge) > 0:
                merge_ids = merge_ids[valid_merge]
                merge_candidates[key] = merge_ids
            else:
                invalid_keys.append(key)
        if len(invalid_keys) > 0:
            for key in invalid_keys:
                already_merged = np.delete(already_merged, np.argwhere(already_merged==key))
                del merge_candidates[key]
    return merge_candidates, already_merged

def compute_merges_across_chunks(
        data_path: str, 
        experiment_name: str, 
        chunk_name1: str, 
        chunk_name2: str, 
        sorter_name1: str, 
        sorter_name2: str, 
        use_isi: bool=False,
        isi_corr_threshold: float=0.6, 
        sta_corr_threshold: float=0.1) -> Tuple[dict, list]:
    """ Compute the merges across chunks. 
    Parameters:
        data_path: The path to the data.
        experiment_name: The name of the experiment.
        chunk_name1: The name of the first chunk.
        chunk_name2: The name of the second chunk.
        sorter_name1: The name of the first sorter.
        sorter_name2: The name of the second sorter.
        use_isi: Whether to use the ISI for the merge.
        isi_threshold: The threshold for the ISI.
        sta_threshold: The threshold for the STA.
    Returns:
        merge_candidates: The indices of the clusters that should be merged.
        already_merged: The indices of the clusters that have already been merged.
    """
    cell_ids1, ei_raw_matrix1, ei_sub_matrix1, wave_matrix1, electrode_locations1 = compute_chunk_properties(data_path, experiment_name, chunk_name1, sorter_name1)
    cell_ids2, ei_raw_matrix2, ei_sub_matrix2, wave_matrix2, electrode_locations2 = compute_chunk_properties(data_path, experiment_name, chunk_name2, sorter_name2)
    electrodes_original, electrode_shifts, _ = compute_electrode_shifts(electrode_locations1)
    ei_corr_raw = compute_shifted_ei_corr(ei_raw_matrix1, ei_raw_matrix2, electrodes_original, electrode_shifts)
    ei_sub_raw = compute_shifted_ei_corr(ei_sub_matrix1, ei_sub_matrix2, electrodes_original, electrode_shifts, mask_invalid=True)
    wave_corr_matrix = compute_wave_correlation(wave_matrix1, wave_matrix2)
    # Get the correlation between the ISI.
    if use_isi:
        cell_ids_1, isi_1 = compute_chunk_isi(data_path=data_path,experiment_name=experiment_name,chunk_name=chunk_name1, sorter_name=sorter_name1, bin_width=2.0)
        cell_ids_2, isi_2 = compute_chunk_isi(data_path=data_path,experiment_name=experiment_name,chunk_name=chunk_name2, sorter_name=sorter_name2, bin_width=2.0)
        isi_corr_matrix = compute_pairwise_correlation(isi_1, isi_2, mask_invalid=True)
    else:
        isi_corr_matrix = None
    merge_candidates, already_merged = find_merge_candidates(
        cell_ids1, cell_ids2, ei_corr_matrix=np.max(ei_corr_raw,axis=0), 
        ei_sub_corr_matrix=np.max(ei_sub_raw,axis=0), 
        wave_corr_matrix=wave_corr_matrix,
        isi_corr_matrix=isi_corr_matrix,
        isi_threshold=isi_corr_threshold)
    merge_candidates, already_merged = evaluate_merge_sta(data_path, experiment_name, chunk_name1, chunk_name2, sorter_name1, sorter_name2, merge_candidates, already_merged, sta_threshold)
    if os.path.exists(os.path.join(data_path, experiment_name, chunk_name1, sorter_name1, sorter_name1 + '.sta')) and os.path.exists(os.path.join(data_path, experiment_name, chunk_name2, sorter_name2, sorter_name2 + '.sta')):
        invalid_keys = list()
        for key, merge_ids in tqdm(merge_candidates.items(), 'Evaluating merge candidates'):
            cell_sta = get_stas_from_ids(data_path, experiment_name, chunk_name1, sorter_name1, [key])
            merge_sta = get_stas_from_ids(data_path, experiment_name, chunk_name2, sorter_name2, merge_ids)
            sta_corr = compute_sta_correlation_for_merge(cell_sta, merge_sta)
            # Valid merges are only those above a certain threshold.
            valid_merge = np.where(sta_corr > sta_corr_threshold)[0]
            merge_ids = np.array(merge_ids).astype(int)
            if len(valid_merge) > 0:
                merge_ids = merge_ids[valid_merge]
                merge_candidates[key] = merge_ids
            else:
                invalid_keys.append(key)
        if len(invalid_keys) > 0:
            for key in invalid_keys:
                del merge_candidates[key]
    return merge_candidates, already_merged

def write_classification_txt(data_dir1: str, data_dir2: str, file_name1: str, file_name2: str, merge_candidates: dict):
    """ Write the classification txt file.
    Parameters:
        data_dir1: The path to the first data directory.
        data_dir2: The path to the second data directory.
        file_name1: The name of the first file.
        file_name2: The name of the second file.
        merge_candidates: The merge candidates.
    """
    source_txt = os.path.join(data_dir1, file_name1 + '.classification.txt')
    out_txt = os.path.join(data_dir2, file_name2 + '.classification.txt')
    cell_ids2 = get_cell_ids_from_chunk(data_dir2, file_name2)
    source_types = dict()
    with open(source_txt, 'r') as f:
        lines = f.readlines()
    for line in lines:
        s = line.split(' ')
        c_id = int(s[0])
        c_type = s[-1]
        source_types[c_id] = c_type
    out_types = dict()
    for key, value in merge_candidates.items():
        if key in source_types:
            for v in value:
                out_types[v] = source_types[key]
    # Write the output.
    with open(out_txt, 'w') as f:
        for cell_id in cell_ids2:
            if cell_id in out_types:
                f.write(str(cell_id) + ' ' + out_types[cell_id])
            else:
                f.write(str(cell_id) + '  All/unmatched\n')

def write_best_match_txt(data_dir1: str, data_dir2: str, file_name1: str, file_name2: str, best_matches: dict):
    """ Write the classification txt file.
    Parameters:
        data_dir1: The path to the first data directory.
        data_dir2: The path to the second data directory.
        file_name1: The name of the first file.
        file_name2: The name of the second file.
        best_matches: The best match dictionary.
    """
    source_txt = os.path.join(data_dir1, file_name1 + '.classification.txt')
    out_txt = os.path.join(data_dir2, file_name2 + '.classification.txt')
    source_types = dict()
    with open(source_txt, 'r') as f:
        lines = f.readlines()
    for line in lines:
        s = line.split(' ')
        c_id = int(s[0])
        c_type = s[-1]
        source_types[c_id] = c_type
    # Write the output.
    with open(out_txt, 'w') as f:
        for key, value in best_matches.items():
            match_str = source_types[value]
            f.write(str(key) + ' ' + match_str)

def collapse_sta_time_dimension(cell_sta: np.ndarray) -> np.ndarray:
    cell_sta = zero_sta_low_freq(cell_sta, threshold_frequency=2.0)
    time_course = calculate_time_course_count(cell_sta, max_stixels=3)
    space_map = compute_spatial_map(cell_sta, time_course)
    if np.mean(time_course[3:7]) < 0:
        space_map *= -1
    return space_map

def compute_sta_correlation(
    sta_1: np.ndarray, 
    sta_2: np.ndarray, 
    time_pts: int=31,
    collapse_time: bool=True) -> np.ndarray:
    assert sta_1.shape[2] == sta_2.shape[2] and sta_1.shape[3] == sta_2.shape[3], 'Spatial STA dimensions do not match!'
    sta_1 = sta_1[...,(0,2)]
    sta_2 = sta_2[...,(0,2)]
    # Replace NaNs with zeros.
    sta_1[np.isnan(sta_1)] = 0.0
    sta_2[np.isnan(sta_2)] = 0.0
    # Reverse the time axis.
    sta_1 = sta_1[:,::-1,:,:,:]
    sta_2 = sta_2[:,::-1,:,:,:]
    if collapse_time:
        space_map_1 = np.zeros((sta_1.shape[0], sta_1.shape[2], sta_1.shape[3], sta_1.shape[4]))
        for ii in range(sta_1.shape[0]):
            space_map_1[ii,...] = collapse_sta_time_dimension(sta_1[ii,...])
        space_map_2 = np.zeros((sta_2.shape[0], sta_2.shape[2], sta_2.shape[3], sta_2.shape[4]))
        for ii in range(sta_2.shape[0]):
            space_map_2[ii,...] = collapse_sta_time_dimension(sta_2[ii,...])
        sta_corr_matrix = compute_pairwise_correlation(space_map_1.reshape((space_map_1.shape[0],-1)), space_map_2.reshape((space_map_2.shape[0],-1)))
    else:
        nt_1 = sta_1.shape[1]
        nt_2 = sta_2.shape[1]
        time_pts = np.min([time_pts, nt_1, nt_2])
        sta_1 = sta_1[:,:time_pts,:,:,:]
        sta_2 = sta_2[:,:time_pts,:,:,:]
        sta_corr_matrix = compute_pairwise_correlation(sta_1.reshape((sta_1.shape[0],-1)), sta_2.reshape((sta_2.shape[0],-1)))
    return sta_corr_matrix



def get_electrode_locations(data_dir: str, file_name: str):
    """ Get the electrode locations for the MEA. """
    v_table = vl.load_vision_data(data_dir, file_name, include_ei=True)
    # Get the electrode map.
    electrode_locations = v_table.get_electrode_map()
    return electrode_locations

def create_ei_dict(ei: np.ndarray, spike_counts: np.ndarray, cluster_id: np.ndarray) -> dict:
    """ Create a dictionary of the EI data. 
    Parameters:
        ei: The matrix of the EI maps.
        spike_counts: The spike counts for each cell.
        cluster_id: The cluster ids.
    Returns:
        writeable_ei_by_cell_id: The dictionary of the EI data.
    """
    spike_counts = spike_counts.astype(np.int32)
    writeable_ei_by_cell_id = dict()
    for i in range(ei.shape[0]):
        ei_matrix = ei[i,...]
        writeable_ei = vw.WriteableEIData(ei_matrix, np.zeros(ei_matrix.shape), spike_counts[i])
        writeable_ei_by_cell_id[cluster_id[i]] = writeable_ei
    # Sort by key
    writeable_ei_by_cell_id = dict(sorted(writeable_ei_by_cell_id.items()))
    return writeable_ei_by_cell_id


class Matching(object):
    def __init__(self) -> None:
        pass

    def compute_correlation_between_chunks(self, cell_ids: np.ndarray, ei_raw_matrix: np.ndarray, ei_sub_matrix: np.ndarray, wave_matrix: np.ndarray, 
        cell_ids2: np.ndarray, ei_raw_matrix2: np.ndarray, ei_sub_matrix2: np.ndarray, wave_matrix2: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """ Compute the correlation matrices. 
        Parameters:
            ei_raw_matrix: The matrix of the raw EI maps.
            ei_sub_matrix: The matrix of the EI maps with the two highest electrodes masked.
            wave_matrix: The matrix of the centered waveforms.
        Returns:
            ei_corr_matrix: The correlation matrix of the EI maps.
            ei_sub_corr_matrix: The correlation matrix of the EI subunits.
            wave_corr_matrix: The correlation matrix of the waveforms.
        """
        # Compute the raw EI correlation matrix.
        ei_corr_matrix = np.zeros((len(cell_ids), len(cell_ids2))).astype(np.float32)
        for i in tqdm(range(len(cell_ids)), desc='Computing EI correlation matrix'):
            for j in range(i, len(cell_ids2)):
                corr = np.corrcoef(ei_raw_matrix[i,:], ei_raw_matrix2[j,:])[0,1]
                ei_corr_matrix[i,j] = corr
                # ei_corr_matrix[j,i] = corr

        # Compute subset EI correlation matrix.
        ei_sub_corr_matrix = np.zeros((len(cell_ids), len(cell_ids2))).astype(np.float32)
        for i in tqdm(range(len(cell_ids)), desc='Computing subset EI correlation matrix'):
            for j in range(i, len(cell_ids2)):
                corr = ma.corrcoef(ma.masked_invalid(ei_sub_matrix[i,:]), ma.masked_invalid(ei_sub_matrix2[j,:]))[0,1]
                ei_sub_corr_matrix[i,j] = corr
        # ei_sub_corr_matrix[j,i] = corr

        # Compute the wave correlation matrix.
        wave_corr_matrix = np.zeros((len(cell_ids), len(cell_ids2))).astype(np.float32)
        for i in tqdm(range(len(cell_ids)), desc='Computing waveform correlation matrix'):
            for j in range(i, len(cell_ids2)):
                corr = np.corrcoef(wave_matrix[i,:], wave_matrix2[j,:])[0,1]
                wave_corr_matrix[i,j] = corr
                # wave_corr_matrix[j,i] = corr
        return ei_corr_matrix, ei_sub_corr_matrix, wave_corr_matrix

class Deduplication(object):
    """ Class for deduplicating clusters."""
    def __init__(self, 
                data_path: str, 
                experiment_name: str, 
                sorter_name: str, 
                merge_thresholds: list=None,
                sort_path: str=None, 
                violation_tolerance: float=0.001,
                recompute_rf_fits: bool=True,
                rewrite_sta: bool=False,
                sample_rate: float=20000.0,
                isi_threshold: float=1.5):
        self.data_path = data_path
        self.experiment_name = experiment_name
        self.sorter_name = sorter_name
        self.violation_tolerance = violation_tolerance
        self.recompute_rf_fits = recompute_rf_fits
        self.rewrite_sta = rewrite_sta
        self.electrode_locations = None
        self.sample_rate = sample_rate
        self.isi_threshold = isi_threshold
        self.sta_corr_matrix = None
        self.set_merge_thresholds(merge_thresholds)
        if sort_path is None:
            self.sort_path = data_path
        else:
            self.sort_path = sort_path
        

    def get_electrode_locations(self, data_dir: str, file_name: str):
        """ Get the electrode locations for the MEA. """
        self.electrode_locations = get_electrode_locations(data_dir, file_name)

    def set_merge_thresholds(self, merge_thresholds: list=None):
        """ Set the thresholds for merging clusters. 
        The thresholds are a list of lists, where each sublist contains the thresholds for merging a cluster with the cluster in the corresponding position.
        """
        if merge_thresholds is None:
            self.merge_thresholds = [[0.93,0.93,0.9], [0.99,0.91,0.96], [0.92,0.99,0.5], [0.99,0.5,0.99]]
        else:
            self.merge_thresholds = merge_thresholds
    
    def save_log_file(self, chunk_name: str, merges: dict):
        """ Save a log file of the merges. 
        Parameters:
            chunk_name: The name of the chunk to be analyzed.
            merges: The dictionary of merges.
        """
        file_path = os.path.join(self.sort_path, self.experiment_name, chunk_name + '_' + self.sorter_name + '_merge_log.p')
        with open(file_path, 'wb') as handle:
            pickle.dump({'merges': merges}, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    def process_chunks(self, chunk_name: Union[str,list]):
        """ Process a chunk of data.
        Parameters:
            chunk_name: The name of the chunk to be processed. If a list is provided, then each element of the list is a chunk to be processed.
        """
        if isinstance(chunk_name, str):
            chunk_name = [chunk_name]
        for chunk in chunk_name:
            self.process_chunk(chunk)
    
    def compute_sta_correlations(self, data_dir: str, sorter_name: str, time_pts: int=21):
        """ Compute the STA correlations for the chunk. 
        Parameters:
            data_dir: The path to the data directory.
            sorter_name: The name of the sorter.
        """
        _, sta_matrix = get_chunk_stas(data_dir, sorter_name)
        self.sta_corr_matrix = compute_sta_correlation(sta_1=sta_matrix, sta_2=sta_matrix, time_pts=time_pts, collapse_time=True)
        # Set diagonal values to zero
        np.fill_diagonal(self.sta_corr_matrix, 0.0)
    
    def process_chunk(self, chunk_name: str, online_analysis: bool=False):
        """ Process a chunk of data. 
        Parameters:
            chunk_name: The name of the chunk to be processed.
        """
        # Get the list of files in the chunk.
        if online_analysis:
            chunk_files = [chunk_name]
            file_path = os.path.join(self.data_path, self.experiment_name, chunk_name)
            bottom_level = chunk_name
        else:
            file_path = os.path.join(self.data_path, self.experiment_name, chunk_name, self.sorter_name)
            bottom_level = self.sorter_name
            chunk_file = os.path.join(self.sort_path,self.experiment_name, self.experiment_name + '_' + chunk_name + '.txt')
            assert os.path.isfile(chunk_file), 'Chunk file does not exist in path: {}'.format(chunk_file)
            chunk_files = self.read_chunk_files(chunk_file)

        cell_ids, ei_raw_matrix, ei_sub_matrix, wave_matrix, self.electrode_locations = compute_chunk_properties(data_dir=file_path, file_name=bottom_level)
        ei_corr_matrix, ei_sub_corr_matrix, wave_corr_matrix = self.compute_correlation_matrices(cell_ids, ei_raw_matrix, ei_sub_matrix, wave_matrix)
        # Use the STA if it exists.
        if not online_analysis and os.path.exists(os.path.join(file_path, bottom_level + '.sta')):
            print('Using STA for merge evaluation.')
            self.compute_sta_correlations(data_dir=file_path, sorter_name=bottom_level)
        merge_candidates, already_merged = self.find_merge_candidates(cell_ids=cell_ids, ei_corr_matrix=ei_corr_matrix, ei_sub_corr_matrix=ei_sub_corr_matrix, wave_corr_matrix=wave_corr_matrix)
        
        if merge_candidates is None:
            print('No merge candidates found.')
            return None
        merge_cand, already_merged = self.evaluate_merge_candidates(data_dir=file_path, file_name=bottom_level, merge_candidates=merge_candidates, already_merged=already_merged)
        if merge_cand is None:
            print('All potential merges violated spike refractory tolerance.')
            return None
        # Save a log file of the merges.
        self.save_log_file(chunk_name=chunk_name, merges=merge_cand.copy())
        # First merge the main STA/chunk files.
        merges = self.fill_merge_dictionary(merges=merge_cand.copy(), already_merged=already_merged, cell_ids=cell_ids)
        if online_analysis:
            self.merge_clusters(file_name=chunk_name, merges=merges, use_sta=False, online_analysis=online_analysis)
        else:
            self.merge_clusters(file_name=chunk_name, merges=merges, use_sta=True)

        # Try to merge the individual files in the chunk.
        if not online_analysis:
            for file_name in tqdm(chunk_files, desc='Merging individual files in chunk: ' + chunk_name):
                try:
                    cell_ids = self.get_cell_ids_from_vision(file_name)
                    # Make sure that at least one of the merge keys is in the cell_id list; otherwise, skip the file.
                    if cell_ids is None:
                        print('No cell IDs found for file: {}'.format(file_name))
                    else:
                        if any(f in cell_ids for f in merge_cand.keys()):
                            merges = self.fill_merge_dictionary(merges=merges, already_merged=already_merged, cell_ids=cell_ids)
                            self.merge_clusters(file_name=file_name, merges=merges, use_sta=False)
                        else:
                            print('No merge candidates found for file: {}'.format(file_name))
                except:
                    pass
    
    def get_cell_ids_from_vision(self, file_name: str):
        if ('chunk' in file_name):
            file_path = os.path.join(self.data_path, self.experiment_name, file_name, self.sorter_name)
            bottom_level = self.sorter_name
        else:
            file_path = os.path.join(self.sort_path, self.experiment_name, file_name, self.sorter_name)
            bottom_level = file_name
        try:
            v_table = vl.load_vision_data(file_path, bottom_level, include_neurons=True)
            cell_ids = v_table.get_cell_ids()
        except:
            return None
        return cell_ids
    
    def read_chunk_files(self, chunk_file: str):
        """ Read the data files in the chunk file. 
        Parameters:
            chunk_file: Full path to the chunk txt file.
        Returns:
            file_names: The list of files in the chunk.
        """
        file_names = list()
        with open(chunk_file, 'r') as f:    
            chunk_files = f.readlines()
        chunk_files = chunk_files[0].strip()
        matches = re.findall('data(\d+)',chunk_files)
        if len(matches) > 0:
            for match in matches:
                file_names.append('data' + match)
            return file_names
        else:
            return None

    def center_wave(self, wave: np.ndarray, wave_len: int=201, wave_center: int=100, offset: int=None) -> np.ndarray:
        """ Center a waveform by padding the front and back with zeros. Places the minimum at the wave_center of a total wave_len length
        Parameters:
            wave: The waveform to be centered.
            wave_len: The length of the centered waveform.
            wave_center: The index of the center of the waveform.
            offset: The index of the minimum of the waveform.
        Returns:
            wave: The centered waveform.
        """
        if offset is None:
            offset = np.argmin(wave)
        pad_front = wave_center - offset
        pad_back = wave_len - len(wave) - pad_front

        if pad_front >= 0 and pad_back >= 0:
            wave = np.pad(wave, (pad_front, pad_back), mode='constant')
        elif pad_front >=0:
            wave = np.pad(wave, (pad_front, 0), mode='constant')
            wave = wave[:pad_back]
        elif pad_back >=0:
            wave = np.pad(wave, (0, pad_back), mode='constant')
            wave = wave[-1*pad_front:]
        else:
            wave = wave[-1*pad_front:pad_back]
        return wave

    def find_merge_candidates(self, cell_ids: np.ndarray, ei_corr_matrix: np.ndarray, ei_sub_corr_matrix: np.ndarray, wave_corr_matrix: np.ndarray) -> Tuple[dict, list]:
        """ Find the clusters that should be merged. 
        Parameters:
            ei_corr_matrix: The correlation matrix of the EI maps.
            ei_sub_corr_matrix: The correlation matrix of the EI subunits.
            wave_corr_matrix: The correlation matrix of the waveforms.
        Returns:
            merge_candidates: The indices of the clusters that should be merged.
            already_merged: The indices of the clusters that have already been merged.
        """
        # Compute the evaluation matrix.
        evaluations = list()
        for threshs in self.merge_thresholds:
            tests = list()
            tests.append(ei_corr_matrix > threshs[0])
            tests.append(ei_sub_corr_matrix > threshs[1])
            tests.append(wave_corr_matrix > threshs[2])
            if self.sta_corr_matrix is not None:
                tests.append(self.sta_corr_matrix > 0.1)
            evaluation = np.array(np.logical_and.reduce(tests, axis=0))
            evaluations.append(evaluation)
        evaluations_t = np.array(np.logical_or.reduce(evaluations, axis=0))

        # Find the merge candidates.
        merge_c = np.nonzero(evaluations_t)[0]
        # if len(merge_c) == 0:
        #     return None, None
        merge_candidates = dict()
        already_merged = list()
        if len(merge_c) > 0:
            for id in merge_c:
                if id in already_merged:
                    continue
                merge_idx = np.nonzero(evaluations_t[id])[0]
                merge_list = list()
                merge_list.append(cell_ids[id])
                for i in merge_idx:
                    if i not in already_merged:
                        already_merged.append(i)
                        merge_list.append(cell_ids[i])
                if len(merge_list) > 0:
                    merge_candidates[cell_ids[id]] = merge_list
                    already_merged.append(id)
        already_merged = sorted(cell_ids[already_merged])
        return merge_candidates, already_merged

    def find_merges_across_chunks(self, cell_ids: np.ndarray, cell_ids2: np.ndarray, ei_corr_matrix: np.ndarray, ei_sub_corr_matrix: np.ndarray, wave_corr_matrix: np.ndarray) -> Tuple[dict, list]:
        """ Find the clusters that should be merged. 
        Parameters:
            ei_corr_matrix: The correlation matrix of the EI maps.
            ei_sub_corr_matrix: The correlation matrix of the EI subunits.
            wave_corr_matrix: The correlation matrix of the waveforms.
        Returns:
            merge_candidates: The indices of the clusters that should be merged.
            already_merged: The indices of the clusters that have already been merged.
        """
        # Compute the evaluation matrix.
        evaluations = list()
        for threshs in self.merge_thresholds:
            tests = list()
            tests.append(ei_corr_matrix > threshs[0])
            tests.append(ei_sub_corr_matrix > threshs[1])
            tests.append(wave_corr_matrix > threshs[2])
            evaluation = np.array(np.logical_and.reduce(tests, axis=0))
            evaluations.append(evaluation)
        evaluations_t = np.array(np.logical_or.reduce(evaluations, axis=0))

        # Find the merge candidates.
        merge_c = np.nonzero(evaluations_t)[0]
        # if len(merge_c) == 0:
        #     return None, None
        merge_candidates = dict()
        already_merged = list()
        if len(merge_c) > 0:
            for id in merge_c:
                if id in already_merged:
                    continue
                merge_idx = np.nonzero(evaluations_t[id])[0]
                merge_list = list()
                merge_list.append(cell_ids2[id])
                for i in merge_idx:
                    if i not in already_merged:
                        already_merged.append(i)
                        merge_list.append(cell_ids2[i])
                if len(merge_list) > 0:
                    merge_candidates[cell_ids[id]] = merge_list
                    already_merged.append(id)
        already_merged = sorted(cell_ids2[already_merged])
        return merge_candidates, already_merged

    def evaluate_merge_sta(self, 
            data_dir: str,
            file_name: str, 
            merge_candidates: dict, 
            already_merged: np.ndarray):
        """ Evaluate the merge candidates. 
        Parameters:
            chunk_name: The name of the chunk to be analyzed.
            merge_candidates: The dictionary of merges.
            already_merged: The list of already merged clusters.
            sample_rate: The sample rate of the data.
        Returns:
            merges: The dictionary of merges.
            already_merged: The list of already merged clusters.
        """
        v_table = vl.load_vision_data(data_dir, file_name, include_neurons=True)
        # v_table = vl.load_vision_data(os.path.join(self.data_path, self.experiment_name, chunk_name, self.sorter_name), self.sorter_name, include_neurons=True)
        merges = dict() # The dictionary of merges.
        found_good_merge = False
        for key, merge_ids in merge_candidates.items():
            found_good_merge_for_cell = False
            spike_matrix = np.zeros(len(merge_ids), dtype=object)
            spike_counts = np.zeros(len(merge_ids), dtype=int)
            for merge_num, merge_id in enumerate(merge_ids):
                try:
                    spike_times = v_table.get_spike_times_for_cell(merge_id)/sample_rate*1000.0
                    spike_counts[merge_num] = len(spike_times)
                    spike_matrix[merge_num] = spike_times
                except Exception as e:
                    print('Error: {}'.format(e))
                    spike_counts[merge_num] = 0
                    spike_matrix[merge_num] = []
            # Sort the spike counts in descending order.
            spike_order = np.argsort(spike_counts)[::-1]
            # The cell with the greatest number of spikes becomes the dominant cell.
            dominant_idx = spike_order[0]
            good_merges = list()
            good_merges.append(merge_ids[dominant_idx])
            # Remove the dominant cell from the list of merge candidates.
            spike_order = spike_order[1:]
            for merge_num in spike_order:
                # Check if the cell is a good merge candidate.
                if spike_counts[merge_num] > 0:
                    if self.check_merge_candidate(spike_matrix[dominant_idx], spike_matrix[merge_num]):
                        good_merges.append(merge_ids[merge_num])
                        found_good_merge = True
                        found_good_merge_for_cell = True
                    else: # If it is not a good merge candidate, remove it from the list of merge candidates.
                        already_merged = np.delete(already_merged, np.argwhere(already_merged==merge_ids[merge_num]))
            if found_good_merge_for_cell:
                merges[merge_ids[dominant_idx]] = good_merges
            else:
                already_merged = np.delete(already_merged, np.argwhere(already_merged==merge_ids[dominant_idx]))
        if not found_good_merge:
            return None, None
        return merges, already_merged
    
    def evaluate_merge_candidates(self, 
            data_dir: str,
            file_name: str, 
            merge_candidates: dict, 
            already_merged: np.ndarray, 
            sample_rate: float=20000.0):
        """ Evaluate the merge candidates. 
        Parameters:
            chunk_name: The name of the chunk to be analyzed.
            merge_candidates: The dictionary of merges.
            already_merged: The list of already merged clusters.
            sample_rate: The sample rate of the data.
        Returns:
            merges: The dictionary of merges.
            already_merged: The list of already merged clusters.
        """
        v_table = vl.load_vision_data(data_dir, file_name, include_neurons=True)
        # v_table = vl.load_vision_data(os.path.join(self.data_path, self.experiment_name, chunk_name, self.sorter_name), self.sorter_name, include_neurons=True)
        merges = dict() # The dictionary of merges.
        found_good_merge = False
        for key, merge_ids in merge_candidates.items():
            found_good_merge_for_cell = False
            spike_matrix = np.zeros(len(merge_ids), dtype=object)
            spike_counts = np.zeros(len(merge_ids), dtype=int)
            for merge_num, merge_id in enumerate(merge_ids):
                try:
                    spike_times = v_table.get_spike_times_for_cell(merge_id)/sample_rate*1000.0
                    spike_counts[merge_num] = len(spike_times)
                    spike_matrix[merge_num] = spike_times
                except Exception as e:
                    print('Error: {}'.format(e))
                    spike_counts[merge_num] = 0
                    spike_matrix[merge_num] = []
            # Sort the spike counts in descending order.
            spike_order = np.argsort(spike_counts)[::-1]
            # The cell with the greatest number of spikes becomes the dominant cell.
            dominant_idx = spike_order[0]
            good_merges = list()
            good_merges.append(merge_ids[dominant_idx])
            # Remove the dominant cell from the list of merge candidates.
            spike_order = spike_order[1:]
            for merge_num in spike_order:
                # Check if the cell is a good merge candidate.
                if spike_counts[merge_num] > 0:
                    if self.check_merge_candidate(spike_matrix[dominant_idx], spike_matrix[merge_num]):
                        good_merges.append(merge_ids[merge_num])
                        found_good_merge = True
                        found_good_merge_for_cell = True
                    else: # If it is not a good merge candidate, remove it from the list of merge candidates.
                        already_merged = np.delete(already_merged, np.argwhere(already_merged==merge_ids[merge_num]))
            if found_good_merge_for_cell:
                merges[merge_ids[dominant_idx]] = good_merges
            else:
                already_merged = np.delete(already_merged, np.argwhere(already_merged==merge_ids[dominant_idx]))
        if not found_good_merge:
            return None, None
        return merges, already_merged

    def check_merge_candidate(self, dominant_spikes: np.ndarray, merge_spikes: np.ndarray) -> bool:
        """ Check if the merge candidate is a good merge. 
        Parameters:
            dominant_spikes: The spike times of the dominant cell.
            merge_spikes: The spike times of the merge candidate.
        Returns:
            True if the merge candidate is a good merge, False otherwise.
        """
        min_isi = np.ceil(self.isi_threshold / 1000.0 * self.sample_rate).astype(float)
        # Get the isi violations for the dominant cell.
        good_merge = False
        combined_spikes = np.concatenate((dominant_spikes, merge_spikes))
        combined_spikes = np.sort(combined_spikes)
        violations_combined = np.nansum(self.get_interspike_intervals(combined_spikes, bin_width=0.5, max_time=self.isi_threshold))
        good_merge = (violations_combined < self.violation_tolerance)
        if good_merge:
            return good_merge
        else:
            min_diff = np.array([np.min(np.abs(dominant_spikes.astype(float) - spike.astype(float))) for spike in merge_spikes])
            good_merge = np.any(min_diff > min_isi)
        return good_merge
    
    def merge_spike_trains(self, dominant_spikes: np.ndarray, merge_spikes: np.ndarray) -> np.ndarray:
        """ Merge the spike trains of the dominant cell and the merge candidate. 
        Parameters:
            dominant_spikes: The spike times of the dominant cell. 
            merge_spikes: The spike times of the merge candidate.
        Returns:
            The merged spike train.
        """
        min_isi = np.ceil(self.isi_threshold / 1000.0 * self.sample_rate).astype(float)
        min_diff = np.array([np.min(np.abs(dominant_spikes.astype(float) - spike.astype(float))) for spike in merge_spikes])
        if not np.any(min_diff <= min_isi):
            return np.sort(np.concatenate((dominant_spikes, merge_spikes)))
        good_idx = np.where(min_diff > min_isi)[0]
        return np.sort(np.concatenate((dominant_spikes, merge_spikes[good_idx])))

    def fill_merge_dictionary(self, merges: dict, already_merged: np.ndarray, cell_ids: np.ndarray) -> dict:
        """ Fill in the merge dictionary. 
        Parameters:
            merges: The dictionary of merges.
            already_merged: The list of already merged clusters.
            cell_ids: The list of cell ids.
        Returns:
            merges: The dictionary of merges.
        """
        # Check whether the merge dictionary entries are in cell_ids.
        merges = dict((k, merges[k]) for k in cell_ids if k in merges)
        # Fill in the rest of the merger dictionary.
        for id in cell_ids:
            if id not in already_merged:
                merges[id] = [id]
        # Sort the dictionary by the keys.
        merges = dict(sorted(merges.items()))
        return merges
    
    def compute_correlation_matrices(self, cell_ids: np.ndarray, ei_raw_matrix: np.ndarray, ei_sub_matrix: np.ndarray, wave_matrix: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """ Compute the correlation matrices. 
        Parameters:
            ei_raw_matrix: The matrix of the raw EI maps.
            ei_sub_matrix: The matrix of the EI maps with the two highest electrodes masked.
            wave_matrix: The matrix of the centered waveforms.
        Returns:
            ei_corr_matrix: The correlation matrix of the EI maps.
            ei_sub_corr_matrix: The correlation matrix of the EI subunits.
            wave_corr_matrix: The correlation matrix of the waveforms.
        """
        # Compute the raw EI correlation matrix.
        ei_corr_matrix = np.zeros((len(cell_ids), len(cell_ids))).astype(np.float32)
        for i in tqdm(range(len(cell_ids)), desc='Computing EI correlation matrix'):
            for j in range(i, len(cell_ids)):
                if i != j:
                    corr = np.corrcoef(ei_raw_matrix[i,:], ei_raw_matrix[j,:])[0,1]
                    ei_corr_matrix[i,j] = corr
                    ei_corr_matrix[j,i] = corr

        # Compute subset EI correlation matrix.
        ei_sub_corr_matrix = np.zeros((len(cell_ids), len(cell_ids))).astype(np.float32)
        for i in tqdm(range(len(cell_ids)), desc='Computing subset EI correlation matrix'):
            for j in range(i, len(cell_ids)):
                if i != j:
                    corr = ma.corrcoef(ma.masked_invalid(ei_sub_matrix[i,:]), ma.masked_invalid(ei_sub_matrix[j,:]))[0,1]
                    ei_sub_corr_matrix[i,j] = corr
                    ei_sub_corr_matrix[j,i] = corr

        # Compute the wave correlation matrix.
        wave_corr_matrix = np.zeros((len(cell_ids), len(cell_ids))).astype(np.float32)
        for i in tqdm(range(len(cell_ids)), desc='Computing waveform correlation matrix'):
            for j in range(i, len(cell_ids)):
                if i != j:
                    corr = np.corrcoef(wave_matrix[i,:], wave_matrix[j,:])[0,1]
                    wave_corr_matrix[i,j] = corr
                    wave_corr_matrix[j,i] = corr
        return ei_corr_matrix, ei_sub_corr_matrix, wave_corr_matrix

    def merge_clusters(self, file_name: str, merges: dict, use_sta: bool = False, online_analysis: bool=False):
        """ Merge the clusters in the file. 
        Parameters:
            file_name: The name of the file to be merged. This method treats chunk files and regular data files differently.
            merges: The dictionary of merges.
            use_sta: Whether to merge the STA data (applies only to chunk files).    
        """
        if online_analysis:
            file_path = os.path.join(self.sort_path, self.experiment_name, file_name)
            bottom_level = file_name
        elif ('chunk' in file_name) or use_sta:
            file_path = os.path.join(self.data_path, self.experiment_name, file_name, self.sorter_name)
            bottom_level = self.sorter_name
        else:
            file_path = os.path.join(self.sort_path, self.experiment_name, file_name, self.sorter_name)
            bottom_level = file_name
        
        if not self.rewrite_sta:
            use_sta = False

        # Load the vision table.
        if use_sta:
            v_table = vl.load_vision_data(file_path, bottom_level, include_params=True, include_ei=True, include_neurons=True, include_sta=True, include_runtimemovie_params=True)
        else:
            v_table = vl.load_vision_data(file_path, bottom_level, include_ei=True, include_neurons=True)
        # Get the array id.
        array_id = vl.EIReader(file_path, bottom_level).array_id

        spikes_by_cell_id = dict() # Spike time dictionary for the Neurons file.
        spike_counts = np.zeros((len(merges),1))
        cell_count = 0
        for i, (cell_id, merge_list) in tqdm(enumerate(merges.items()), desc='Merging clusters for {}'.format(file_name)):
            spike_count = 0
            try:
                for merge_count, merge_id in enumerate(merge_list):
                    if merge_count == 0:
                        spike_times = v_table.get_spike_times_for_cell(merge_id)
                    else:
                        spike_times = self.merge_spike_trains(spike_times, v_table.get_spike_times_for_cell(merge_id))
                        # spike_times = np.concatenate((spike_times, v_table.get_spike_times_for_cell(merge_id)))
                    ei_tmp = v_table.get_ei_for_cell(merge_id)
                    if use_sta:
                        try:
                            sta_tmp = v_table.get_sta_for_cell(merge_id)
                            acf_tmp = v_table.get_acf_numpairs_for_cell(merge_id)
                        except:
                            pass
                    if cell_count == 0:
                        ei = np.zeros((len(merges),ei_tmp.ei.shape[0],ei_tmp.ei.shape[1]))
                        if use_sta:
                            sta = np.zeros((len(merges),sta_tmp.red.shape[2],sta_tmp.red.shape[0],sta_tmp.red.shape[1],3))
                            isi = np.zeros((len(merges),acf_tmp.shape[0]))
                    ei[i,:,:] += ei_tmp.ei * ei_tmp.n_spikes
                    if use_sta:
                        try: 
                            isi[i,:] += acf_tmp * ei_tmp.n_spikes
                            sta[i,:,:,:,0] += np.transpose(sta_tmp.red, (2,0,1)) * ei_tmp.n_spikes
                            sta[i,:,:,:,1] += np.transpose(sta_tmp.green, (2,0,1)) * ei_tmp.n_spikes
                            sta[i,:,:,:,2] += np.transpose(sta_tmp.blue, (2,0,1)) * ei_tmp.n_spikes
                        except:
                            pass
                    spike_count += ei_tmp.n_spikes
                    cell_count += 1
                spikes_by_cell_id[cell_id] = np.sort(spike_times)
            except:
                pass
            if use_sta:
                if np.max(np.abs(sta[i,...])) > 0:
                    sta[i,...] /= np.max(np.abs(sta[i,...]))
            if spike_count > 0:
                ei[i,:,:] /= spike_count
                spike_counts[i] = spike_count
                if use_sta:
                    isi[i,:] /= spike_count
            else:
                spike_counts[i] = 0

        spike_counts = spike_counts.astype(np.float32).ravel()
        ei = ei.astype(np.float32)
        cluster_id = np.array(list(merges.keys())).astype(int)
        left_samples = ei_tmp.nl_points
        right_samples = ei_tmp.nr_points
        if use_sta:
            isi = isi.astype(np.float32)
            sta = sta.astype(np.float32)
            sta = sta[:,::-1,:,:,:] # Flip the STA time axis
            sta[:,-1,:,:,:] = 0.0 # Set the last frame to zero 
            movie_params = v_table.get_runtimemovie_params()
            self.write_sta_file(file_path=file_path, sta=sta, cluster_id=cluster_id, micronsPerStixel=movie_params.micronsPerStixelX, refreshPeriod=movie_params.refreshPeriod)
            # Get the STA fit for the cells.
            if self.recompute_rf_fits:
                timecourse_matrix, x0, y0, sigma_x, sigma_y, theta = self.compute_rf_parameters(sta.astype(np.float64))
                timecourse_matrix = timecourse_matrix[:,::-1,:] # Flip the STA time axis
            else:
                timecourse_matrix, x0, y0, sigma_x, sigma_y, theta = self.collect_rf_parameters(v_table, cluster_id)
                y0 = movie_params.height - y0
            
            # Write the new Params file.
            self.write_params_file(file_path=file_path, cluster_id=cluster_id, timecourse_matrix=timecourse_matrix, isi=isi, spike_counts=spike_counts, x0=x0, y0=y0, sigma_x=sigma_x, sigma_y=sigma_y, theta=theta, isi_binning=0.5)

        # params = vl.ParametersFileReader(file_path,self.sorter_name)
        # f_names = params.get_all_field_names()
        # params.col_row_to_arbitrary_data
        # idx = np.argwhere(np.array(f_names) == 'isiBinning')[0][0]
        
        # Write the new EI file.
        self.write_ei_file(file_path=file_path, bottom_level=bottom_level, writeable_ei_by_cell_id=self.create_ei_dict(ei=ei, spike_counts=spike_counts, cluster_id=cluster_id), 
                       left_samples=left_samples, right_samples=right_samples, array_id=array_id)
        # Write the new Neurons file.
        self.write_neurons_file(file_path=file_path, bottom_level=bottom_level, spikes_by_cell_id=spikes_by_cell_id)
        # self.write_ei_file(file_name=file_name, bottom_level=bottom_level, writeable_ei_by_cell_id=self.create_ei_dict(ei=ei, spike_counts=spike_counts, cluster_id=cluster_id), 
        #                left_samples=left_samples, right_samples=right_samples, array_id=array_id)
        # # Write the new Neurons file.
        # self.write_neurons_file(file_name=file_name, bottom_level=bottom_level, spikes_by_cell_id=spikes_by_cell_id)
        
    def compute_rf_parameters(self, sta: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """ Compute the RF parameters for the STA from scratch. 
        Parameters:
            sta: The STA to be analyzed.
            
        Returns:
            timecourse_matrix: The timecourse matrix.
            x0: The x coordinate for the RF center.
            y0: The y coordinate for the RF center.
            sigma_x: The standard deviation of the RF along the x-axis.
            sigma_y: The standard deviation of the RF along the y-axis.
            theta: The orientation of the RF (in degrees).
        """
        timecourse_matrix, _, _, hull_parameters, _, _, _ = lnp.compute_spatiotemporal_maps(sta)
        x0 = hull_parameters[:, 0]
        y0 = sta.shape[2] - hull_parameters[:, 1]
        sigma_x = hull_parameters[:, 2]
        sigma_y = hull_parameters[:, 3]
        theta = hull_parameters[:, 4]
        return timecourse_matrix, x0, y0, sigma_x, sigma_y, theta

    def collect_rf_parameters(self, v_table: vl.VisionCellDataTable, cluster_id: Union[np.ndarray,list]) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """ Collect the RF parameters from the Vision table. 
        Parameters:
            v_table: The Vision table.
            cluster_id: The cluster ids to be analyzed.
        Returns:
            timecourse_matrix: The timecourse matrix.
            x0: The x coordinate for the RF center.
            y0: The y coordinate for the RF center.
            sigma_x: The standard deviation of the RF along the x-axis.
            sigma_y: The standard deviation of the RF along the y-axis.
            theta: The orientation of the RF (in degrees).
        """
        x0 = np.zeros(len(cluster_id))
        y0 = np.zeros(len(cluster_id))
        sigma_x = np.zeros(len(cluster_id))
        sigma_y = np.zeros(len(cluster_id))
        theta = np.zeros(len(cluster_id))
        for ii, cell_id in enumerate(cluster_id):
            try:
                sta_fit = v_table.get_stafit_for_cell(cell_id)
                x0[ii] = sta_fit.center_x
                y0[ii] = sta_fit.center_y
                sigma_x[ii] = sta_fit.std_x
                sigma_y[ii] = sta_fit.std_y
                theta[ii] = sta_fit.rot
                tc_red = v_table.get_data_for_cell(cell_id, 'RedTimeCourse')
                tc_green = v_table.get_data_for_cell(cell_id, 'GreenTimeCourse')
                tc_blue = v_table.get_data_for_cell(cell_id, 'BlueTimeCourse')
                if ii == 0:
                    timecourse_matrix = np.zeros((len(cluster_id), len(tc_red), 3))
                timecourse_matrix[ii,:,0] = tc_red
                timecourse_matrix[ii,:,1] = tc_green
                timecourse_matrix[ii,:,2] = tc_blue
            except:
                pass
        return timecourse_matrix, x0, y0, sigma_x, sigma_y, theta

    def write_params_file(self, file_path: str, cluster_id: Union[np.ndarray,list], timecourse_matrix: np.ndarray, isi: np.ndarray, spike_counts: np.ndarray, x0: np.ndarray, y0: np.ndarray, sigma_x: np.ndarray, sigma_y: np.ndarray, theta: np.ndarray, isi_binning: float=0.5):
        """ Write the Params file. 
        Parameters:
            file_path: The path to the file.
            cluster_id: The cluster ids to be analyzed.
            timecourse_matrix: The timecourse matrix.
            isi: The interspike interval.
            spike_counts: The number of spikes for each cluster.
            x0: The x coordinate for the RF center.
            y0: The y coordinate for the RF center.
            sigma_x: The standard deviation of the RF along the x-axis.
            sigma_y: The standard deviation of the RF along the y-axis.
            theta: The orientation of the RF (in degrees).
            isi_binning: The bin width for the ISI histogram.
        """
        wr = ParamsWriter(filepath = os.path.join(file_path, self.sorter_name) + '.params', cluster_id=np.array(cluster_id))
        wr.write(timecourse_matrix=timecourse_matrix, isi=isi, spike_count=spike_counts, x0=x0, y0=y0, sigma_x=sigma_x, sigma_y=sigma_y, theta=theta, isi_binning=isi_binning)
    
    def write_sta_file(self, file_path: str, sta: np.ndarray, cluster_id: Union[np.ndarray,list], micronsPerStixel: float, refreshPeriod: float):
        """ Write the STA file. 
        Parameters:
            file_path: The path to the file.
            sta: The STA to be analyzed.
            cluster_id: The cluster ids to be analyzed.
            micronsPerStixel: The number of microns per stixel.
            refreshPeriod: The refresh period of the movie.
        """
        with STAWriter(filepath = os.path.join(file_path, self.sorter_name) + '.sta') as sw:
            sw.write(sta=sta, cluster_id=np.array(cluster_id), stixel_size=micronsPerStixel, frame_refresh=refreshPeriod)
    
    def write_ei_file(self, file_path: str, bottom_level: str, writeable_ei_by_cell_id: dict, left_samples: int, right_samples: int, array_id: int):
        with vw.EIWriter(analysis_folder_path=file_path, dataset_name=bottom_level, left_samples=left_samples, right_samples=right_samples, array_id=array_id, overwrite_existing=True) as ei_writer:
            ei_writer.write_eis_by_cell_id(writeable_ei_by_cell_id)
    # def write_ei_file(self, file_name: str, bottom_level: str, writeable_ei_by_cell_id: dict, left_samples: int, right_samples: int, array_id: int):
    #     file_path = os.path.join(self.data_path, self.experiment_name, file_name, self.sorter_name)
    #     with vw.EIWriter(analysis_folder_path=file_path, dataset_name=bottom_level, left_samples=left_samples, right_samples=right_samples, array_id=array_id, overwrite_existing=True) as ei_writer:
    #         ei_writer.write_eis_by_cell_id(writeable_ei_by_cell_id)

    def create_ei_dict(self, ei: np.ndarray, spike_counts: np.ndarray, cluster_id: np.ndarray) -> dict:
        spike_counts = spike_counts.astype(np.int32)
        writeable_ei_by_cell_id = dict()
        for i in range(ei.shape[0]):
            ei_matrix = ei[i,...]
            writeable_ei = vw.WriteableEIData(ei_matrix, np.zeros(ei_matrix.shape), spike_counts[i])
            writeable_ei_by_cell_id[cluster_id[i]] = writeable_ei
        # Sort by key
        writeable_ei_by_cell_id = dict(sorted(writeable_ei_by_cell_id.items()))
        return writeable_ei_by_cell_id

    def write_neurons_file(self, file_path: str, bottom_level: str, spikes_by_cell_id: dict):
        with vl.NeuronsReader(analysis_folder_path=file_path, dataset_name=bottom_level) as nr:
            num_samples = nr.n_samples
            ttl_times = nr.get_TTL_times()
        with vw.NeuronsFileWriter(analysis_write_path=file_path, dset_name=bottom_level) as nfw:
            nfw.write_neuron_file(spike_times_by_cell_id=spikes_by_cell_id, ttl_times=ttl_times, n_samples_total=num_samples)
    # def write_neurons_file(self, file_name: str, bottom_level: str, spikes_by_cell_id: dict):
    #     file_path = os.path.join(self.data_path, self.experiment_name, file_name, self.sorter_name)
    #     with vl.NeuronsReader(analysis_folder_path=file_path, dataset_name=bottom_level) as nr:
    #         num_samples = nr.n_samples
    #         ttl_times = nr.get_TTL_times()
    #     with vw.NeuronsFileWriter(analysis_write_path=file_path, dset_name=bottom_level) as nfw:
    #         nfw.write_neuron_file(spike_times_by_cell_id=spikes_by_cell_id, ttl_times=ttl_times, n_samples_total=num_samples)
    
    def get_interspike_intervals(self, spike_times: np.ndarray, bin_width: float=0.5, max_time: float=300.0) -> np.ndarray:
        bin_edges = np.arange(0,max_time+bin_width,bin_width)
        # Compute the interspike interval
        if len(spike_times) > 1:
            isi_tmp = np.diff(spike_times / 1000.0 * self.sample_rate)
            isi = np.histogram(isi_tmp,bins=bin_edges)[0]/float(len(isi_tmp))
        else:
            isi = np.zeros((len(bin_edges)-1,)).astype(int)
        return isi

    def cross_chunk_candidates(self, chunk_name: str, chunk_name2) -> dict:
        cell_ids, ei_raw_matrix, ei_sub_matrix, wave_matrix, self.electrode_locations = self.compute_chunk_properties(chunk_name)
        cell_ids2, ei_raw_matrix2, ei_sub_matrix2, wave_matrix2, self.electrode_locations = self.compute_chunk_properties(chunk_name2)
        ei_corr_matrix, ei_sub_corr_matrix, wave_corr_matrix = self.compute_correlation_between_chunks(cell_ids, ei_raw_matrix, ei_sub_matrix, wave_matrix, cell_ids2, ei_raw_matrix2, ei_sub_matrix2, wave_matrix2)

        evaluations = list()
        for threshs in self.merge_thresholds:
            tests = list()
            tests.append(ei_corr_matrix > threshs[0])
            tests.append(ei_sub_corr_matrix > threshs[1])
            tests.append(wave_corr_matrix > threshs[2])
            evaluation = np.array(np.logical_and.reduce(tests, axis=0))
            evaluations.append(evaluation)
        evaluations_t = np.array(np.logical_or.reduce(evaluations, axis=0))

        # Find the merge candidates.
        merge_c = np.nonzero(evaluations_t)[0]
        # if len(merge_c) == 0:
        #     return None, None
        merge_candidates = dict()
        if len(merge_c) > 0:
            for id in merge_c:
                merge_idx = np.nonzero(evaluations_t[id])[0]
                merge_list = list()
                # merge_list.append(cell_ids[id])
                for i in merge_idx:
                    merge_list.append(cell_ids2[i])
                if len(merge_list) > 0:
                    merge_candidates[cell_ids[id]] = merge_list
        return merge_candidates


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Deduplicate an experiment across a single Preparation/Mount')
    parser.add_argument('experimentName', type=str, help='Experiment name (e.g. 20230328C)')
    parser.add_argument('-a','--algorithm', default='kilosort2', type=str, help='Sorting algorithm used (yass or kilosort2)')
    parser.add_argument('-c','--chunk_names', nargs='+', default=None, type=str, help='Chunks to deduplicate (e.g., chunk1 chunk2 chunk3)')
    parser.add_argument('-r','--rewrite_sta', default=False, action='store_true', help='Rewrite the STA file')
    parser.add_argument('-o','--online_analysis', default=False, action='store_true', help='Is this online analysis')

    args = parser.parse_args()

    if type(args.chunk_names) is str:
        chunk_names = [args.chunk_names]
    else:
        chunk_names = args.chunk_names

    # Get the output file path.
    data_path, _ = cfg.get_data_paths()

    experiment_name = args.experimentName

    sorter_name = args.algorithm # 'kilosort2' or 'yass'
    rewrite_sta = args.rewrite_sta

    self = Deduplication(data_path=data_path, experiment_name=experiment_name, sorter_name=sorter_name, rewrite_sta=rewrite_sta)

    if args.online_analysis:
        self.process_chunk(chunk_name=chunk_names[0], online_analysis=True)
    else:
        self.process_chunks(chunk_name=chunk_names)


# import numpy as np
# from deduplication import Deduplication
# import config as cfg

# experiment_name='20231026C'
# sorter_name='kilosort4'
# chunk_names = ['chunk4b']
# rewrite_sta=True
# data_path, _ = cfg.get_data_paths()
# self = Deduplication(data_path=data_path, experiment_name=experiment_name, sorter_name=sorter_name, rewrite_sta=rewrite_sta)
# self.process_chunks(chunk_name=chunk_names)
