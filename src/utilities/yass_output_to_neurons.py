
import visionwriter as vw
from visionwriter import WriteableEIData
import bin2py

import argparse
import os
from progressbar import *

import numpy as np

import lib.spsort_utils as s_utils
import lib.harray_utils as harray_utils

from typing import List, Dict, Union

from lib import kilosort_driver_utils, spsort_utils

def compute_ei(bin_file, spike_times, left_samples=20, right_samples=40, n_channels=512, dtype=np.dtype('int16')):
    '''
    read waveforms from recording
    '''

    n_times = left_samples + right_samples + 1

    # read all channels
    channels = np.arange(n_channels)

    # ***** LOAD RAW RECORDING *****
    ei_matrix = np.zeros((n_times, len(channels)),'float32')

    skipped_idx = []
    total_size = n_times*n_channels
    # spike_times are the centers of waveforms
    spike_times_shifted = spike_times - left_samples
    offsets = spike_times_shifted.astype('int64') * dtype.itemsize * n_channels
    spike_count = 0
    with open(bin_file, "rb") as fin:
        for ctr, spike in enumerate(spike_times_shifted):
            try:
                fin.seek(offsets[ctr], os.SEEK_SET)
                wf = np.fromfile(fin,
                                    dtype=dtype,
                                    count=total_size)
                ei_matrix += wf.reshape(
                    n_times, n_channels)[:,channels]
                spike_count += 1
            except:
                skipped_idx.append(ctr)
    if spike_count > 1:
        ei_matrix /= spike_count
    fin.close()

    return WriteableEIData(ei_matrix, np.zeros(ei_matrix.shape), spike_count)

def get_litke_triggers(bin_path, RW_BLOCKSIZE=2000000, TTL_THRESHOLD=-1000):
    epoch_starts = []
    epoch_ends = []
    with bin2py.PyBinFileReader(bin_path, chunk_samples=RW_BLOCKSIZE, is_row_major=True) as pbfr:
        array_id = pbfr.header.array_id
        n_samples = pbfr.length
        for start_idx in range(0, n_samples, RW_BLOCKSIZE):
            n_samples_to_get = min(RW_BLOCKSIZE, n_samples - start_idx)
            samples = pbfr.get_data_for_electrode(0, start_idx, n_samples_to_get)
            # Find the threshold crossings at the beginning and end of each epoch.
            below_threshold = (samples < TTL_THRESHOLD)
            above_threshold = np.logical_not(below_threshold)
            # Epoch starts.
            above_to_below_threshold = np.logical_and.reduce([
                above_threshold[:-1],
                below_threshold[1:]
            ])
            trigger_indices = np.argwhere(above_to_below_threshold) + start_idx
            epoch_starts.append(trigger_indices[:, 0])
            below_to_above_threshold = np.logical_and.reduce([
                below_threshold[:-1],
                above_threshold[1:]
            ])
            trigger_indices = np.argwhere(below_to_above_threshold) + start_idx
            epoch_ends.append(trigger_indices[:, 0])
    epoch_starts = np.concatenate(epoch_starts, axis=0)
    epoch_ends = np.concatenate(epoch_ends, axis=0)
    return epoch_starts, epoch_ends, array_id, n_samples

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Convert Kilosort output to Vision .neurons format')
    parser.add_argument('yass_npy', type=str, help='path to yass output .npy, must already exist')
    parser.add_argument('yass_toplevel', type=str, help='path to YASS toplevel temporary folder')
    parser.add_argument('vision_path', type=str, help='folder to write Vision files, must already exist')
    parser.add_argument('vision_dset_name', type=str, help='name of dataset that we want to write, i.e. data000')
    parser.add_argument('-d', '--datapath', nargs='?', type=str, help='Path to raw data for computing triggers.')

    args = parser.parse_args()


    # deal with TTL and n_samples here, since we need to reach into the raw data to get it
    # also get the electrode map figured out and generate the .globals file here
    n_samples = -1
    ttl_times = np.array([1, ])  # temporary placeholder

    print("Loading YASS spike times")
    spike_times_array = np.load(args.yass_npy)  # this has shape (n, 2) where n is the number of spikes,
    neuron_time_offset = 0

    if args.datapath == None:
        # first column is is the spike time, second column is the unit id
        n_samples = np.max(spike_times_array[:, 0])
        array_id = 503 # Hard code for now.
    else:
        ttl_times, _, array_id, n_samples = get_litke_triggers(args.datapath, RW_BLOCKSIZE=2000000, TTL_THRESHOLD=-1000)

    # print("Finding TTL times")
    # trigger_file_path = os.path.join(args.yass_toplevel, kilosort_driver_utils.TEMPORARY_TRIGGER_TIMES)
    # array_id, n_samples, ttl_times, neuron_time_offset = spsort_utils.load_trigger_pickle_file(trigger_file_path)

    print("Writing globals file for array id {0}".format(array_id))
    with vw.GlobalsFileWriter(args.vision_path, args.vision_dset_name) as gfw:
        gfw.write_simplified_litke_array_globals_file(array_id,
                                                        0,
                                                        0,
                                                        'YASS converted',
                                                        '',
                                                        0,
                                                        n_samples)

    n_spikes = spike_times_array.shape[0]
    spike_times_by_cell_id = {}  # type: Dict[int, List[int]]
    for i in range(n_spikes):

        spike_time, unit_id = spike_times_array[i, 0], spike_times_array[i, 1]

        if unit_id not in spike_times_by_cell_id:
            spike_times_by_cell_id[unit_id] = []
        spike_times_by_cell_id[unit_id].append(spike_time)
    spikes_by_cell_id_np = {}  # type: Dict[int, np.ndarray]

    #
    widgets = ['Computing spike times: ', Percentage(), ' ', Bar(marker='*',
               left='[',right=']'),' ', ETA()]
    pbar = ProgressBar(widgets=widgets, maxval=len(spike_times_by_cell_id))
    pbar.start()
    cell_count = 0
    for cell_id, spike_list in spike_times_by_cell_id.items():
        spikes_by_cell_id_np[cell_id + 1] = np.array(spike_list) + neuron_time_offset
        pbar.update(cell_count)
        cell_count += 1
        # we add 1 because Vision/MATLAB requires that real cells start at index 1
        # (MATLAB does 1-based indexing)
    pbar.finish()

    print("Writing neurons file")
    spikes_by_cell_id = {}
    print("Found {0} cells".format(len(spikes_by_cell_id_np)))

    with vw.NeuronsFileWriter(args.vision_path, args.vision_dset_name) as nfw:
        nfw.write_neuron_file(spikes_by_cell_id_np, ttl_times, n_samples)

    # Get the electrical images.
    bin_path = os.path.join(args.yass_toplevel, args.vision_dset_name + '.bin')
    print(bin_path)

    # Create an EI Writer.
#    left_samples = 20
#    right_samples = 40
#    ei_writer = vw.EIWriter(args.vision_path, args.vision_dset_name, left_samples, right_samples, array_id, overwrite_existing=True)

#    writeable_ei_by_cell_id = dict()
    # Loop through and get the EI.
#    widgets = ['Computing EIs: ', Percentage(), ' ', Bar(marker='*',
#               left='[',right=']'),' ', ETA()]
#    pbar = ProgressBar(widgets=widgets, maxval=len(spike_times_by_cell_id))
#    pbar.start()
#    cell_count = 0
#    for cell_id, spike_list in spike_times_by_cell_id.items():
#        writeable_ei_by_cell_id[cell_id + 1] = compute_ei(bin_path, np.array(spike_list), left_samples, right_samples)
#        pbar.update(cell_count)
#        cell_count += 1
#    pbar.finish()

#    ei_writer.write_eis_by_cell_id(writeable_ei_by_cell_id) 
    print("Done")


# import bin2py
# RW_BLOCKSIZE = 200000
# TTL_THRESHOLD = -1000
# TTL_CHANNEL = 0

# epoch_starts = []
# epoch_ends = []

# bin_path = '/Users/michaelmanookin/Documents/Data/rawdata/MEA_Data/data/raw/20220531C/data026/'

# with bin2py.PyBinFileReader(bin_path, chunk_samples=RW_BLOCKSIZE, is_row_major=True) as pbfr:
#     start_sample = 0
#     # end_sample = pbfr.length if args.end_samples is None else min(pbfr.length, args.end_samples)
#     for start_idx in range(0, pbfr.length, RW_BLOCKSIZE):
#         n_samples_to_get = min(RW_BLOCKSIZE, pbfr.length - start_idx)
#         samples = pbfr.get_data_for_electrode(TTL_CHANNEL, start_idx, n_samples_to_get)
#         # Find the threshold crossings at the beginning and end of each epoch.
#         below_threshold = (samples < TTL_THRESHOLD)
#         above_threshold = np.logical_not(below_threshold)
#         # Epoch starts.
#         above_to_below_threshold = np.logical_and.reduce([
#             above_threshold[:-1],
#             below_threshold[1:]
#         ])
#         trigger_indices = np.argwhere(above_to_below_threshold) + start_idx
#         epoch_starts.append(trigger_indices[:, 0])
#         below_to_above_threshold = np.logical_and.reduce([
#             below_threshold[:-1],
#             above_threshold[1:]
#         ])
#         trigger_indices = np.argwhere(below_to_above_threshold) + start_idx
#         epoch_ends.append(trigger_indices[:, 0])
# epoch_starts = np.concatenate(epoch_starts, axis=0)
# epoch_ends = np.concatenate(epoch_ends, axis=0)