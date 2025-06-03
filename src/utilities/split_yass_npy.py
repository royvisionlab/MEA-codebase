import argparse
import os

import numpy as np

from collections import namedtuple

StartEndPair = namedtuple('StartEndPair', ['start', 'end'])

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Split combined YASS output .npy into individual .npy')
    parser.add_argument('yass_npy', type=str, help='combined YASS output')
    parser.add_argument('join_metadata', type=str, help='joining metadata csv')
    parser.add_argument('output_path', type=str, help='output path, must be folder that exists already')

    args = parser.parse_args()

    print("Loading YASS spike times")
    spike_times_array = np.load(args.yass_npy) # this has shape (n, 2) where n is the number of spikes,
        # first column is is the spike time, second column is the unit id

    # figure out how many output files there are, and what the time splits are
    print("Splitting .npy file")
    path_keys, dset_lengths = [], []
    with open(args.join_metadata, 'r') as metadata_csv:
        for line in metadata_csv.readlines():
            data_line = line.strip('\n').split(',')
            path_keys.append(data_line[0])
            dset_lengths.append(int(data_line[1]))

    start_end_interval = []
    cc = 0
    for i in range(0, len(dset_lengths)):
        start = cc
        end = cc + dset_lengths[i]
        start_end_interval.append(StartEndPair(start, end))
        cc = end
    print(start_end_interval)

    spike_times_only = spike_times_array[:,0]
    for fcount, interval in enumerate(start_end_interval):
        start, end = interval.start, interval.end

        data_file = path_keys[fcount].split('/')[-1]

        after_start = (spike_times_only >= start)
        before_end = (spike_times_only < end)

        in_interval_slice = np.logical_and(after_start, before_end)

        big_array = spike_times_array[in_interval_slice,:]
        big_array[:,0] = big_array[:,0] - start

        npy_write_path = os.path.join(args.output_path, data_file + '.npy')

        np.save(npy_write_path, big_array)
