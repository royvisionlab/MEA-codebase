import visionwriter as vw
import bin2py

import numpy as np
import h5py

import os
import argparse

import lib.harray_utils as harray_utils
import lib.kilosort_driver_utils as kilosort_driver_utils
import lib.spsort_utils as spsort_utils
from yass_output_to_neurons import get_litke_triggers

SPIKE_TIMES_FILENAME = 'spike_times.npy'
SPIKE_IDENTITY_FILENAME = 'spike_clusters.npy'
CLUSTER_QUALITY_FILENAME = 'cluster_KSLabel.tsv'


def build_cluster_quality_dict(filepath):
    cluster_quality_by_id = {}
    with open(filepath, 'r') as cluster_quality_file:
        cluster_quality_file.readline()

        remaining_lines = cluster_quality_file.readlines()

        for line in remaining_lines:
            data_list = line.strip('\n').split('\t')
            cluster_quality_by_id[int(data_list[0])+1] = data_list[1]

    return cluster_quality_by_id


def extract_ttl_times (raw_data_path, ttl_threshold, n_samples=None):

    ttl_times = []

    with bin2py.PyBinFileReader(raw_data_path, chunk_samples=10000) as pbfr:

        if n_samples is None:
            n_samples = pbfr.length
        ttl_samples = pbfr.get_data_for_electrode(0, 0, n_samples)
        below_threshold = (ttl_samples < -ttl_threshold)

        j = 0
        while j < ttl_samples.shape[0]:
            while j < ttl_samples.shape[0] and not below_threshold[j]:
                j += 1

            # now we've reached true, or the end
            if j < ttl_samples.shape[0] and below_threshold[j]:
                interval_start = j
                while j < ttl_samples.shape[0] and below_threshold[j]:
                    j += 1

                ttl_times.append(interval_start)

        array_id = pbfr.header.array_id

    return np.array(ttl_times).astype(np.int), n_samples, array_id


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Convert Kilosort output to Vision .neurons format')
    parser.add_argument('kilosort_spike_path', type=str, help='folder containing Kilosort outputs, must already exist')
    parser.add_argument('kilosort_ttl_path', type=str, help='folder containing TTL times pickle, must already exist')
    parser.add_argument('vision_path', type=str, help='path to write Vision files, must already exist')
    parser.add_argument('vision_dset_name', type=str, help='name of dataset that we want to write, i.e. data000')
    parser.add_argument('-l', '--litke', action='store_true', help='Use if original raw data is Litke raw data')
    parser.add_argument('-k', '--hier', action='store_true', help='Use if original raw data is Hierlemann h5')
    parser.add_argument('-q', '--quality', action='store_true', help='Include only clusters that are tagged good by kilosort')
    parser.add_argument('-d', '--datapath', nargs='?', type=str, help='Path to raw data for computing triggers.')

    args = parser.parse_args()

    spike_times_filepath = os.path.join(args.kilosort_spike_path, SPIKE_TIMES_FILENAME)
    spike_identity_filepath = os.path.join(args.kilosort_spike_path, SPIKE_IDENTITY_FILENAME)
    cluster_quality_dict_filepath = os.path.join(args.kilosort_spike_path, CLUSTER_QUALITY_FILENAME)

    spike_times_vector = np.load(spike_times_filepath)
    spike_identity_vector = np.load(spike_identity_filepath)
    quality_by_cluster_id = build_cluster_quality_dict(cluster_quality_dict_filepath)

    n_spikes = spike_times_vector.shape[0]

    # deal with TTL and n_samples here, since we need to reach into the raw data to get it
    # also get the electrode map figured out and generate the .globals file here
    if args.datapath == None:
        if np.ndim(spike_times_vector) == 1:
            n_samples = np.max(spike_times_vector)
        else:
            n_samples = np.max(spike_times_vector[:, 0])
        array_id = 504
        ttl_times = np.array([1, ]) # temporary placeholder
    else:
        ttl_times, _, array_id, n_samples = get_litke_triggers(args.datapath, RW_BLOCKSIZE=2000000, TTL_THRESHOLD=-1000)
    
    neuron_time_offset = 0
    # array_id = 3501
    
    if not args.litke and not args.hier:
        assert False, "Unsupported raw data type specified in arguments"
    elif args.litke:

        trigger_file_path = os.path.join(args.kilosort_ttl_path, args.vision_dset_name+'.p')
        # array_id, n_samples, ttl_times, neuron_time_offset = spsort_utils.load_trigger_pickle_file(trigger_file_path)

        # If the file exists, we need to remove it.
        if os.path.exists(args.vision_path + args.vision_dset_name + '.globals'):
            os.remove(args.vision_path + args.vision_dset_name + '.globals')
        print("Writing globals file for array id {0}".format(array_id))
        with vw.GlobalsFileWriter(args.vision_path, args.vision_dset_name) as gfw:
            gfw.write_simplified_litke_array_globals_file(array_id & 0xFFF, # FIXME get rid of this after we figure out what happened with 120um
                                                          0,
                                                          0,
                                                          'Kilosort converted',
                                                          '',
                                                          0,
                                                          n_samples)

    elif args.hier:
        with h5py.File(args.rawdata_path, 'r') as h5_data:
            n_samples = h5_data[harray_utils.DATA_KEY].shape[1]

            electrode_geometry_dict, ordered_amp_id = harray_utils.generate_mapping_coordinates_from_h5_mapping(h5_data)

            electrode_config_ordered_no_ttl = np.zeros((ordered_amp_id.shape[0], 2))
            for i, amp_id in enumerate(ordered_amp_id):
                electrode_config_ordered_no_ttl[i,0] = electrode_geometry_dict[amp_id][2]
                electrode_config_ordered_no_ttl[i,1] = electrode_geometry_dict[amp_id][3]
            electrode_config_ordered_no_ttl = electrode_config_ordered_no_ttl.astype(np.float)

            print("Writing globals file for Hierlemann array")
            with vw.GlobalsFileWriter(args.vision_path, args.vision_dset_name) as gfw:
                gfw.write_simplified_reconfigurable_array_globals_file(0,
                                                                       0,
                                                                       'Kilosort converted',
                                                                       '',
                                                                       0,
                                                                       harray_utils.HIERLEMANN_FREQ,
                                                                       n_samples,
                                                                       electrode_config_ordered_no_ttl,
                                                                       harray_utils.HIERLEMANN_PITCH)

            trigger_file_path = os.path.join(args.kilosort_ttl_path, args.vision_dset_name+'.p')
            array_id, n_samples, ttl_times, neuron_time_offset = spsort_utils.load_trigger_pickle_file(trigger_file_path)

    print("Writing neurons file")
    spikes_by_cell_id = {}
    for i in range(n_spikes):

        if np.ndim(spike_times_vector) == 1:
            spike_time = spike_times_vector[i]
            spike_id = spike_identity_vector[i] + 1
        else:
            spike_time = spike_times_vector[i,0]
            spike_id = spike_identity_vector[i,0] + 1
        # we add 1 because Vision/MATLAB requires that real cells start at index 1
        # (MATLAB does 1-based indexing)

        if args.quality:
            if quality_by_cluster_id[spike_id] == 'good':
                if spike_id not in spikes_by_cell_id:
                    spikes_by_cell_id[spike_id] = []

                spikes_by_cell_id[spike_id].append(spike_time)
        else:
            if spike_id not in spikes_by_cell_id:
                spikes_by_cell_id[spike_id] = []

            spikes_by_cell_id[spike_id].append(spike_time)

    spikes_by_cell_id_np = {}
    for cell_id, spike_list in spikes_by_cell_id.items():
        spikes_by_cell_id_np[cell_id] = np.array(spike_list) + neuron_time_offset

    print("Found {0} cells".format(len(spikes_by_cell_id_np)))

    # If the file exists, we need to remove it.
    if os.path.exists(args.vision_path + args.vision_dset_name + '.neurons'):
        os.remove(args.vision_path + args.vision_dset_name + '.neurons')
    with vw.NeuronsFileWriter(args.vision_path, args.vision_dset_name) as nfw:
        nfw.write_neuron_file(spikes_by_cell_id_np, ttl_times, n_samples)

    print("Done")
