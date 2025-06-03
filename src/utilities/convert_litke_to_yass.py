import bin2py
import electrode_map
import electrode_map_503

import tqdm

import numpy as np

import argparse
import os
import socket
import warnings

from lib import kilosort_driver_utils, spsort_utils

RW_BLOCKSIZE = 100000
LITKE_FREQ = 20000
TTL_THRESHOLD = 1000

SPIKE_SIZE_MS = 3 # From Peter.
ASSETS_FOLDER = os.path.dirname(os.path.abspath(__file__))
YASS_GEOM_FILE = 'geom.txt'
YASS_CONFIG_FILE = 'config.yaml'
YASS_CONFIG_NN_TRAIN = 'config_train_nn.yaml'
YASS_CONFIG_POST_TRAIN = 'config_post_train.yaml'
CONFIG_FILE_PARENT = os.path.join(ASSETS_FOLDER, 'assets')
SPIKE_TRAIN_DIR = 'tmp/output/spike_train.npy'
NN_DIR = 'tmp/nn_train'

# Modified compute parameters for peggyo.
MASTER_CONFIG_TEMPLATE = 'yass_config_template.yaml'

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='convert Litke array data into YASS input format')

    parser.add_argument('binpath', type=str, help='Litke dataset binary path')
    parser.add_argument('yass_folder_path', type=str, help='Folder to put the output files in')
    parser.add_argument('ds_name_out', type=str, help='name of output dataset')
    parser.add_argument('-w', '--writebinary', action='store_true', help='convert and save binary data files')
    parser.add_argument('-t','--train',action='store_true',help='retrain neural network on data set and re-run yass')
    parser.add_argument('-y', '--yassconfig', type=str, help='specify custom YASS config file')
    parser.add_argument('-s', '--start_samples', type=int, help='sample number to start at', default=0)
    parser.add_argument('-m', '--end_samples', type=int, help='maximum number of samples to convert')

    args = parser.parse_args()

    # Check if we need to use the jacked up electrode map.
    path_split = args.yass_folder_path.split('/')
    if path_split[-1] == '':
        exp_num = int(path_split[-2][:8])
    else:
        exp_num = int(path_split[-1][:8])
    use_jacked_map = (exp_num <= 20230221)

    # Determine whether user provided none, name, or full path to config.
    if args.yassconfig is not None:

        # Use the full path provided by the user, if it exists.
        if os.path.exists(args.yassconfig):
            config_file_path = args.yassconfig

        # Default to assets folder otherwise. 
        else:
            config_file_path = os.path.join(CONFIG_FILE_PARENT,args.yassconfig)
    else:
        config_file_path = os.path.join(CONFIG_FILE_PARENT,
                                        MASTER_CONFIG_TEMPLATE)

    if not os.path.exists(config_file_path):
        raise OSError('Specified YASS config file not found.')

    with bin2py.PyBinFileReader(args.binpath, chunk_samples=RW_BLOCKSIZE, is_row_major=True) as pbfr:

        n_channels_data = pbfr.num_electrodes
        start_sample = args.start_samples
        end_sample = pbfr.length if args.end_samples is None else min(pbfr.length, args.end_samples)

        # now get the electrode map coordinates
        # this is pretty trivial for the Litke system
        # there doesn't appear to be the option to ignore disconnected electrodes with yass
        # so we need to identify them now and cut them out of the data
        array_id = pbfr.array_id
        if use_jacked_map:
            print('using jacked electrode map for exp {0}'.format(exp_num))
            emap_array = electrode_map_503.get_litke_array_coordinates_by_array_id()
        else:
            emap_array = electrode_map.get_litke_array_coordinates_by_array_id(pbfr.array_id)

        yass_neighbor_pitch = spsort_utils.kilosort_variance_distance_array_id(array_id)

        connected_electrode_geometry_list_ordered = []
        included_electrode_id = []
        disconnected_electrode_set = electrode_map.get_disconnected_electrode_set_by_array_id(pbfr.array_id)
        for i in range(emap_array.shape[0]):
            if i not in disconnected_electrode_set:
                # valid connected electrode
                included_electrode_id.append(i)
                connected_electrode_geometry_list_ordered.append((emap_array[i,0], emap_array[i,1]))

        output_file_path = os.path.join(args.yass_folder_path, '{0}.bin'.format(args.ds_name_out))

        if args.writebinary:
            print("Writing {0} channels, {1} samples of data, starting at sample {2}".format(n_channels_data,
                                                                                             end_sample - start_sample,
                                                                                             start_sample))

            print("Converting binary data...")
            with open(output_file_path, 'wb') as yass_binfile_output, \
                    tqdm.tqdm(total = end_sample // RW_BLOCKSIZE) as pbar:

                # In order to support the ability to run YASS for arbitrary
                # sections of the data (i.e. spike sort a chunk of the data starting at sample M
                # and ending at sample N without screwing up the triggers, we loop through
                # the data in two loops:

                # Loop 1 corresponds to the chunk of data before the start sample M. In this loop
                #   all we need to do is keep track of the trigger times
                # Loop 2 corresponds to the chunk of data in the spike-sorting region (between M and N)
                #   here we do everything, keeping track of trigger times and writing raw binary data
                #   for YASS to ingest


                ttl_times_buffer = []

                ## LOOP 1 #####################################################################
                for start_idx in range(0, start_sample, RW_BLOCKSIZE):

                    n_samples_to_get = min(RW_BLOCKSIZE, end_sample - start_idx)

                    samples_with_ttl = pbfr.get_data(start_idx, n_samples_to_get)
                    ttl_samples = samples_with_ttl[0,:]

                    ######## TTL STUFF #######################################
                    # deal with the TTL stuff, figure out where we cross the threshold
                    below_threshold = (ttl_samples < -TTL_THRESHOLD)
                    above_threshold = np.logical_not(below_threshold)

                    below_threshold_to_above_threshold = np.logical_and.reduce([
                        below_threshold[:-1],
                        above_threshold[1:]
                    ])

                    trigger_indices = np.argwhere(below_threshold_to_above_threshold) + start_idx
                    ttl_times_buffer.append(trigger_indices[:, 0])

                    pbar.update(1)

                ## LOOP 2 ######################################################################
                # loop over the data section, dump from bin file to little-endian int16 binary
                # we assume Intel, Linux, etc. so we don't have to do anything special here
                for start_idx in range(start_sample, end_sample, RW_BLOCKSIZE):

                    n_samples_to_get = min(RW_BLOCKSIZE, end_sample - start_idx)

                    samples_with_ttl = pbfr.get_data(start_idx, n_samples_to_get)
                    samples_without_ttl = samples_with_ttl[1:,:]

                    # samples_connected_channels = samples_without_ttl[included_electrode_id,:]
                    samples_connected_channels = samples_without_ttl
                    samples_to_write = np.ascontiguousarray(samples_connected_channels.T.astype(np.int16))

                    #samples_without_ttl.T.tofile(ks_binfile_output)
                    samples_to_write.tofile(yass_binfile_output)

                    ######## TTL STUFF #######################################
                    # deal with the TTL stuff, figure out where we cross the threshold
                    ttl_samples = -samples_with_ttl[0,:]
                    below_threshold = (ttl_samples < TTL_THRESHOLD)
                    above_threshold = np.logical_not(below_threshold)

                    below_threshold_to_above_threshold = np.logical_and.reduce([
                        below_threshold[:-1],
                        above_threshold[1:]
                    ])

                    trigger_indices = np.argwhere(below_threshold_to_above_threshold) + start_idx
                    ttl_times_buffer.append(trigger_indices[:, 0])

                    pbar.update(1)

                ttl_times = np.concatenate(ttl_times_buffer, axis=0)

                trigger_file_path = os.path.join(args.yass_folder_path, kilosort_driver_utils.TEMPORARY_TRIGGER_TIMES)
                spsort_utils.save_trigger_pickle_file(trigger_file_path, array_id, end_sample,
                                                      ttl_times, neuron_spike_time_offset=start_sample)

    # also generate the yass electrode map
    # the yass electrode geometry maps is as follows
    # "x y\n" where x and y are separated by space
    # newline for each electrode
    # in the same order as recorded
    print("Writing geom.txt")
    output_geom_path = os.path.join(args.yass_folder_path, YASS_GEOM_FILE)
    with open(output_geom_path, 'w') as geom_file:
        for x, y in connected_electrode_geometry_list_ordered:
            geom_file.write('{0} {1}\n'.format(x, y))

    # now generate the yass config, for the pretraining.
    print("Writing config yaml")
    output_config_path = os.path.join(args.yass_folder_path, YASS_CONFIG_FILE)
    with open(output_config_path, 'w+') as yass_config_file, \
            open(config_file_path, 'r') as master_ks_template:
        master_template_as_str = master_ks_template.read()
        yass_config_file.write(master_template_as_str.format(
                              '{0}.bin'.format(args.ds_name_out),
                              LITKE_FREQ,
                              len(connected_electrode_geometry_list_ordered),
                              yass_neighbor_pitch,
                              "detect.pt","denoise.pt","",""))

    # If indicated, write the configs for pre/post training.
    if args.train:
        print("Writing NN pre-training config yaml")
        output_config_path = os.path.join(args.yass_folder_path,
                                          YASS_CONFIG_NN_TRAIN)
        with open(output_config_path, 'w+') as yass_config_file,\
                open(config_file_path ,'r') as master_ks_template:
             master_template_as_str = master_ks_template.read()

             # Set the path to the generated spikes file, write to config.
             spike_train_path = os.path.join(args.yass_folder_path,
                                             SPIKE_TRAIN_DIR)

             yass_config_file.write(master_template_as_str.format(
                               '{0}.bin'.format(args.ds_name_out),
                                LITKE_FREQ,
                                len(connected_electrode_geometry_list_ordered),
                                yass_neighbor_pitch,
                                "detect.pt","denoise.pt",
                                spike_train_path,SPIKE_SIZE_MS))

        print("Writing NN post-training config yaml")
        output_config_path = os.path.join(args.yass_folder_path,
                                          YASS_CONFIG_POST_TRAIN)
        with open(output_config_path, 'w+') as yass_config_file,\
                open(config_file_path, 'r') as master_ks_template:
             master_template_as_str = master_ks_template.read()

             # Set the path to the generated NN outputs.
             detect_path = os.path.join(args.yass_folder_path,NN_DIR,
                                        'detect.pt')
             denoise_path = os.path.join(args.yass_folder_path,NN_DIR,
                                         'denoise.pt')
             yass_config_file.write(master_template_as_str.format(
                               '{0}.bin'.format(args.ds_name_out),
                                LITKE_FREQ,
                                len(connected_electrode_geometry_list_ordered),
                                yass_neighbor_pitch,
                                detect_path,denoise_path,
                                "",""))
