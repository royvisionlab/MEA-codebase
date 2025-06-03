import bin2py
import electrode_map
import electrode_map_503

import tqdm

import numpy as np

import argparse
import os

from typing import List, Dict, Union

RW_BLOCKSIZE = 100000

LITKE_FREQ = 20000

ASSETS_FOLDER = os.path.dirname(os.path.abspath(__file__))


YASS_GEOM_FILE = 'geom.txt'
YASS_CONFIG_FILE = 'config.yaml'
# JOINED_RECORDING_METADATA_CSV = 'yass_join.csv'
MASTER_CONFIG_FILE_PATH=os.path.join(ASSETS_FOLDER, 'assets', 'yass_config_template.yaml')

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='join multiple Litke datasets into single YASS dataset')

    parser.add_argument('binpath', type=str, nargs='+', help='Litke dataset binary path, can chain multiple args')
    parser.add_argument('yass_folder_path', type=str, help='Folder to put the output files in')
    parser.add_argument('ds_name_out', type=str, help='name of output dataset')
    parser.add_argument('-b', '--binary', action='store_true', help='convert and save binary data files')
    parser.add_argument('-m', '--mean', action='store_true', help='Subtract the average across the channels')

    args = parser.parse_args()

    JOINED_RECORDING_METADATA_CSV = args.ds_name_out + '.csv'

    source_data_paths = args.binpath

    # Check if we need to use the jacked up electrode map.
    path_split = args.yass_folder_path.split('-')
    path_split = ''.join(path_split[-1])
    path_split = path_split.split('/')
    if path_split[-1] == '':
        exp_num = int(path_split[-2][:8])
    else:
        exp_num = int(path_split[-1][:8])
    use_jacked_map = (exp_num <= 20230221)

    if args.binary:

        output_file_path = os.path.join(args.yass_folder_path, '{0}.bin'.format(args.ds_name_out))
        with open(output_file_path, 'wb') as yass_binfile_output:

            for binpath in source_data_paths:

                with bin2py.PyBinFileReader(binpath, chunk_samples=RW_BLOCKSIZE, is_row_major=True) as pbfr:

                    n_channels_data = pbfr.num_electrodes
                    n_samples = pbfr.length

                    # now get the electrode map coordinates
                    # this is pretty trivial for the Litke system
                    # there doesn't appear to be the option to ignore disconnected electrodes with yass
                    # so we need to identify them now and cut them out of the data
                    if use_jacked_map: # Jacked up map.
                        print('using jacked electrode map for exp {0}'.format(exp_num))
                        emap_array = electrode_map_503.get_litke_array_coordinates_by_array_id()
                    else:
                        print('using normal electrode map for exp {0}'.format(exp_num))
                        emap_array = electrode_map.get_litke_array_coordinates_by_array_id(pbfr.array_id)

                    connected_electrode_geometry_list_ordered = []
                    included_electrode_id = []
                    disconnected_electrode_set = electrode_map.get_disconnected_electrode_set_by_array_id(pbfr.array_id)
                    print(disconnected_electrode_set)
                    # for i in range(emap_array.shape[0]):
                    #     if i not in disconnected_electrode_set:
                    #         # valid connected electrode
                    #         included_electrode_id.append(i)
                    #         connected_electrode_geometry_list_ordered.append((emap_array[i,0], emap_array[i,1]))
                    # Mike testing something.
                    for i in range(emap_array.shape[0]):
                        included_electrode_id.append(i)
                        connected_electrode_geometry_list_ordered.append((emap_array[i,0], emap_array[i,1]))

                    print("Writing {0} channels, {1} samples of data from {2}".format(n_channels_data, n_samples, binpath))

                    # loop over the data section, dump from bin file to little-endian int16 binary
                    # we assume Intel, Linux, etc. so we don't have to do anything special here
                    print("Converting binary data for {0}...".format(binpath))
                    with tqdm.tqdm(total = n_samples // RW_BLOCKSIZE) as pbar:
                        for start_idx in range(0, n_samples, RW_BLOCKSIZE):

                            n_samples_to_get = min(RW_BLOCKSIZE, n_samples - start_idx)

                            samples_with_ttl = pbfr.get_data(start_idx, n_samples_to_get)
                            samples_without_ttl = samples_with_ttl[1:,:]

                            # Subtract the mean across the channels.
                            # if args.mean:
                            #     samples_without_ttl = samples_without_ttl - samples_without_ttl.mean(axis=0, keepdims=True)

                            samples_connected_channels = samples_without_ttl[included_electrode_id,:]
                            samples_to_write = np.ascontiguousarray(samples_connected_channels.T.astype(np.int16))

                            samples_to_write.tofile(yass_binfile_output)
                            pbar.update(1)


    # now go through the bin files again
    # to write the metadata that we need to undo the join
    # and to generate the YASS configs
    rec_lengths = [] # type: List[int]
    for fcount, binpath in enumerate(source_data_paths):
        # for the first file write all the configs

        with bin2py.PyBinFileReader(binpath, chunk_samples=RW_BLOCKSIZE, is_row_major=True) as pbfr:
            # the yass electrode geometry maps is as follows
            # "x y\n" where x and y are separated by space
            # newline for each electrode
            # in the same order as recorded

            if fcount == 1:

                is_512_board = electrode_map.is_litke_512_board(pbfr.array_id)

                print("Writing geom.txt")
                output_geom_path = os.path.join(args.yass_folder_path, YASS_GEOM_FILE)
                with open(output_geom_path, 'w') as geom_file:
                    for x, y in connected_electrode_geometry_list_ordered:
                        geom_file.write('{0} {1}\n'.format(x, y))

                # now generate the yass config
                print("Writing config yaml")
                output_config_path = os.path.join(args.yass_folder_path, YASS_CONFIG_FILE)
                with open(output_config_path, 'w+') as yass_config_file, \
                        open(MASTER_CONFIG_FILE_PATH, 'r') as master_ks_template:


                    yass_neighbor_pitch = 70 if is_512_board else 35

                    master_template_as_str = master_ks_template.read()
                    yass_config_file.write(master_template_as_str.format('{0}.bin'.format(args.ds_name_out),
                                                                         LITKE_FREQ,
                                                                         len(connected_electrode_geometry_list_ordered),
                                                                         yass_neighbor_pitch))

            # get the length of the recording
            rec_lengths.append(pbfr.length)

    join_metadata_path = os.path.join(args.yass_folder_path, JOINED_RECORDING_METADATA_CSV)
    with open(join_metadata_path, 'w') as join_file:
        for path, length in zip(source_data_paths, rec_lengths):
            join_file.write('{0},{1}\n'.format(path, length))
