import bin2py
import electrode_map

import tqdm

import numpy as np

import lib.spsort_utils as spsort_utils
import lib.kilosort_driver_utils as kilosort_driver_utils

import argparse
import os
import sys # Mike added this

RW_BLOCKSIZE = 100000

TTL_THRESHOLD = 1000

LITKE_FREQ = 20000

# Mike added this
print(sys.path)
# sys.path.append('C:\Users\manoo\Documents\GitRepos\Chichilnisky-Lab\artificial-retina-software-pipeline\utilities')

if __name__ == '__main__':
    '''
    Kilosort does not have a place for the TTL channel, we have to detect TTL events when generating the .neurons file
        from the output of Kilosort
    '''

    parser = argparse.ArgumentParser(description='convert Litke array data into Kilosort/Spyking Circus input format')

    parser.add_argument('binpath', type=str, help='Litke dataset binary path')
    parser.add_argument('ks_folder_path', type=str, help='Folder to put the output files in')
    parser.add_argument('ds_name_out', type=str, help='name of output dataset')
    parser.add_argument('-w', '--writebinary', action='store_true', help='convert and save binary data files')
    parser.add_argument('-k', '--kilosort', action='store_true', help='generate Kilosort mapping file')
    parser.add_argument('-s', '--start_samples', type=int, help='sample number to start at', default=0)
    parser.add_argument('-m', '--end_samples', type=int, help='sample number to end at', default=None)
    parser.add_argument('-d', '--matlabdriverpath', type=str,
                        help='Optional path to MATLAB scripts for running Kilosort. Must specify -k as well')
    parser.add_argument('-r', '--streamingdata', action='store_true', help='Use if running on completed streaming data')
    parser.add_argument('-v', '--version', type=int, default=2, help='Kilosort version, either 2 or 3. Default 2')

    args = parser.parse_args()



    # need a little bit of path wizardry to deal with the possibility of streaming data
    # if the data is streamed data and the path is specified as /Volumes/Stream/2020-08-1-0/data000
    # we need to add a .bin file extension to the end, and use PyBinFileReader in file mode rather than
    # folder mode
    # if not streaming data, do the obvious thing
    binfile_read_path = args.binpath if not args.streamingdata else '{0}.bin'.format(args.binpath)

    with bin2py.PyBinFileReader(binfile_read_path, chunk_samples=RW_BLOCKSIZE) as pbfr:

        array_id = pbfr.array_id

        n_channels_data = pbfr.num_electrodes

        start_sample = args.start_samples
        end_sample = pbfr.length if args.end_samples is None else min(pbfr.length, args.end_samples)

        print("Writing {0} channels, {1} samples of data, starting at sample {2}".format(n_channels_data,
                                                                                         end_sample - start_sample,
                                                                                         start_sample))

        output_file_path = os.path.join(args.ks_folder_path, '{0}.dat'.format(args.ds_name_out))
        print("Converting binary data...")
        if args.writebinary:
            with open(output_file_path, 'wb') as ks_binfile_output, \
                    tqdm.tqdm(total=end_sample // RW_BLOCKSIZE) as pbar:

                # In order to support the ability to run Kilosort for arbitrary
                # sections of the data (i.e. spike sort a chunk of the data starting at sample M
                # and ending at sample N without screwing up the triggers, we loop through
                # the data in two loops:

                # Loop 1 corresponds to the chunk of data before the start sample M. In this loop
                #   all we need to do is keep track of the trigger times
                # Loop 2 corresponds to the chunk of data in the spike-sorting region (between M and N)
                #   here we do everything, keeping track of trigger times and writing raw binary data
                #   for Kilosort to ingest

                ttl_times_buffer = []

                ## LOOP 1 ####################################################################
                for start_idx in range(0, start_sample, RW_BLOCKSIZE):
                    n_samples_to_get = min(RW_BLOCKSIZE, start_sample - start_idx)

                    # just load everything, since otherwise we have to jump around on disk
                    # and loading just the TTL channel doesn't save any time
                    samples_with_ttl = pbfr.get_data(start_idx, n_samples_to_get)

                    print(samples_with_ttl.shape)

                    ttl_samples = samples_with_ttl[:, 0]

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

                ## LOOP 2 #####################################################################
                # loop over the data section, dump from bin file to little-endian int16 binary
                # we assume Intel, Linux, etc. so we don't have to do anything special here
                for start_idx in range(start_sample, end_sample, RW_BLOCKSIZE):
                    ####### write samples #####################################
                    n_samples_to_get = min(RW_BLOCKSIZE, end_sample - start_idx)

                    samples_with_ttl = pbfr.get_data(start_idx, n_samples_to_get)

                    ttl_samples = samples_with_ttl[:, 0]
                    samples_without_ttl = np.ascontiguousarray(samples_with_ttl[:, 1:]).astype(np.int16)

                    samples_without_ttl.tofile(ks_binfile_output)

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

                pbar.close()

                ttl_times = np.concatenate(ttl_times_buffer, axis=0)

                trigger_file_path = os.path.join(args.ks_folder_path, kilosort_driver_utils.TEMPORARY_TRIGGER_TIMES)
                spsort_utils.save_trigger_pickle_file(trigger_file_path, array_id, end_sample,
                                                      ttl_times, neuron_spike_time_offset=start_sample)

    # now generate the electrode map coordinates data structure
    # this is pretty trivial for the Litke system
    emap_array = electrode_map.get_litke_array_coordinates_by_array_id(array_id)
    geometry_dict = {}
    electrode_id_in_order = []
    disconnected_electrode_set = electrode_map.get_disconnected_electrode_set_by_array_id(array_id)
    for i in range(emap_array.shape[0]):
        geometry_dict[i] = (i, i, emap_array[i, 0], emap_array[i, 1])
        electrode_id_in_order.append(i)

    if args.kilosort:

        config_matfile_path = os.path.join(args.ks_folder_path, '{0}_map.mat'.format(args.ds_name_out))

        spsort_utils.generate_kilosort_mapping_file(geometry_dict, electrode_id_in_order,
                                                    os.path.join(args.ks_folder_path, config_matfile_path),
                                                    LITKE_FREQ,
                                                    disconnected_idx=disconnected_electrode_set)

        # generate the MATLAB run scripts for Kilosort if desired
        if args.matlabdriverpath:

            output_file_abs_path = os.path.abspath(args.ks_folder_path)
            config_matfile_abs_path = os.path.abspath(config_matfile_path)

            master_file_path = os.path.abspath(
                os.path.join(args.matlabdriverpath, kilosort_driver_utils.MASTER_KILOSORT_FNAME))
            config_file_path = os.path.abspath(
                os.path.join(args.matlabdriverpath, kilosort_driver_utils.KILOSORT_CONFIG_FNAME))

            if args.version == 2:
                with open(kilosort_driver_utils.MASTER_KILOSORT_TEMPLATE_PATH, 'r') as master_ks_template:
                    master_template_as_str = master_ks_template.read()
            elif args.version == 3:
                with open(kilosort_driver_utils.MASTER_KILOSORT3_TEMPLATE_PATH, 'r') as master_ks_template:
                    master_template_as_str = master_ks_template.read()
            else:
                assert False, "kilosort version must be either 2 or 3"

            with open(master_file_path, 'w+') as ks_driver_file:
                ks_driver_file.write(master_template_as_str.format(output_file_abs_path,
                                                                   config_file_path,
                                                                   config_matfile_abs_path,
                                                                   n_channels_data))

            with open(config_file_path, 'w+') as ks_config_file, \
                    open(kilosort_driver_utils.KILOSORT_CONFIG_TEMPLATE_PATH, 'r') as ks_config_template:

                ks_config_template_as_str = ks_config_template.read()

                variance_distance = spsort_utils.kilosort_variance_distance_array_id(array_id)

                ks_config_file.write(ks_config_template_as_str.format(config_matfile_abs_path,
                                                                      LITKE_FREQ,
                                                                      variance_distance))

    print("Done")
