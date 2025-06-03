
import visionloader as vl
import visionwriter as vw
import os
import argparse
import numpy as np


if __name__ == '__main__':
    '''
    Average and save all of the EI files that were combined for YASS spike sorting.
    Example: 
        python merge_ei.py /usr/share/pool/SortedData/20220712C/ -c chunk1 -file data003 data004
    '''

    parser = argparse.ArgumentParser(description='Average EI files..')

    parser.add_argument('rootdir', type=str, help='Path to sorted spikes and .ei files')
    parser.add_argument('-f','--file', nargs='+', default=None, type=str, help='name of data file(s) to analyze i.e. data003 (default: None)')
    parser.add_argument('-c','--chunk', default=None, type=str, help='chunk name')
    parser.add_argument('-a','--algorithm', default='yass', type=str, help='sorting algorithm (default: yass)')

    args = parser.parse_args()

    rootdir = args.rootdir
    chunk = args.chunk
    algorithm = args.algorithm

    if type(args.file) is str:
        file_names = [args.file]
    else:
        file_names = args.file

    writeable_ei_by_cell_id = dict()
    for count, data_file in enumerate(file_names):
        ei_path = os.path.join(rootdir, data_file, algorithm)
        eir = vl.EIReader(ei_path, data_file)

        electrode_map = eir.get_electrode_map()

        if count == 0:
            array_id = eir.array_id
        eis = eir.get_all_eis_by_cell_id()
        eir.close()

        for ei_cell_id, readable_ei in eis.items():
            ei_matrix = readable_ei.ei
            ei_error = readable_ei.ei_error
            n_spikes = readable_ei.n_spikes

            if count == 0:
                left_samples = readable_ei.nl_points
                right_samples = readable_ei.nr_points
                count += 1
            
            if ei_cell_id in writeable_ei_by_cell_id:
                writeable_ei = writeable_ei_by_cell_id[ei_cell_id]
                ei_matrix = (n_spikes*ei_matrix + writeable_ei.n_spikes*writeable_ei.ei_matrix)/(n_spikes + writeable_ei.n_spikes)
                ei_error = np.sqrt(ei_matrix**2/n_spikes + writeable_ei.ei_matrix**2/writeable_ei.n_spikes)/2.0
                n_spikes = n_spikes + writeable_ei.n_spikes
                
            else:
                writeable_ei = vw.WriteableEIData(ei_matrix, ei_error, n_spikes)
            writeable_ei_by_cell_id[ei_cell_id] = writeable_ei
        print(data_file)

    # Sort by key
    writeable_ei_by_cell_id = dict(sorted(writeable_ei_by_cell_id.items()))

    # Write out the combined EI
    ei_writer = vw.EIWriter(os.path.join(args.rootdir, args.chunk, args.algorithm), args.algorithm, left_samples, right_samples, array_id, True)
    ei_writer.write_eis_by_cell_id(writeable_ei_by_cell_id)