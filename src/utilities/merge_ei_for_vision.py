
import visionloader as vl
import visionwriter as vw
import os
import argparse
import numpy as np
from scipy.io import savemat
import pickle

# path1 = '/gscratch/retina/data/sorted/20221006C/noise/kilosort2/'
# path2 = '/gscratch/retina/data/sorted/20221006C/natimg/kilosort2/'

# eir = vl.EIReader(path1, 'noise')
# eis = eir.get_all_eis_by_cell_id()
# eir.close()

# eir = vl.EIReader(path2, 'natimg')
# eis2 = eir.get_all_eis_by_cell_id()
# eir.close()

# def compute_corr_coeff(ei1, ei2):
#     return np.sum(ei1 * ei2) / np.sqrt(np.sum(ei1*ei1) * np.sum(ei2*ei2))

# # Remove the n-strongest electrodes.
# def compute_ei_subset(raw_ei1, raw_ei2, num=2):
#     amps1 = np.sum(raw_ei1, axis=1)
#     amps2 = np.sum(raw_ei2, axis=1)
#     # Indices for two strongest electrodes.
#     ind = np.argpartition(amps1,-num)[-num:]
#     np.append(ind, np.argpartition(amps2,-num)[-num:])
#     ei1 = raw_ei1.copy()
#     ei1 = np.delete(ei1,ind,axis=0)
#     ei2 = raw_ei2.copy()
#     ei2 = np.delete(ei2,ind,axis=0)
#     return np.sum(ei1 * ei2) / np.sqrt(np.sum(ei1*ei1) * np.sum(ei2*ei2)) # Return the r^2

# raw_corr = np.zeros((len(eis), len(eis2)))
# sub_corr = np.zeros((len(eis), len(eis2)))

# sta_count = 0
# for sta_cell_id, readable_ei in eis.items():
#     sta_matrix = readable_ei.ei
#     # Load in the EIs from the other stimulus and compare.
#     ei_count = 0
#     for ei_cell_id, other_ei in eis2.items():
#         ei_matrix = other_ei.ei
#         # Raw correlation
#         raw_corr[sta_count,ei_count] = compute_corr_coeff(sta_matrix, ei_matrix)
#         # Subset correlation.
#         sub_corr[sta_count,ei_count] = compute_ei_subset(sta_matrix, ei_matrix)
#         ei_count += 1
#     sta_count += 1

def map_wn_to_ns_cellids(ns_vcd,ns_cellids,wn_vcd,
                         wn_cellids,celltypes_dict,
                         corr_dict = {'on parasol': .95,'off parasol': .95,
                                       'on midget': .95,'off midget': .95,
                                       'a1': .95},
                         mask=False,n_sigmas=None):
    """
    Maps WN to NS EIs according to a threshold value of correlation. Computes
    EI power over space, both with and without masking (user choice). Does a
    pass over the NS cellids and finds the corresponding WN cell. If none is
    found the cell doesn't get mapped (does not appear in the dictionary).
    Also writes a field for the normalized x,y locations for each RF.
    Parameters:
        ns_vcd: natural scenes vision data object
        ns_cellids: natural scenes cellids to map
        wn_cellids: white noise cellids to map
        celltypes_dict: dictionary mapping white noise cell ids to celltype.
    """

    channel_noise = wn_vcd.channel_noise

    # Initialize a dictionary and loop over the cells.
    cellids_dict = dict()

    for key in ['wn_to_ns','ns_to_wn','celltypes']:
        cellids_dict[key] = dict()

    for wn_cell in wn_cellids:

        # Get the cell type and write as well.
        celltype = celltypes_dict[wn_cell].lower()

        # Hardcode these for now TODO FIXME
        if "on" in celltype and "parasol" in celltype:
            celltype = 'on parasol'
        elif "off" in celltype and "parasol" in celltype:
            celltype = "off parasol"
        elif "on" in celltype and "midget" in celltype:
            celltype = 'on midget'
        elif "off" in celltype and "midget" in celltype:
            celltype = 'off midget'
        elif "a1" in celltype:
            celltype = 'a1'
        else:
            continue

        # If masking, only look at the significant indices. 
        wn_cell_ei = wn_vcd.get_ei_for_cell(wn_cell).ei

        if mask and n_sigmas is not None:
            sig_inds = np.argwhere(np.abs(np.amin(wn_cell_ei,axis=1))
                                   > n_sigmas * channel_noise).flatten()
            wn_cell_ei_power = np.zeros(wn_cell_ei.shape[0])
            wn_cell_ei_power[sig_inds] = np.sum(wn_cell_ei[sig_inds,:]**2,
                                                axis=1)
        else:
            wn_cell_ei_power = np.sum(wn_cell_ei**2,axis=1)

        corrs = []

        for ns_cell in ns_cellids:
            ns_cell_ei = ns_vcd.get_ei_for_cell(ns_cell).ei

            if mask and n_sigmas is not None:
                sig_inds = np.argwhere(np.abs(np.amin(ns_cell_ei,axis=1))
                                       > n_sigmas * channel_noise).flatten()
                ns_cell_ei_power = np.zeros(ns_cell_ei.shape[0])
                ns_cell_ei_power[sig_inds] = np.sum(ns_cell_ei[sig_inds,:]**2,
                                                axis=1)
            else:
                ns_cell_ei_power = np.sum(ns_cell_ei**2,axis=1)

            corr = np.corrcoef(wn_cell_ei_power,ns_cell_ei_power)[0,1]
            corrs.append(corr)

        # Take the cell with the largest correlation.
        if np.max(corrs) < corr_dict[celltype]:
            continue

        max_ind = np.argmax(np.asarray(corrs))
        cellids_dict['wn_to_ns'][wn_cell] = ns_cellids[max_ind]
        cellids_dict['ns_to_wn'][ns_cellids[max_ind]] = wn_cell
        cellids_dict['celltypes'][wn_cell] = celltype

        # Once the cell has been mapped, remove it (hack) FIXME.
        ns_cellids.remove(ns_cellids[max_ind]) 



if __name__ == '__main__':
    '''
    Average and save all of the EI files that were combined for YASS spike sorting.
    Example: 
        python merge_ei_for_vision.py /usr/share/pool/SortedData/20220712C/ noise /media/mike/NVME/sort/20220712C/noise.csv yass
    '''

    parser = argparse.ArgumentParser(description='Loop through the root directory and average all EI files.')

    parser.add_argument('rootdir', type=str, help='Path to sorted spikes and .ei files')
    parser.add_argument('dataset_name', type=str, help='Name you want to give to the ei file')
    parser.add_argument('csv_path', type=str, help='Path to csv metadata file for merged files.')
    parser.add_argument('sortDir', type=str, default='yass', help='yass (default) or kilosort2') 

    args = parser.parse_args()

    path_keys, dset_lengths = [], []
    with open(args.csv_path, 'r') as metadata_csv:
        for line in metadata_csv.readlines():
            data_line = line.strip('\n').split(',')
            path_keys.append(data_line[0])
            dset_lengths.append(int(data_line[1]))


    writeable_ei_by_cell_id = dict()
    for count, value in enumerate(path_keys):
        data_file = value.split('/')[-1]
        ei_path = os.path.join(args.rootdir, data_file, args.sortDir)
        # ei_path = '/usr/share/pool/SortedData/20220712C/' + data_file + '/' + args.sortDir + '/' 
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

            # This is a test to see if the channels are misordered in the file.
            # ei_matrix[1:,:] = ei_matrix[:-1,:]

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
    ei_writer = vw.EIWriter(args.rootdir, args.dataset_name, left_samples, right_samples, array_id, True)
    ei_writer.write_eis_by_cell_id(writeable_ei_by_cell_id)

    # Convert the EI's to a matrix and save to MAT and Pickle
    E = np.zeros((len(writeable_ei_by_cell_id), ei_matrix.shape[0], ei_matrix.shape[1]))
    count = np.zeros((len(writeable_ei_by_cell_id),1))
    ct=0
    for key in writeable_ei_by_cell_id:
        myei = writeable_ei_by_cell_id.get(key)
        E[ct,...] = myei.ei_matrix
        count[ct] = myei.n_spikes
        ct += 1

    # MAT file
    mdic = {'E': E, 'count': count, 'electrode_map': electrode_map}
    savemat(os.path.join(args.rootdir, args.dataset_name + '_' + args.sortDir + '_ei.mat'), mdic)

    # Pickle file
    with open(os.path.join(args.rootdir, args.dataset_name + '_' + args.sortDir + '_ei.p'), 'wb') as handle:
        pickle.dump(mdic, handle, protocol=pickle.HIGHEST_PROTOCOL)