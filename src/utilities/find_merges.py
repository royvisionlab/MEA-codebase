import visionloader as vl
import os
import scipy.io as sio
from scipy.io import savemat
import numpy as np
from tqdm import tqdm
import torch

# path1 = '/gscratch/retina/data/sorted/20221006C/noise/kilosort2/'
# path2 = '/gscratch/retina/data/sorted/20221006C/natimg/kilosort2/'

# eir = vl.EIReader(path1, 'noise')
# eis = eir.get_all_eis_by_cell_id()
# eir.close()

# eir = vl.EIReader(path2, 'natimg')
# eis2 = eir.get_all_eis_by_cell_id()
# eir.close()



def compute_corr_coeff(ei1, ei2):
    return np.sum(ei1 * ei2) / np.sqrt(np.sum(ei1*ei1) * np.sum(ei2*ei2))

# Remove the n-strongest electrodes.
def compute_ei_subset(raw_ei1, raw_ei2, num=2):
    amps1 = np.sum(raw_ei1, axis=1)
    amps2 = np.sum(raw_ei2, axis=1)
    # Indices for two strongest electrodes.
    ind = np.argpartition(amps1,-num)[-num:]
    np.append(ind, np.argpartition(amps2,-num)[-num:])
    ei1 = raw_ei1.copy()
    ei1 = np.delete(ei1,ind,axis=0)
    ei2 = raw_ei2.copy()
    ei2 = np.delete(ei2,ind,axis=0)
    return np.sum(ei1 * ei2) / np.sqrt(np.sum(ei1*ei1) * np.sum(ei2*ei2)) # Return the r^2

# Compute the correlation between two ei's.

def compute_ei_power(ei, channel_noise, n_sigmas=None):
    if n_sigmas is not None:
        sig_inds = np.argwhere(np.abs(np.amin(ei,axis=1))
                            > n_sigmas * channel_noise).flatten()
        ei_power = np.zeros(ei.shape[0])
        ei_power[sig_inds] = np.sum(ei[sig_inds,:]**2,
                                            axis=1)
    else:
        ei_power = np.sum(ei**2,axis=1)
    return ei_power

def threshold_ei_channel_noise(ei, channel_noise, n_sigmas=2.0):
    if n_sigmas is not None:
        non_sig_inds = np.argwhere(np.std(ei,axis=1)
                            < n_sigmas * channel_noise).flatten()
        ei[non_sig_inds,:] = 0.0
    return ei

def threshold_ei_std(ei, n_sigmas=1.0):
    # Compute the standard dev for the EI across electrodes.
    this_std = np.nanstd(ei)
    if n_sigmas is not None:
        non_sig_inds = np.argwhere(np.std(ei,axis=1)
                            < n_sigmas*this_std).flatten()
        ei[non_sig_inds,:] = 0.0
    return ei

def threshold_ei_all_electrodes(ei, n_sigmas=1.0):
    this_std = np.nanstd(ei)
    if n_sigmas is not None:
        ei[np.abs(ei) < this_std*n_sigmas] = 0.0
    return ei

analysis_path = '/Volumes/Samsung870/analysis/20221123C/'
dataset_name = 'chunk1'
# ref_vcd = vl.load_vision_data(analysis_path + 'chunk2', 'chunk2', include_ei=True, include_noise=True)
# test_vcd = vl.load_vision_data(analysis_path + 'chunk1', 'chunk1', include_ei=True, include_noise=True)

ref_vcd = vl.load_vision_data(analysis_path + 'chunk2', 'chunk2', include_ei=True, include_noise=True, include_sta=True)
test_vcd = vl.load_vision_data(analysis_path + 'chunk1', 'chunk1', include_ei=True, include_noise=True, include_sta=True)

# Try to pull the triggers.
vcd = vl.load_vision_data('/Volumes/Data1TB/20221123C/data002/yass','data002', include_noise=True)

ref_noise = ref_vcd.channel_noise
test_noise = test_vcd.channel_noise

# Get the cell IDs.
ref_ids = ref_vcd.get_cell_ids()
test_ids = test_vcd.get_cell_ids()

cell_idx = dict()
for ref_cell in ref_ids:
    # ref_power = compute_ei_power(ref_vcd.get_ei_for_cell(ref_cell).ei, ref_noise)
    ref_ei = ref_vcd.get_ei_for_cell(ref_cell).ei
    corr = []
    for test_cell in test_ids:
        test_ei = test_vcd.get_ei_for_cell(test_cell).ei
        # test_power = compute_ei_power(test_vcd.get_ei_for_cell(test_cell).ei, test_noise)
        # corr.append(np.corrcoef(ref_power, test_power)[0,1])
        corr.append(np.corrcoef(ref_ei.flatten(), test_ei.flatten())[0,1])
    corr = np.array(corr)
    corr[np.isnan(corr)] = 0.0
    corr[corr == np.inf] = 0.0
    max_ind = np.argmax(corr)
    cell_idx[ref_cell] = test_ids[max_ind]

# Compare STAs.
sta_idx = dict()
for ref_cell in ref_ids:
    ref_sta = ref_vcd.get_sta_for_cell(ref_cell)
    corr = []
    for test_cell in test_ids:
        test_sta = test_vcd.get_sta_for_cell(test_cell)
        corr.append(
            (np.corrcoef(ref_sta.red[:,:,25:28].flatten(), test_sta.red[:,:,25:28].flatten())[0,1] + 
            np.corrcoef(ref_sta.blue[:,:,25:28].flatten(), test_sta.blue[:,:,25:28].flatten())[0,1] )/2.0 )
        # corr.append( np.maximum(
        #     np.corrcoef(ref_sta.red[:,:,25:28].flatten(), test_sta.red[:,:,25:28].flatten())[0,1], 
        #     np.corrcoef(ref_sta.blue[:,:,25:28].flatten(), test_sta.blue[:,:,25:28].flatten())[0,1] ) )
    corr = np.array(corr)
    corr[np.isnan(corr)] = 0.0
    corr[corr == np.inf] = 0.0
    max_ind = np.argmax(corr)
    sta_idx[ref_cell] = test_ids[max_ind]

ei_cells = np.zeros((len(ref_ids,)))
sta_cells = np.zeros((len(ref_ids,)))
for idx, ref_cell in enumerate(ref_ids):
    ei_cells[idx] = cell_idx[ref_cell]
    sta_cells[idx] = sta_idx[ref_cell]

ref_acf = ref_vcd.get_acf_numpairs_for_cell(ref_cell)

rootdir = '/gscratch/retina/data/sorted/20220712C/'
sortDir = 'yass'
csv_path = '/gscratch/retina/data/sorted/20220712C/all.csv'

# Read the noise file.
eir = vl.EIReader('/gscratch/retina/data/sorted/20220712C/data021/yass/', 'data021')
noise_ei = eir.get_all_eis_by_cell_id()
eir.close()
noise_keys = np.array(list(noise_ei.keys())).astype(int)
n_ei = np.zeros((np.max(noise_keys),noise_ei[noise_keys[0]].ei.shape[0],noise_ei[noise_keys[0]].ei.shape[1]))
for noise_cell_id, ns_ei in noise_ei.items():
    n_ei[noise_cell_id-1,...] = ns_ei.ei

path_keys, dset_lengths = [], []
with open(csv_path, 'r') as metadata_csv:
    for line in metadata_csv.readlines():
        data_line = line.strip('\n').split(',')
        path_keys.append(data_line[0])
        dset_lengths.append(int(data_line[1]))

corr = np.zeros((len(path_keys),np.max(noise_keys),3500))

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu') #torch.device('cpu')
        
n_ei = torch.tensor(n_ei).to(device).to(torch.float32)

for count, value in enumerate(path_keys):
    data_file = value.split('/')[-1]
    ei_path = os.path.join(rootdir, data_file, sortDir)
    # ei_path = '/usr/share/pool/SortedData/20220712C/' + data_file + '/' + args.sortDir + '/' 
    eir = vl.EIReader(ei_path, data_file)
    eis = eir.get_all_eis_by_cell_id()
    eir.close()
    ei_keys = np.array(list(eis.keys())).astype(int)
    ei_matrix = np.zeros((3500,512,61))
    for ei_cell_id, readable_ei in eis.items():
        ei_matrix[ei_cell_id-1,...] = readable_ei.ei
    corr_tmp = np.zeros((np.max(noise_keys),3500))
    corr_tmp = torch.tensor(corr_tmp).to(device).to(torch.float32)
    ei_matrix = torch.tensor(ei_matrix).to(device).to(torch.float32)
    for i in tqdm(range(len(noise_keys))):
        for e_key in ei_keys:
            corr_tmp[noise_keys[i]-1,e_key-1] = torch.sum( torch.mul(n_ei[noise_keys[i]-1,...], ei_matrix[e_key-1,...]) ) / torch.sqrt(torch.mul(torch.sum( torch.mul(n_ei[noise_keys[i]-1,...], n_ei[noise_keys[i]-1,...]) ), torch.sum( torch.mul(ei_matrix[e_key-1,...], ei_matrix[e_key-1,...]) )))
            # corr_tmp = torch.sum( torch.mul(n_ei[noise_keys[i]-1,...], ei_matrix[e_key-1,...]) )
    # Set nans to zero
    corr_tmp[corr_tmp != corr_tmp] = 0.0
    corr[count,...] = corr_tmp.detach().cpu().numpy() 

filepath = '/gscratch/retina/data/sorted/20220712C/ei_corr.mat'
mdic = {'corr': corr, 'path_keys': path_keys}
savemat(filepath, mdic)

for i in tqdm(range(noise_keys)):
    for e_count in range(ei_keys):
        e_key = ei_keys[count]
        corr_tmp = np.sum( n_ei[noise_keys[i]-1,...] * ei_matrix[e_key-1,...] )
        corr[count,noise_keys[i]-1,e_key-1] = corr_tmp

for i in tqdm(range(noise_keys)):

for offset in offset_range:

    with torch.no_grad():
        sta = torch.matmul(stimulus[:,:n_frames-offset],
                        binned_spikes[offset:,:])

        # Index appropriately to avoid nans and normalize by number of spikes.
        nonzero_inds = torch.nonzero(torch.sum(binned_spikes[offset:,:],
                                            dim=0),as_tuple=False)
        sta[:,nonzero_inds] /= torch.sum(binned_spikes[offset:,:],
                                        dim=0)[nonzero_inds]

    sta = sta.detach().cpu().numpy() 

# for count, value in enumerate(path_keys):
#     data_file = value.split('/')[-1]
#     ei_path = os.path.join(rootdir, data_file, sortDir)
#     # ei_path = '/usr/share/pool/SortedData/20220712C/' + data_file + '/' + args.sortDir + '/' 
#     eir = vl.EIReader(ei_path, data_file)
#     eis = eir.get_all_eis_by_cell_id()
#     eir.close()
#     for noise_cell_id, n_ei in noise_ei.items():
#         for ei_cell_id, readable_ei in eis.items():
#             ei_matrix = readable_ei.ei
#             corr[count,noise_cell_id-1,ei_cell_id] = compute_corr_coeff(n_ei.ei, readable_ei.ei)





# # Load the noise file
# mdic = sio.loadmat('/Users/michaelmanookin/Documents/Data/rawdata/MEA_Data/Analysis/20221006C/noise_yass_ei.mat')
# eis = mdic['E']
# mdic = sio.loadmat('/Users/michaelmanookin/Documents/Data/rawdata/MEA_Data/Analysis/20221006C/natimg_yass_ei.mat')
# eis2 = mdic['E']

# # raw_corr = np.zeros((len(eis), len(eis2)))
# # sub_corr = np.zeros((len(eis), len(eis2)))

# exclude_pts = [0,1,2,3,4,5]
# corr = np.zeros((len(exclude_pts),eis.shape[0], eis2.shape[0])) 

# for sta_count in range(eis.shape[0]):
#     sta_matrix = eis[sta_count,...]
#     # Load in the EIs from the other stimulus and compare.
#     for ei_count in range(eis2.shape[0]):
#         ei_matrix = eis2[ei_count,...]
#         for ex_count, exclude_pt in enumerate(exclude_pts):
#             if exclude_pt == 0:
#                 # Raw correlation
#                 corr[ex_count,sta_count,ei_count] = compute_corr_coeff(sta_matrix, ei_matrix)
#             else:
#                 # Subset correlation.
#                 corr[ex_count,sta_count,ei_count] = compute_ei_subset(sta_matrix, ei_matrix, exclude_pt)

# # Test self-similarity.
# corr_noise = np.zeros((len(exclude_pts),eis.shape[0], eis.shape[0])) 
# for sta_count in range(eis.shape[0]):
#     sta_matrix = eis[sta_count,...]
#     # Load in the EIs from the other stimulus and compare.
#     for ei_count in range(eis.shape[0]):
#         ei_matrix = eis[ei_count,...]
#         for ex_count, exclude_pt in enumerate(exclude_pts):
#             if exclude_pt == 0:
#                 # Raw correlation
#                 corr_noise[ex_count,sta_count,ei_count] = compute_corr_coeff(sta_matrix, ei_matrix)
#             else:
#                 # Subset correlation.
#                 corr_noise[ex_count,sta_count,ei_count] = compute_ei_subset(sta_matrix, ei_matrix, exclude_pt)


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