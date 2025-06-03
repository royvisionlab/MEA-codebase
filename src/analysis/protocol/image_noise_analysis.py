

import numpy as np
import platform
import argparse
from symphony_data import Dataset, Analysis
import hdf5storage
import torch

def compute_sta(seed):
    screen_size = np.array([800.0,600.0])
    stimulus_size = np.floor(screen_size/6.6/2.0)*2
    noiseMatrix = imgaussfilt(obj.noiseStream.randn(size(obj.imagePatchMatrix)), obj.noiseFilterSD)
    # Seed the random number generator.
    np.random.seed( seed )
    return None

def compute_sta_torch(self, stimulus, binned_spikes, stride, depth=30, norm=False):
    """
    Computes STA by convolving the spiketrain with the stimulus. This
    implementation uses PyTorch. It doesn't use GPUs because of the memory
    intensiveness, but it is extremely fast relative to numpy even on CPU.

    Parameters:
        stimulus: stimulus tensor of size t, y, x, color (float32 in contrast [-1,1]).
        binned_spikes: spiketrain matrix of size t,cells
        stride: the stride of convolution (interval)
        depth: number of lags for convolution
        norm: boolean whether to normalize STA for nice plotting
    Returns:
        STA tensor, indexed by cell,t,y,x,color channel
    """

    # Intialize STA tensor.
    num_cells = binned_spikes.shape[1]
    n_frames,height,width,n_phosphors = stimulus.shape
    sta_tensor = np.zeros((num_cells,depth,height,width,n_phosphors))

    # Vectorize stimulus.
    stimulus = np.reshape(stimulus,(n_frames,
                                    height * width * n_phosphors)).astype(np.float32).T
    
    # Initialize STA, the offset range, and compute.
    offset_cnt = 0
    sta_tensor = np.zeros((num_cells,depth,height,width,n_phosphors))
    offset_range = np.arange(0,(depth * stride),stride)

    # Initialize progress bar.
    widgets = ['Computing STAs: ', Percentage(), ' ', Bar(marker='*',
            left='[',right=']'),' ', ETA()]
    pbar = ProgressBar(widgets=widgets, maxval=depth)
    pbar.start()

    # Set up torch stuff.
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu') #torch.device('cpu')
    
    binned_spikes = torch.tensor(binned_spikes).to(device).to(torch.float32)
    stimulus = torch.tensor(stimulus).to(device).to(torch.float32)

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
        sta_tensor[:,offset_cnt,...] = np.reshape(sta.T,
                                                (num_cells,height,
                                                width,n_phosphors))
        pbar.update(offset_cnt)
        offset_cnt +=1

    pbar.finish()

    if device.type == 'cuda':
        torch.cuda.empty_cache()
    
    # Clear some memory
    del sta
    gc.collect()

    if not norm:
        return sta_tensor

    return self.normalize_sta(sta_tensor)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Natural image flash analysis')
    parser.add_argument('experimentName', type=str, help='Name of experiment (eg. 20220531C)')
    parser.add_argument('-a','--algorithm', default='yass', type=str, help='Sorting algorithm used (yass or kilosort2)')
    parser.add_argument('-p','--protocol', default='NaturalImageFlashPlusNoise', type=str, help='Name of stimulus protocol to analyze')
    parser.add_argument('-g','--group', default=None, type=str, help='Search string for EpochGroup (optional)')
    parser.add_argument('-b','--bin_rate', default=1000.0, type=float, help='Bin rate for spikes')
    parser.add_argument('-r','--sample_rate', default=20000.0, type=float, help='Data sample rate')

    args = parser.parse_args()

    # filepath = '/Users/michaelmanookin/Documents/Data/rawdata/MEA_Data/Analysis/' + args.experimentName + '_sta.mat'
    if (platform.node() == 'maverick'):
        filepath = '/home/mike/ftp/files/' + args.experimentName + '_' + args.algorithm + '_ImageNoise.mat'
    else:
        filepath = '/gscratch/retina/data/sorted/' + args.experimentName + '/' + args.experimentName + '_' + args.algorithm + '_ImageNoise.mat'

    # Get the avg spike response for each stimulus condition.
    d = Dataset(args.experimentName)

    param_names = ['noiseSeed', 'imagePatchIndex', 'currentPatchLocation','flashTime','apertureDiameter',
            'noiseFilterSD','noiseContrast','numNoiseRepeats','linearizeCones','grateBarSize','WeberConstant',
            'maxIntensity','currentNoiseContrast']

    spike_array, cluster_id, params, unique_params, pre_pts, stim_pts, tail_pts = d.get_spike_times_and_parameters(
        args.protocol, args.group, param_names, sort_algorithm=args.algorithm, bin_rate=args.bin_rate, sample_rate=args.sample_rate)

    # spike_dict, cluster_id, params, unique_params, pre_pts, stim_pts, tail_pts = d.get_spike_count_and_parameters(
    #     args.protocol, args.group, param_names, sort_algorithm=args.algorithm, bin_rate=args.bin_rate, sample_rate=args.sample_rate)

    # flashTime = 200 % ms 
    # apertureDiameter = 200 % um
    # noiseFilterSD = 2 % pixels
    # noiseContrast = 1;
    # numNoiseRepeats = 5;

    # spike_array = np.zeros((len(spike_dict), spike_dict[cluster_id[0]].shape[0]), dtype=object)
    # # Loop through and fill the matrix.
    # for count, key in enumerate(spike_dict):
    #     for i in range(spike_array.shape[1]):
    #         spike_array[count,i] = spike_dict[key][i]
    
    noiseSeed = np.array(params['noiseSeed'])
    imagePatchIndex = np.array(params['imagePatchIndex'])
    currentPatchLocation = np.array(params['currentPatchLocation'])
    flashTime = np.array(params['flashTime'])
    apertureDiameter = np.array(params['apertureDiameter'])
    noiseFilterSD = np.array(params['noiseFilterSD'])
    noiseContrast = np.array(params['noiseContrast'])
    numNoiseRepeats = np.array(params['numNoiseRepeats'])

    # spike_array = np.zeros((len(spike_dict), spike_dict[cluster_id[0]].shape[0], spike_dict[cluster_id[0]].shape[1]))
    # for count, key in enumerate(spike_dict):
    #     spike_array[count,...] = spike_dict[key]

    # spike_array = spike_array.astype(int)

    # Save out a MAT file.
    analysis = Analysis()

    mdic = {'spike_array': spike_array, 'cluster_id': cluster_id, 'pre_pts': pre_pts, 'stim_pts': stim_pts, 'tail_pts': tail_pts, 
        'noiseSeed': noiseSeed, 'imagePatchIndex': imagePatchIndex, 'currentPatchLocation': currentPatchLocation, 'flashTime': flashTime,
        'apertureDiameter': apertureDiameter, 'noiseFilterSD': noiseFilterSD, 'noiseContrast': noiseContrast, 'numNoiseRepeats': numNoiseRepeats}

    # analysis.save_mat(filepath, mdic)
    hdf5storage.savemat( filepath, mdic, format=7.3, matlab_compatible=True, compress=False )



# ids = list()
# for groupCount, group in enumerate(protocol['group']):
#     nBlocks = len(group['block'])
#     # Grab the epoch blocks
#     for blockCount, block in enumerate(group['block']):

#         dataFile = self.M.get_data_file_from_block(block)

#         # Try to load a data file.
#         filepath = os.path.join(self.M.sortpath, self.experimentName, dataFile, sort_algorithm)

#         self.load_vision_table(filepath, dataFile)
#         ids.append(self.cell_ids)


# count = 0
# empty_array = np.zeros((1076, len(epoch_frames)), dtype=object)
# for count, epoch in enumerate(block['epoch']):

#     parameters = epoch['parameters']

#     duration=(parameters['preTime']+parameters['stimTime']+parameters['tailTime'])
#     for cell_count, cell in enumerate(self.cell_ids):
#         spike_times = self.vcd.get_spike_times_for_cell(cell) / sample_rate * 1000 # ms
#         spike_times = spike_times - epoch_frames[count][0]
#         spike_times = spike_times[np.where((spike_times >= 0) * (spike_times < duration))[0]]
#         empty_array[cell_count,count] = spike_times

