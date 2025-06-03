
import numpy as np
import platform
import argparse
from collections import Counter
from symphony_data import Dataset, Analysis

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Object motion dots analysis')
    parser.add_argument('experimentName', type=str, help='Name of experiment')
    parser.add_argument('-a','--algorithm', default='yass', type=str, help='Sorting algorithm used (yass or kilosort2)')
    parser.add_argument('-s','--search', default='ChromaticSpot', type=str, help='Name of stimulus protocol to analyze')
    parser.add_argument('-g','--group', default=None, type=str, help='Search string for EpochGroup (optional)')
    parser.add_argument('-b','--bin_rate', default=100.0, type=float, help='Bin rate for spikes')
    parser.add_argument('-r','--sample_rate', default=20000.0, type=float, help='Data sample rate')
    parser.add_argument('-f','--file', default=None, type=str, help='name of single data file to analyze i.e. data003 (default: None)')

    args = parser.parse_args()

    if (platform.node() == 'mike'):
        filepath = '/home/mike/ftp/files/' + args.experimentName + '_' + args.algorithm +  '_Flashes.mat'
    else:
        filepath = '/gscratch/retina/data/sorted/' + args.experimentName + '/' + args.experimentName + '_' + args.algorithm + '_Flashes.mat'

    # Get the avg spike response for each stimulus condition.
    d = Dataset(args.experimentName)

    param_names = ['chromaticClass','contrast']

    spike_dict, cluster_id, params, unique_params, pre_pts, stim_pts, tail_pts = d.get_spike_rate_and_parameters(
        args.search, args.group, param_names, sort_algorithm=args.algorithm, file_name=args.file, bin_rate=args.bin_rate, sample_rate=args.sample_rate)

    
    classes = list()
    for i in range(len(params['chromaticClass'])):
        if params['chromaticClass'][i] == 'achromatic':
            if params['contrast'][i] < 0:
                classes.append('off')
            else:
                classes.append('on')
        else:
            classes.append(params['chromaticClass'][i])
    
    u_classes = list(Counter(classes).keys())

    avg_psth = np.zeros((len(cluster_id),len(u_classes), spike_dict[cluster_id[0]].shape[1]), dtype=np.float32)
    avg_response = np.zeros((len(cluster_id),len(u_classes)), dtype=np.float32)
    for count, cell in enumerate(spike_dict):
        spikes = spike_dict[cell]
        for const_count, constant in enumerate(u_classes):
            idx = np.where(np.array(classes) == constant)[0]
            avg_psth[count,const_count,:] = np.mean(spikes[idx,:], axis=0)
            # Get the spike rate during the stimulus.
            spike_r = np.zeros((len(idx),))
            for i, value in enumerate(idx):
                pts = range(pre_pts[value].astype(int),pre_pts[value].astype(int)+stim_pts[value].astype(int))
                spike_r[i] = np.sum(spikes[value, pts]) / (float(len(pts))/args.bin_rate)
            avg_response[count,const_count] = np.mean(spike_r)

    # Save out a MAT file.
    analysis = Analysis()

    mdic = {'avg_psth': avg_psth, 'cluster_id': cluster_id, 'u_classes': u_classes, 'avg_response': avg_response, 
        'pre_pts': np.unique(pre_pts), 'stim_pts': np.unique(stim_pts), 'tail_pts': np.unique(tail_pts), 'bin_rate': args.bin_rate}

    analysis.save_mat(filepath, mdic)

    
        


