"""
CRF protocol data extraction and analysis.

Primary output is '4Hz_amp' which is array of shape [n_clusters, n_contrast_levels]
Values are the peak-trough amplitude of the 4Hz component of the avg PSTH for 1 cycle of CRF.
This is equivalent to 2*amp of a 4Hz sine wave fit to the PSTH.

Clusters are indexed by 'cluster_id'

'mean_psth_cycle' is array of shape [n_clusters, n_bins_per_cycle, n_contrast_levels]

Optionally, entire PSTH timecourse can be saved with --savefull True
'avg_psth' is array of shape [n_clusters, n_bins, n_contrast_levels]
This will increase output file size a lot.

Usage:
    python crf_analysis.py 20221216C -f data009 -a kilosort2
"""

import numpy as np
from scipy.fftpack import fft
import argparse
from symphony_data import Dataset, Analysis
import pickle
import config as cfg
import os
import pickle

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='CRF analysis')
    parser.add_argument('experimentName', type=str, help='Name of experiment (eg. 20220531C)')
    parser.add_argument('-f', '--filenames', nargs='+', type=str, help='List of input datafiles')
    parser.add_argument('--savefull', default=False, type=bool, help='Save full PSTH timecourse (default: False)')
    parser.add_argument('-a','--algorithm', default='yass', type=str, help='Sorting algorithm used (yass or kilosort2)')
    parser.add_argument('-s','--search', default='ContrastResponseGrating', type=str, help='Name of stimulus protocol to analyze')
    parser.add_argument('-g','--group', default=None, type=str, help='Search string for EpochGroup (optional)')
    parser.add_argument('-b','--bin_rate', default=1000.0, type=float, help='Bin rate for spikes')
    parser.add_argument('-r','--sample_rate', default=20000.0, type=float, help='Data sample rate')

    args = parser.parse_args()

    # Get the data paths
    SORT_PATH, JSON_PATH, OUTPUT_PATH = cfg.get_data_paths()
    
    # Set filepath for output file
    str_outputfile = args.experimentName + '_' + args.algorithm + '_'
    if isinstance(args.filenames, list):
        for str_file in args.filenames:
            str_outputfile+= str_file + '_'
        str_outputfile+= 'crf.p'
    else:
        str_outputfile+= args.filenames +  '_crf.p'
    
    filepath = os.path.join(OUTPUT_PATH, str_outputfile)
    print('Saving output to: ' + filepath)

    print('Loading data...')
    # Get the avg spike response for each stimulus condition.
    d = Dataset(args.experimentName)

    param_names = ['contrast', 'barWidth', 'temporalFrequency', 'orientation', 
                   'apertureClass', 'spatialClass', 'temporalClass', 'chromaticClass']
    
    spike_dict, cluster_id, params, unique_params, pre_pts, stim_pts, tail_pts = d.get_spike_rate_and_parameters(
        args.search, args.group, param_names, sort_algorithm=args.algorithm, bin_rate=args.bin_rate, 
        sample_rate=args.sample_rate, file_name=args.filenames)
    
    if len(unique_params['temporalFrequency']) > 1:
        raise ValueError('More than one temporal frequency in data.')
    else:
        print('Temporal frequency: ' + str(unique_params['temporalFrequency'][0]) + ' Hz')

    # Compute avg PSTH of single cycle
    n_bins_per_cycle = int(args.bin_rate / params['temporalFrequency'][0])

    n_unique_contrasts = len(np.unique(params['contrast']))
    arr_mean_psth_cycle = np.zeros((len(cluster_id), n_bins_per_cycle, n_unique_contrasts), dtype=np.float32)
    avg_psth = np.zeros((len(cluster_id), spike_dict[cluster_id[0]].shape[1], n_unique_contrasts), dtype=np.float32)
    
    print('Computing average PSTH of single cycle for {} contrasts'.format(n_unique_contrasts))
    for i, c in enumerate(cluster_id):
        for j, contrast in enumerate(np.unique(params['contrast'])):
            arr_psth = np.mean(spike_dict[c][params['contrast'] == contrast, :], axis=0)
            avg_psth[i, :, j] = arr_psth
            arr_mean_psth_cycle[i, :, j] = arr_psth.reshape((-1, n_bins_per_cycle)).mean(axis=0)

    # Compute 4Hz amp
    def compute_4Hz_amp(arr_psth_cycle: np.ndarray):
        """Computes peak-trough amplitude of 4Hz component of PSTH.

        Args:
            arr_psth_cycle (np.ndarray): 1D array of 1 cell's avg PSTH for 1 cycle of CRF.

        Returns:
            _type_: _description_
        """
        n_bins_per_cycle = arr_psth_cycle.shape[0]
        arr_ft = np.abs(fft(arr_psth_cycle)[:n_bins_per_cycle//2]) # Magnitude of positive frequencies
        arr_ft = 2/n_bins_per_cycle*arr_ft # Normalize by number of bins, x2 for pos and neg frequencies
        return arr_ft[1] * 2 # x2 for peak-trough amp
    
    arr_4Hz_amp = np.zeros((len(cluster_id), 3))
    for cell_idx in range(len(cluster_id)):
        for contrast_idx in range(n_unique_contrasts):
            arr_4Hz_amp[cell_idx, contrast_idx] = compute_4Hz_amp(arr_mean_psth_cycle[cell_idx, :, contrast_idx])


    print('Saving data...')
    # Save out a MAT file.
    analysis = Analysis()

    mdic = {'params': params, 'unique_params': unique_params, 'mean_psth_cycle': arr_mean_psth_cycle,
            '4Hz_amp': arr_4Hz_amp, 'cluster_id': cluster_id, 'pre_pts': np.unique(pre_pts), 
            'stim_pts': np.unique(stim_pts), 'tail_pts': np.unique(tail_pts), 'bin_rate': args.bin_rate}
    if args.savefull:
        mdic['avg_psth'] = avg_psth

    with open(filepath, 'wb') as out_file:
        pickle.dump(mdic, out_file)
    print('Saved to ' + filepath)
