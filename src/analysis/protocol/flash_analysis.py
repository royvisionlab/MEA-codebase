import os
import numpy as np
import config as cfg
from progressbar import *
import argparse
from symphony_data import Dataset, Analysis
import hdf5storage

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Convert Kilosort output to Vision .neurons format')
    parser.add_argument('experimentName', type=str, help='folder containing Kilosort outputs, must already exist')
    parser.add_argument('-a','--algorithm', default='yass', type=str, help='Sorting algorithm used (yass or kilosort2)')
    parser.add_argument('-p','--protocol', default='ChromaticSpot', type=str, help='protocol identifier')
    parser.add_argument('-g','--group', default=None, type=str, help='EpochGroup identifier')
    parser.add_argument('-f','--file', nargs='+', default=None, type=str, help='Data file identifier')
    parser.add_argument('-b','--bin_rate', default=100.0, type=float, help='Bin rate for spikes')
    parser.add_argument('-r','--sample_rate', default=20000.0, type=float, help='Data sample rate')

    args = parser.parse_args()

    # Get the output file path.
    SORT_PATH, _ = cfg.get_data_paths()
    filepath = SORT_PATH + args.experimentName + '/' + args.experimentName + '_' + args.algorithm +  '_Flash.mat'

    d = Dataset(args.experimentName)

    param_names = ['chromaticClass','contrast']

    spike_array, cluster_id, params, unique_params, pre_pts, stim_pts, tail_pts = d.get_spike_times_and_parameters(
        args.protocol, args.group, param_names, sort_algorithm=args.algorithm, file_name=args.file, bin_rate=args.bin_rate, sample_rate=args.sample_rate)

    # Unpack the parameters.
    contrast = np.array(params['contrast'])
    chromaticClass = np.array(params['chromaticClass'])
    
    mdic = {'spike_array': spike_array, 'cluster_id': cluster_id, 'pre_pts': pre_pts, 'stim_pts': stim_pts, 'tail_pts': tail_pts, 
        'contrast': contrast, 'chromaticClass': chromaticClass}

    # analysis.save_mat(filepath, mdic)
    
    if os.path.exists(filepath+'.mat'):
        os.remove(filepath+'.mat')
    hdf5storage.savemat( file_name=filepath+'.mat', mdict=mdic, format=7.3, matlab_compatible=True, compress=False )


# from movingbar_analysis import MovingBarAnalysis
# experimentName = '20220712C'
# a = MovingBarAnalysis(experimentName, sort_algorithm='yass')
# avg_spikes, cluster_id, u_contrasts, u_orientations = a.get_spikes()
