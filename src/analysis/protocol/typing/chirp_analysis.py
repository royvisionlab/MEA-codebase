import numpy as np
import platform
import argparse
import os
from symphony_data import Dataset, Analysis

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Chirp analysis')
    parser.add_argument('experimentName', type=str, help='Name of experiment')
    parser.add_argument('-f', '--filenames', nargs='+', type=str, help='List of input datafiles')
    parser.add_argument('-a','--algorithm', default='kilosort2', type=str, help='Sorting algorithm used (yass or kilosort2)')
    parser.add_argument('-s','--search', default='Chirp', type=str, help='Name of stimulus protocol to analyze')
    parser.add_argument('-g','--group', default=None, type=str, help='Search string for EpochGroup (optional)')
    parser.add_argument('-b','--bin_rate', default=100.0, type=float, help='Bin rate for spikes')
    parser.add_argument('-r','--sample_rate', default=20000.0, type=float, help='Data sample rate')

    args = parser.parse_args()

    # Set the filepath for the output MAT file.

    str_outputfile = args.experimentName + '_' + args.algorithm + '_'
    if isinstance(args.filenames, list):
        for str_file in args.filenames:
            str_outputfile+= str_file + '_'
        str_outputfile+= 'chirp.mat'
    else:
        str_outputfile+= args.filenames +  '_chirp.mat'

    if (platform.node() == 'maverick'):
        filepath = os.path.join('/home/mike/ftp/files/', str_outputfile)
    elif platform.node() == 'Riekes-MacBook-Pro.local': # For Vyom's laptop
        filepath = os.path.join('/Users/riekelabbackup/Desktop/Output/', str_outputfile)
    # If on hyak, save to experiment dir
    elif 'mmfs1' in os.getcwd() or 'gscratch' in os.getcwd():
        filepath = os.path.join('/gscratch/retina/data/sorted/', args.experimentName, str_outputfile)
    # If on other personal computer, save to desktop
    else:
        filepath = os.path.join(os.path.expanduser('~'), 'Desktop', str_outputfile)

    # Get the avg spike response for each stimulus condition.
    d = Dataset(args.experimentName)

    param_names = ['stepTime','frequencyTime','contrastTime','interTime','stepContrast','frequencyContrast',
        'frequencyMin','frequencyMax','contrastMin','contrastMax','contrastFrequency'] # Just need to grab it all.

    spike_dict, cluster_id, params, unique_params, pre_pts, stim_pts, tail_pts = d.get_spike_rate_and_parameters(
        args.search, args.group, param_names, sort_algorithm=args.algorithm, bin_rate=args.bin_rate, 
        sample_rate=args.sample_rate, file_name=args.filenames)

    # Compute PSTH across all epochs for cells
    avg_psth = np.zeros((len(cluster_id), spike_dict[cluster_id[0]].shape[1]), dtype=np.float32)
    for count, cell in enumerate(spike_dict):
        avg_psth[count,:] = np.mean(spike_dict[cell], axis=0)

    # Save out a MAT file.
    analysis = Analysis()

    mdic = {'params': unique_params, 'avg_psth': avg_psth, 'cluster_id': cluster_id, 'pre_pts': np.unique(pre_pts), 
        'stim_pts': np.unique(stim_pts), 'tail_pts': np.unique(tail_pts), 'bin_rate': args.bin_rate}

    analysis.save_mat(filepath, mdic)