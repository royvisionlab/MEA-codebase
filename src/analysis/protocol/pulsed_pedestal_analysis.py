import numpy as np
import platform
import argparse
from symphony_data import Dataset, Analysis
import pickle
import os

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Pulsed Pedestal analysis')
    parser.add_argument('experimentName', type=str, help='Name of experiment')
    parser.add_argument('-f', '--filenames', nargs='+', type=str, help='List of input datafiles')
    parser.add_argument('-a','--algorithm', default='yass', type=str, help='Sorting algorithm used (yass or kilosort2)')
    parser.add_argument('-s','--search', default='PulsedPedestal', type=str, help='Name of stimulus protocol to analyze')
    parser.add_argument('-g','--group', default=None, type=str, help='Search string for EpochGroup (optional)')
    parser.add_argument('-b','--bin_rate', default=100.0, type=float, help='Bin rate for spikes')
    parser.add_argument('-r','--sample_rate', default=20000.0, type=float, help='Data sample rate')

    args = parser.parse_args()

    str_outputfile = args.experimentName + '_' + args.algorithm + '_'
    if isinstance(args.filenames, list):
        for str_file in args.filenames:
            str_outputfile+= str_file + '_'
        str_outputfile+= '_pulsed_pedestal.p'
    else:
        str_outputfile+= args.filenames +  '_pulsed_pedestal.p'

    if (platform.node() == 'maverick'):
        filepath = os.path.join('/home/mike/ftp/files/', str_outputfile)
    elif platform.node() == 'Riekes-MacBook-Pro.local':
        filepath = os.path.join('/Users/riekelabbackup/Desktop/Output/', str_outputfile)
    else:
        filepath = os.path.join('/gscratch/retina/vyomr/data/sorted/', str_outputfile)
    
    # Create directory in case it doesn't exist
    os.makedirs(os.path.dirname(filepath), exist_ok=True)

    # Get the avg spike response for each stimulus condition.
    d = Dataset(args.experimentName)

    param_names = ['contrasts','flashColor','flashTime','gridWidth','interpulseInterval','numberOfAverages',
        'postFlashTime','preFlashTime','preTime','sampleRate','separationSize',
         'stimContrast', 'stimTime', 'stixelSize', 'tailTime', 'testSquareIdx', 'testStimContrast'] # Just need to grab it all.

    spike_dict, cluster_id, params, unique_params, pre_pts, stim_pts, tail_pts = d.get_spike_rate_and_parameters(
        args.search, args.group, param_names, sort_algorithm=args.algorithm, bin_rate=args.bin_rate, sample_rate=args.sample_rate,
        file_name=args.filenames)

    # Save out a pickle file.
    analysis = Analysis()

    mdic = {'params': params, 'unique_params': unique_params, 'spike_dict': spike_dict, 'cluster_id': cluster_id,
     'pre_pts': pre_pts, 'stim_pts': stim_pts, 'tail_pts': tail_pts, 'bin_rate': args.bin_rate}

    with open(filepath, 'wb') as handle:
        pickle.dump(mdic, handle, protocol=pickle.HIGHEST_PROTOCOL)