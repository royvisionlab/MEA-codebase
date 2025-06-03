
import numpy as np
import platform
import argparse
from symphony_data import Dataset, Analysis
import hdf5storage

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Object motion dots analysis')
    parser.add_argument('experimentName', type=str, help='Name of experiment')
    parser.add_argument('-a','--algorithm', default='yass', type=str, help='Sorting algorithm used (yass or kilosort2)')
    parser.add_argument('-p','--protocol', default='EyeMovement', type=str, help='Name of stimulus protocol to analyze')
    parser.add_argument('-g','--group', default=None, type=str, help='Search string for EpochGroup (optional)')
    parser.add_argument('-b','--bin_rate', default=1000.0, type=float, help='Bin rate for spikes')
    parser.add_argument('-r','--sample_rate', default=20000.0, type=float, help='Data sample rate')

    args = parser.parse_args()

    # filepath = '/Users/michaelmanookin/Documents/Data/rawdata/MEA_Data/Analysis/' + args.experimentName + '_sta.mat'
    if (platform.node() == 'maverick'):
        filepath = '/home/mike/ftp/files/' + args.experimentName + '_' + args.algorithm + '_EyeTraj.mat'
    else:
        filepath = '/gscratch/retina/data/sorted/' + args.experimentName + '/' + args.experimentName + '_' + args.algorithm + '_ImageNoise.mat'

    # Get the avg spike response for each stimulus condition.
    d = Dataset(args.experimentName)

    param_names = ['D', 'backgroundIntensity', 'currentBackgroundScale', 'currentImageName', 'currentP0', 'randomSeed', 'currentStimSet', 'currentImageSet', 'xTraj', 'yTraj', 'apertureDiameter', 'patchMean']

    spike_array, cluster_id, params, unique_params, pre_pts, stim_pts, tail_pts = d.get_spike_times_and_parameters(
        args.protocol, args.group, param_names, sort_algorithm=args.algorithm, bin_rate=args.bin_rate, sample_rate=args.sample_rate)
    
    D = np.array(params['D'])
    backgroundIntensity = np.array(params['backgroundIntensity'])
    currentBackgroundScale = np.array(params['currentBackgroundScale'])
    currentImageName = params['currentImageName']
    currentP0 = np.array(params['currentP0'])
    randomSeed = np.array(params['randomSeed'])
    currentStimSet = params['currentStimSet']
    currentImageSet = params['currentImageSet']
    xTraj = np.array(params['xTraj'])
    yTraj = np.array(params['yTraj'])
    apertureDiameter = np.array(params['apertureDiameter'])
    patchMean = params['patchMean']

    mdic = {'spike_array': spike_array, 'cluster_id': cluster_id, 'pre_pts': pre_pts, 'stim_pts': stim_pts, 'tail_pts': tail_pts, 
        'D': D, 'backgroundIntensity': backgroundIntensity, 'currentBackgroundScale': currentBackgroundScale, 
        'currentImageName': currentImageName, 'currentP0': currentP0, 'randomSeed': randomSeed, 'currentStimSet': currentStimSet,
        'currentImageSet': currentImageSet, 'xTraj': xTraj, 'yTraj': yTraj, 'apertureDiameter': apertureDiameter, 'patchMean': patchMean}

    hdf5storage.savemat( filepath, mdic, format=7.3, matlab_compatible=True, compress=False )