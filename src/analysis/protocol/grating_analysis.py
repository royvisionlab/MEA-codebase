import numpy as np
import platform
import argparse
import os
from symphony_data import Dataset, Analysis

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='GratingDSOS analysis')
    #parser = argparse.ArgumentParser(description='Chirp analysis')
    parser.add_argument('experimentName', type=str, help='Name of experiment')
    parser.add_argument('-f', '--filenames', nargs='+', type=str, help='List of input datafiles')
    parser.add_argument('-a','--algorithm', default='yass', type=str, help='Sorting algorithm used (yass or kilosort2)')
    parser.add_argument('-s','--search', default='GratingDSOS', type=str, help='Name of stimulus protocol to analyze')
    #parser.add_argument('-s','--search', default='Chirp', type=str, help='Name of stimulus protocol to analyze')
    parser.add_argument('-g','--group', default=None, type=str, help='Search string for EpochGroup (optional)')
    parser.add_argument('-b','--bin_rate', default=100.0, type=float, help='Bin rate for spikes')
    parser.add_argument('-r','--sample_rate', default=20000.0, type=float, help='Data sample rate')

    args = parser.parse_args()

    # Set the filepath for the output MAT file.

    str_outputfile = args.experimentName + '_' + args.algorithm + '_'
    if isinstance(args.filenames, list):
        for str_file in args.filenames:
            str_outputfile+= str_file + '_'
        str_outputfile+= 'gratingDSOS.mat'
        #str_outputfile+= 'chirp.mat'
    else:
        str_outputfile+= args.filenames +  '_gratingDSOS.mat'
        #str_outputfile+= args.filenames +  '_chirp.mat'
        

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

    param_names = ['orientation','contrast','temporalFrequency','barWidth']

    #param_names = ['stepTime','frequencyTime','contrastTime','interTime','stepContrast','frequencyContrast',
    #    'frequencyMin','frequencyMax','contrastMin','contrastMax','contrastFrequency'] # Just need to grab it all.
    
    spike_dict, cluster_id, params, unique_params, pre_pts, stim_pts, tail_pts = d.get_spike_rate_and_parameters(
        args.search, args.group, param_names, sort_algorithm=args.algorithm, bin_rate=args.bin_rate, 
        sample_rate=args.sample_rate, file_name=args.filenames)

    # Compute PSTH across all epochs for cells
    #avg_psth = np.zeros((len(cluster_id), spike_dict[cluster_id[0]].shape[1]), dtype=np.float32)
    #for count, cell in enumerate(spike_dict):
    #    avg_psth[count,:] = np.mean(spike_dict[cell], axis=0)

    # Unique orientations.
    u_orientation = unique_params['orientation']
    # Unique contrasts
    u_contrast = unique_params['contrast']
    # Unique temporalFrequency
    u_temporalFrequency = unique_params['temporalFrequency']
    u_barWidth = unique_params['barWidth']

    print(u_barWidth)
    print(len(u_barWidth))
    
    avg_psth = np.zeros((len(cluster_id),len(u_contrast),len(u_barWidth),len(u_temporalFrequency),len(u_orientation),spike_dict[cluster_id[0]].shape[1]), dtype=np.float32)
    avg_std = np.zeros((len(cluster_id),len(u_contrast),len(u_barWidth),len(u_temporalFrequency),len(u_orientation),spike_dict[cluster_id[0]].shape[1]), dtype=np.float32)
    
    for count, cell in enumerate(spike_dict):
        spikes = spike_dict[cell]
        for ct_count, contrast in enumerate(u_contrast):
            for bw_count, barWidth in enumerate(u_barWidth):
                for tf_count, temporalFrequency in enumerate(u_temporalFrequency):
                    for o_count, orientation in enumerate(u_orientation):               
                        idx = np.where((np.array(params['contrast']) == contrast) & (np.array(params['barWidth']) == barWidth) & (np.array(params['temporalFrequency']) == temporalFrequency) & (np.array(params['orientation']) == orientation))[0]
                        avg_psth[count, ct_count, bw_count, tf_count, o_count, :] = np.mean(spikes[idx,:], axis=0)
                        avg_std[count, ct_count, bw_count, tf_count, o_count, :] = np.std(spikes[idx,:], axis=0)

    pre_pts = np.unique(pre_pts).astype(int)
    stim_pts = np.unique(stim_pts).astype(int)
    tail_pts = np.unique(tail_pts).astype(int)

    avg_response = np.mean(avg_psth[:,:,:,:,:,pre_pts[0]:pre_pts[0]+stim_pts[0]], axis=5)
    avg_std = np.mean(avg_std[:,:,:,:,:,pre_pts[0]:pre_pts[0]+stim_pts[0]], axis=5)
    avg_spontaneous = np.mean(avg_psth[:,:,:,:,:,:pre_pts[0]], axis=(4,5))

    # Save out a MAT file.
    analysis = Analysis()

    # Loop through and compute OSI and DSI for each cell/contrast.
    osi = np.zeros((avg_response.shape[0],avg_response.shape[1],avg_response.shape[2],avg_response.shape[3]))
    osi_angle = np.zeros((avg_response.shape[0],avg_response.shape[1],avg_response.shape[2],avg_response.shape[3]))
    dsi = np.zeros((avg_response.shape[0],avg_response.shape[1],avg_response.shape[2],avg_response.shape[3]))
    dsi_angle = np.zeros((avg_response.shape[0],avg_response.shape[1],avg_response.shape[2],avg_response.shape[3]))
    for cell in range(avg_response.shape[0]):
        for ct in range(avg_response.shape[1]):
            for bw in range(avg_response.shape[2]):
                for tf in range(avg_response.shape[3]):
                    osi_global, osi_ang = analysis.compute_OSI(u_orientation, avg_response[cell,ct,bw,tf,:])
                    dsi_global, dsi_ang = analysis.compute_DSI(u_orientation, avg_response[cell,ct,bw,tf,:])
                    osi[cell,ct,bw,tf] = osi_global
                    osi_angle[cell,ct,bw,tf] = osi_ang
                    dsi[cell,ct,bw,tf] = dsi_global
                    dsi_angle[cell,ct,bw,tf] = dsi_ang

    mdic = {'avg_psth': avg_psth, 'avg_response':avg_response,'avg_spontaneous':avg_spontaneous,'orientations': np.array(params['orientation']),
            'avg_std':avg_std, 'cluster_id': cluster_id, 'params': unique_params, 'bin_rate': args.bin_rate, 'barWidth': params['barWidth'],
            'temporalFrequency': np.array(params['temporalFrequency']), 'u_contrast': u_contrast,'u_orientation': u_orientation, 'osi': osi,
            'osi_angle': osi_angle, 'dsi': dsi, 'dsi_angle': dsi_angle, 'pre_pts': np.unique(pre_pts), 'stim_pts': np.unique(stim_pts),
            'spikes':spikes, 'tail_pts': np.unique(tail_pts)}

    #mdic = {'params': unique_params, 'avg_psth': avg_psth, 'cluster_id': cluster_id, 'pre_pts': np.unique(pre_pts), 
    #    'stim_pts': np.unique(stim_pts), 'tail_pts': np.unique(tail_pts), 'bin_rate': args.bin_rate}

    analysis.save_mat(filepath, mdic)