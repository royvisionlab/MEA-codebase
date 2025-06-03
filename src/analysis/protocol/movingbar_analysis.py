
import numpy as np
import config as cfg
import argparse
from symphony_data import Dataset, Analysis

def get_response_from_rf_params(avg_psth: np.ndarray, xy: np.ndarray, orientations: np.ndarray, screen_size: np.ndarray, speed: float, bin_rate: float, pre_pts: float):
    """ Get the response for each cluster based on the receptive field parameters.
    Parameters: 
        avg_psth: Average spike count for each cluster. (np.ndarray)
        xy: Receptive field center for each cluster. (np.ndarray)
        orientations: Orientations for each stimulus condition. (np.ndarray)
        screen_size: Size of the screen in microns. (np.ndarray)
        speed: Speed of the bar in microns per second. (float)
        bin_rate: Bin rate for the spike counts. (float)
        pre_pts: Number of pre-stimulus points. (float)
    Returns:
        response: Response for each cluster. (np.ndarray)
    """
    sample_window = np.arange(-bin_rate*0.1, bin_rate*0.4, 1.0).astype(int)
    # Get the response for each cluster.
    response = np.zeros((avg_psth.shape[0], avg_psth.shape[1], avg_psth.shape[2], avg_psth.shape[3]))
    for ii in range(len(orientations)):
        # Get the start samples.
        t_sample = get_bar_time(xy=xy, screen_size=screen_size, speed=speed, orientation=orientations[ii], bin_rate=bin_rate)
        t_sample = np.floor(t_sample).astype(int) + pre_pts
        for jj in range(avg_psth.shape[0]):
            response[jj,:,ii,:] = np.nanmean(avg_psth[jj,:,ii,:,t_sample[jj]+sample_window],axis=0)
    return response

def get_bar_time(xy: np.ndarray, screen_size: np.ndarray, speed: float, orientation: float, bin_rate: float) -> np.ndarray:
    orientation_rads = orientation * np.pi / 180.0
    a_inv = np.array([[np.cos(orientation_rads),-np.sin(orientation_rads)],
                      [np.sin(orientation_rads),np.cos(orientation_rads)]])
    t_sample = np.zeros((xy.shape[0],))
    for ii in range(xy.shape[0]):
        xy_centered = xy[ii,:] - screen_size / 2.0
        # x_rotate = xy_centered[0]*np.cos(orientation_rads) - xy_centered[1]*np.sin(orientation_rads)
        xy_rotate = np.matmul(xy_centered,a_inv.T)
        # t_sample[ii] = (x_rotate + screen_size[0] / 2.0) / speed * bin_rate
        t_sample[ii] = (xy_rotate[0] + screen_size[0] / 2.0) / speed * bin_rate
    return t_sample

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Analyse moving/oriented bar data.')
    parser.add_argument('experimentName', type=str, help='Experiment name (e.g. 20230802C)')
    parser.add_argument('-a','--algorithm', default='kilosort2', type=str, help='Sorting algorithm used (yass or kilosort2)')
    parser.add_argument('-p','--protocol', default='MovingChromaticBar', type=str, help='protocol name (e.g. MovingChromaticBar)')
    parser.add_argument('-g','--group', default=None, type=str, help='EpochGroup search string/substring (e.g. oriented)')
    parser.add_argument('-f','--file', nargs='+', default=None, type=str, help='Data file identifier(s) (e.g. data007 data008)')
    parser.add_argument('-c','--chunk_name', default=None, type=str, help='Name of chunk containing the file(s) (e.g., chunk1)')
    parser.add_argument('-b','--bin_rate', default=60.0, type=float, help='Bin rate for spikes')
    parser.add_argument('-r','--sample_rate', default=20000.0, type=float, help='Data sample rate')

    args = parser.parse_args()

    file_name = args.file
    experiment_name = args.experimentName
    bin_rate = args.bin_rate
    sort_algorithm = args.algorithm
    group = args.group
    protocol_name = args.protocol
    chunk_name = args.chunk_name

    # Get the output file path.
    SAVE_PATH = cfg.get_save_path()
    filepath = SAVE_PATH + args.experimentName + '_' + args.algorithm +  '_movingbar.mat'

    d = Dataset(experiment_name)

    param_names = ['orientation','contrast','speed','barSize']

    spike_dict, cluster_id, params, unique_params, pre_pts, stim_pts, tail_pts = d.get_spike_rate_and_parameters(
        protocolStr=protocol_name, groupStr=group, param_names=param_names, sort_algorithm=sort_algorithm, file_name=file_name, bin_rate=bin_rate, sample_rate=args.sample_rate)
    
    # Pull the spike time separately.
    spike_times, _, _, _, pre_pts_times, stim_pts_times, tail_pts_times = d.get_spike_times_and_parameters(
            protocolStr=protocol_name, groupStr=group, param_names=param_names, sort_algorithm=sort_algorithm, file_name=file_name, bin_rate=1000.0, sample_rate=args.sample_rate)
    
    # Get the receptive field parameters in microns.
    if chunk_name is not None:
        rf_params = d.get_rf_parameters(file_name = file_name, chunk_name = chunk_name, scale_type = 'microns')
        screen_x, screen_y = d.get_screen_size(file_name = file_name, chunk_name = chunk_name, scale_type = 'microns')

    # Get the spike counts for each cluster.
    psth = np.zeros((len(cluster_id), spike_dict[cluster_id[0]].shape[0], spike_dict[cluster_id[0]].shape[1]))
    for ii in range(len(cluster_id)):
        psth[ii,:,:] = spike_dict[cluster_id[ii]]
    
    params['barSize'] = np.array(params['barSize'])
    # Unique orientations.
    u_orientation = unique_params['orientation']
    # Unique contrasts
    u_contrast = unique_params['contrast']
    # Unique speeds.
    u_speed = unique_params['speed']
    # Unique bar sizes.
    u_width = unique_params['barSize']
    if np.ndim(u_width) > 1:
        u_width = u_width[:,0]
    else: 
        u_width = np.array([u_width[0]])
    
    avg_psth = np.zeros((len(cluster_id),len(u_contrast),len(u_orientation),len(u_speed), len(u_width), spike_dict[cluster_id[0]].shape[1]), dtype=np.float32)
    for count, cell in enumerate(spike_dict):
        spikes = spike_dict[cell]
        for ct_count, contrast in enumerate(u_contrast):
            for o_count, orientation in enumerate(u_orientation):
                for s_count, speed in enumerate(u_speed):
                    for b_count, barSize in enumerate(u_width):
                        idx = np.where((np.array(params['contrast']) == contrast) & (np.array(params['orientation']) == orientation) & (np.array(params['speed']) == speed) & (np.array(params['barSize'][:,0]) == barSize))[0]
                        avg_psth[count, ct_count, o_count, s_count, b_count, :] = np.mean(spikes[idx,:], axis=0)

    pre_pts = np.unique(pre_pts).astype(int)
    stim_pts = np.unique(stim_pts).astype(int)
    tail_pts = np.unique(tail_pts).astype(int)

    # Get the average spike rate during the full bar motion.
    if chunk_name is not None:
        avg_response = np.zeros((len(cluster_id),len(u_contrast),len(u_speed),len(u_orientation),len(u_width)), dtype=np.float32)
        for s_count, speed in enumerate(u_speed):
            avg_response[:,:,s_count,:,:] = get_response_from_rf_params(avg_psth[:,:,:,s_count,:,:],rf_params[:,:2],u_orientation,np.array([screen_x, screen_y]),speed,bin_rate,pre_pts[0])
        # avg_response = np.mean(avg_psth[:,:,:,pre_pts[0]:pre_pts[0]+stim_pts[0]], axis=3)
        # Reorder the axes.
        avg_response = np.transpose(avg_response, (0,1,2,4,3))
    else:
        avg_response = np.mean(avg_psth[:,:,:,:,:,pre_pts[0]:pre_pts[0]+stim_pts[0]], axis=5)
        # Reorder the axes.
        avg_response = np.transpose(avg_response, (0,1,3,4,2))

    # Save out a MAT file.
    analysis = Analysis()

    # Loop through and compute OSI and DSI for each cell/contrast.
    osi = np.zeros((len(cluster_id),len(u_contrast),len(u_speed),len(u_width)))
    osi_angle = np.zeros((len(cluster_id),len(u_contrast),len(u_speed),len(u_width)))
    dsi = np.zeros((len(cluster_id),len(u_contrast),len(u_speed),len(u_width)))
    dsi_angle = np.zeros((len(cluster_id),len(u_contrast),len(u_speed),len(u_width)))
    for cell in range(len(cluster_id)):
        for ct in range(len(u_contrast)):
            for sp in range(len(u_speed)):
                for wd in range(len(u_width)):
                    osi_global, osi_ang = analysis.compute_OSI(u_orientation, avg_response[cell,ct,sp,wd,:])
                    dsi_global, dsi_ang = analysis.compute_DSI(u_orientation, avg_response[cell,ct,sp,wd,:])
                    osi[cell,ct,wd] = osi_global
                    osi_angle[cell,ct,wd] = osi_ang
                    dsi[cell,ct,wd] = dsi_global
                    dsi_angle[cell,ct,wd] = dsi_ang

    mdic = {'avg_psth': avg_psth, 'avg_response': avg_response, 'cluster_id': cluster_id, 'u_contrast': u_contrast, 'u_orientation': u_orientation, 
        'osi': osi, 'osi_angle': osi_angle, 'dsi': dsi, 'dsi_angle': dsi_angle, 'bin_rate': args.bin_rate, 'barSize': params['barSize'], 
        'pre_pts': pre_pts, 'stim_pts': stim_pts, 'tail_pts': tail_pts, 'speed': np.array(params['speed']), 'psth':psth, 'contrast': np.array(params['contrast']),
        'orientation': np.array(params['orientation']), 'speed': np.array(params['speed']), 'width': np.array(params['barSize']),
        'spike_times': spike_times, 'pre_pts_times':pre_pts_times, 'stim_pts_times':stim_pts_times, 'tail_pts_times':tail_pts_times}
    
    if args.chunk_name is not None:
        mdic['rf_params'] = rf_params

    analysis.save_mat(filepath, mdic)
