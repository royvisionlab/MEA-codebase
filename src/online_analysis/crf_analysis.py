import os
import argparse
import numpy as np
from tqdm import tqdm
import visionloader as vl

def get_binned_spikes(vcd, cell_ids, frame_times, sample_rate: float=20000.0):
    '''
    Parameters:
        vcd: Vision data table object
        cells: list of cells to compute STAs for
    Returns:
        A matrix of size frames by cells

    Bins the spike train using the (full resolution) frames of the monitor 
    as bin edges.
    '''
    # Loop through cells and bin spike train.
    binned_spikes = []
    for cell in cell_ids:
        spike_times = vcd.get_spike_times_for_cell(cell) / sample_rate * 1000 # ms
        binned_spikes.append(np.histogram(spike_times,bins=frame_times)[0])
    return np.asarray(binned_spikes).T


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Spatiotemporal noise analysis')
    parser.add_argument('data_path', type=str, help='Name of experiment')
    parser.add_argument('file_name', type=str, help='Output file path (e.g., /path/to/file/out)')
    parser.add_argument('-a','--algorithm', default='vision', type=str, help='Sorting algorithm used (yass or kilosort*)')
    # parser.add_argument('-f','--file_name', nargs='+', default=None, type=str, help='Data file identifier')
    parser.add_argument('-s','--stim_time', default=20000, type=int, help='Duration of stimulus in ms')
    parser.add_argument('-p','--pre_time', default=250, type=int, help='Duration of pre-stimulus period in ms')
    parser.add_argument('-c','--contrasts', nargs='+', default=[5,10,20], type=float, help='Contrasts to analyze')
    parser.add_argument('-t','--temporal_frequency', default=4.0, type=float, help='Temporal frequency of stimulus')

    args = parser.parse_args()
    
    data_path = args.data_path
    file_name = args.file_name
    contrasts = args.contrasts
    pre_time = args.pre_time
    stim_time = args.stim_time
    temporal_frequency = args.temporal_frequency
    

    sample_rate = 20000.0
    frame_rate = 59.941548817817917 # Upper bound on the frame rate to make sure that we don't miss any frames.

    # data_path = '/gscratch/retina/data/sorted/20241016C/data009/vision'
    # file_name = 'vision'
    # contrasts = [5,10,20] # -c
    # pre_time = 250 # -p
    # stim_time = 20000 # -s
    # temporal_frequency = 4.0 # -f

    cycle_bins = np.floor(60.0/temporal_frequency).astype(int) # Number of bins in a cycle
    pre_frames = np.floor(pre_time/1000*frame_rate).astype(int) # Number of frames in the pre-stimulus period
    stim_frames = np.floor(stim_time/1000*frame_rate).astype(int) # Number of frames in the stimulus

    vcd = vl.load_vision_data( data_path, file_name, include_neurons=True )
    cell_ids = sorted(vcd.get_cell_ids())
    epoch_starts = vcd.get_ttl_times()

    # Define the bin edges for the full resolution frames.
    frame_times = np.arange(0,stim_frames+pre_frames+1,1).astype(float)/frame_rate*1000.0

    binned_spikes = np.zeros((len(epoch_starts),len(cell_ids),len(frame_times)-1))
    for ii in range(len(epoch_starts)):
        binned_spikes[ii,...] = get_binned_spikes(vcd, cell_ids, frame_times+float(epoch_starts[ii])/sample_rate*1000.0, sample_rate).T

    binned_spikes *= frame_rate # Convert to Hz
    
    # Get the average binned spikes for the unique contrasts.
    contrasts = np.array(contrasts)
    u_contrasts = np.unique(contrasts)
    if len(u_contrasts) != len(contrasts):
        avg_psth = np.zeros((len(u_contrasts),binned_spikes.shape[1],binned_spikes.shape[2]))
        for ii, contrast in enumerate(u_contrasts):
            idx = np.where(contrasts == contrast)[0]
            avg_psth[ii,...] = np.nanmean(binned_spikes[idx,...],axis=0)
        binned_spikes = avg_psth
    
    # Compute the F1 modulation for each cell.
    F1 = np.zeros((binned_spikes.shape[0],binned_spikes.shape[1]))
    for ii in range(binned_spikes.shape[0]):
        for jj in range(binned_spikes.shape[1]):
            spike_train = binned_spikes[ii,jj,pre_frames:pre_frames+stim_frames]
            spike_train = spike_train[:len(spike_train)//cycle_bins*cycle_bins].reshape(-1,cycle_bins).T
            avg_cycle = np.nanmean(spike_train,axis=1)
            try:
                ft = np.fft.fft(avg_cycle)
                F1[ii,jj] = np.abs(ft)[1]/float(len(avg_cycle))*2.0
            except Exception as e:
                print('Error computing F1 for cell: ', jj, ' in epoch: ', ii)
                print(e)
                F1[ii,jj] = np.max(avg_cycle) - np.min(avg_cycle)

    # Find the top 5% of cells with the highest F1 values
    top_cells = np.argsort(np.nanmean(F1,axis=0))[::-1][:int(len(cell_ids)*0.05)]
    print('F1 modulation for top 5% of cells: ', np.nanmean(F1[:,top_cells], axis=1))
    
    # Find the top 10% of cells with the highest F1 values
    top_cells = np.argsort(np.nanmean(F1,axis=0))[::-1][:int(len(cell_ids)*0.1)]
    print('F1 modulation for top 10% of cells: ', np.nanmean(F1[:,top_cells], axis=1))

    # Log the results to a text file.
    with open(data_path + file_name + '.log', 'w') as f:
        f.write('F1 modulation for top 5% of cells: ' + str(np.nanmean(F1[:,top_cells], axis=1)) + '\n')
        f.write('F1 modulation for top 10% of cells: ' + str(np.nanmean(F1[:,top_cells], axis=1)) + '\n')
        f.write('Contrasts: ' + str(contrasts) + '\n')
        f.write('Temporal frequency: ' + str(temporal_frequency) + '\n')
        f.write('Pre-stimulus time: ' + str(pre_time) + '\n')
        f.write('Stimulus time: ' + str(stim_time) + '\n')
        f.write('Cycle bins: ' + str(cycle_bins) + '\n')
        f.write('Pre frames: ' + str(pre_frames) + '\n')
        f.write('Stim frames: ' + str(stim_frames) + '\n')
        f.write('Frame rate: ' + str(frame_rate) + '\n')
        f.write('Sample rate: ' + str(sample_rate) + '\n')
        f.write('Epoch starts: ' + str(epoch_starts) + '\n')
        f.write('F1: ' + str(F1) + '\n')
        f.write('Top cells: ' + str(top_cells) + '\n')

