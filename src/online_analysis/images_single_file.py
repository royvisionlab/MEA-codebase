
import numpy as np
import visionloader as vl
import config as cfg
import argparse
from symphony_data import Dataset
import lnp
import pickle
import hdf5storage

# Minimal flashed images analysis.
class ImageAnalysisSingle(object):
    def __init__(self, 
            data_path: str,
            file_name: str,
            sort_algorithm: str='kilosort2'):
        self.sort_algorithm = sort_algorithm
        # self.vcd = vl.load_vision_data( os.path.join(data_path,file_name), file_name, include_neurons=True )
        self.vcd = vl.load_vision_data( data_path, file_name, include_neurons=True )
        # Get the cell Ids.
        self.cell_ids = sorted(self.vcd.get_cell_ids())

    def analyze(self, 
                    parameters: dict,
                    seed: int,
                    num_images: int,
                    frame_rate: float=60.31807657):
        
        epoch_starts = self.vcd.get_ttl_times()

        # Loop through each epoch.
        for count, start in enumerate(epoch_starts):
            
            # Need to get rid of pre/post frames.
            preFrames = round(parameters['preTime'] * 1e-3 * frame_rate)
            stimFrames = round(parameters['stimTime'] * 1e-3 * frame_rate)

            frame_times = np.round(start/20000*1000 + np.arange(stimFrames+1)*1000/frame_rate)

            # Get the frames when the stimulus was on.
            idx = (frame_times-frame_times[0]).searchsorted((parameters['preTime']+parameters['stimTime']),'right')
            frame_times = frame_times[preFrames:idx+1]

            # Get the binned spike count.
            binned_spikes = self.get_binned_spikes(frame_times, stride=1)
            if count == 0:
                psth = np.zeros((len(epoch_starts), binned_spikes.shape[0], binned_spikes.shape[1]))
            psth[count,:,:] = binned_spikes
        
        # Swap the axes.
        # psth = psth.transpose((1,0,2))

        # Get the image sequence.
        num_reps = np.ceil(float(len(epoch_starts))/float(num_images))
        sequence = list()
        for jj in range(num_reps):
            sequence.append(np.arange(0,num_images))
        sequence = np.asarray(sequence).flatten()

        if seed > 0:
            np.random.seed(seed)
            idx = np.floor(np.random.rand(len(sequence))*len(sequence)).astype(int)
            sequence = sequence[idx]

        # Get the unique image sequence.
        avg_psth = np.zeros((num_images, psth.shape[1], psth.shape[2]))
        sd_psth = np.zeros((num_images, psth.shape[1], psth.shape[2]))
        for ii in range(num_images):
            idx = np.where(sequence==ii)[0]
            if len(idx) > 0:
                tmp_data = psth[idx,:,:]
                avg_psth[ii,:,:] = np.nanmean(tmp_data, axis=0)
                sd_psth[ii,:,:] = np.nanstd(tmp_data, axis=0)
        avg_psth = avg_psth.transpose((1,0,2))
        sd_psth = sd_psth.transpose((1,0,2))
        return avg_psth, sd_psth

    def get_binned_spikes(self, frame_times, stride=1, sample_rate: float=20000.0):
        '''
        Parameters:
            vcd: Vision data table object
            cells: list of cells to compute STAs for
        Returns:
            A matrix of size frames by cells

        Bins the spike train using the (full resolution) frames of the monitor 
        as bin edges.
        '''

        # If the stride/binsPerFrame is greater than 1, interpolate the frame times.
        if stride > 1:
            frame_times = lnp.get_frame_times(frame_times, stride)

        # Loop through cells and bin spike train.
        binned_spikes = []

        for cell in self.cell_ids:
            spike_times = self.vcd.get_spike_times_for_cell(cell) / sample_rate * 1000 # ms
            binned_spikes.append(np.histogram(spike_times,bins=frame_times)[0])
        
        return np.asarray(binned_spikes).T

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Camo break analysis')
    parser.add_argument('experimentName', type=str, help='Name of experiment')
    parser.add_argument('-a','--algorithm', default='kilosort2', type=str, help='Sorting algorithm used (yass or kilosort2)')
    parser.add_argument('-p','--protocol', default='PresentImages', type=str, help='Name of stimulus protocol to analyze')
    parser.add_argument('-f','--file', nargs='+', default=None, type=str, help='Data file identifier')
    parser.add_argument('-s','--seed', default=None, type=int, help='Name of chunk containing the file(s)')
    parser.add_argument('-n','--num_images', default=None, type=int, help='Name of chunk containing the file(s)')
    parser.add_argument('-b','--bin_rate', default=60.0, type=float, help='Bin rate for spikes')
    parser.add_argument('-r','--sample_rate', default=20000.0, type=float, help='Data sample rate')

    args = parser.parse_args()

    # Get the output file path.
    SORT_PATH, _ = cfg.get_data_paths()
    filepath = SORT_PATH + args.experimentName + '/' + args.experimentName + '_' + args.algorithm +  '_PresentImages'

    parameters = dict()
    parameters['preTime'] = 250.0
    parameters['stimTime'] = 250.0
    parameters['tailTime'] = 250.0

    a = ImageAnalysisSingle(args.data_path, args.file, sort_algorithm=args.algorithm)

    avg_psth, sd_psth = a.analyze(parameters=parameters, seed=args.seed, num_images=args.num_images)

    mdic = {'avg_psth': avg_psth, 'sd_psth': sd_psth, 'parameters': parameters}

    # Save the results.
    with open(filepath+'.pkl','wb') as f:
        pickle.dump(mdic, f, protocol=pickle.HIGHEST_PROTOCOL)