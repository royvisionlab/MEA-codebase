# import struct
# import sys
import os
# from more_itertools import tail
# sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import math
# import platform
import json
import lnp
from scipy.io import savemat
from scipy.signal.windows import gaussian
# from scipy.signal import wiener, filtfilt, butter, freqz
from scipy.ndimage import filters, gaussian_filter
import pickle
import hdf5storage
from sys import getsizeof
from collections import Counter
import warnings
from typing import Union
import config as cfg # Path configuration
from tqdm import tqdm

import visionloader as vl
import visionwriter as vw

os.environ['KMP_DUPLICATE_LIB_OK']='True'


SPIKE_TIMES_FILENAME = 'spike_times.npy'
SPIKE_IDENTITY_FILENAME = 'spike_clusters.npy'
CLUSTER_QUALITY_FILENAME = 'cluster_KSLabel.tsv'
SPIKE_LOCATIONS_FILENAME = 'spike_xy.npy'

def psth(x, bin_rate: float=100.0, filter_sigma: float=6.0+2.0/3.0):
    b = gaussian(39, filter_sigma/1000*bin_rate)
    b /= b.sum()
    y = np.zeros(x.shape)
    for row in range(x.shape[0]):
        y[row,:] = filters.convolve1d(x[row,:], b)
    return y

class Metadata(object):

    def __init__(self, 
                experimentName: str):
        """
        Constructor method
        Parameters:
            experimentName: name of the experiment (e.g. '20230105C')
        """

        # Set the path to the json files according to OS.
        sortpath, datapath = cfg.get_data_paths()
        
        self.datapath = datapath
        self.sortpath = sortpath
        self.sampleRate = 20000.0
        self.metadata = None

        # Load the json file.
        if os.path.exists(os.path.join(datapath, experimentName + '.json')):
            self.loadJSON(os.path.join(datapath, experimentName + '.json'))
        
        if self.metadata is None:
            print('No metadata found for experiment: ' + experimentName)

    def loadJSON(self, filepath: str):
        """
        Loads the json file and stores it in the metadata attribute.
        Parameters:
            filepath: path to the json file.
        """
        # Read the JSON
        with open(filepath) as f:
            self.metadata = json.load(f)
    
    def export_json(self, filepath):
        """
        Exports the metadata to a text file.
        Parameters:
            filepath: path to the json file.
        """
        out_path = filepath.replace('.json','.txt')
        self.loadJSON(filepath=filepath)
        with open(out_path, 'w') as file:
            for protocol in self.metadata['protocol']:
                protocol,_ = self.check_protocol_format(protocol)
                file.write(protocol['label'] + '\n')
                for group in protocol['group']:
                    file.write('  -' + group['label'] + '\n')
                    for block in group['block']:
                        file.write('    ->' + self.get_data_file_from_block(block) + '\n')

    # Get the name of the Litke data file saved for the EpochBlock
    def get_data_file_from_block(self, block: dict):
        """
        Extracts the data file name from the block dictionary.
        Parameters:
            block: Block dictionary (type: dict)
            Returns: Data file name (type: str)
        """
        dataFile = block['dataFile'].split('/')[-2:-1]
        dataFile = dataFile[0]
        return dataFile

    def get_stimulus_parameters(self, protocol: dict, param_names: list):
        """
        Arguments:
            protocol: Protocol dictionary (type: dict)
            param_names: List of parameters to search on (type: list)
        Returns: Dictonary with the parameter values and unique parameter values.
        """
        params = dict()
        # Initialize the dictionary.
        for value in param_names:
            params[value] = list()

        for group in protocol['group']:
            for block in group['block']:
                for epoch in block['epoch']:
                    parameters = epoch['parameters']
                    for value in param_names:
                        if value in parameters:
                            params[value].append(parameters[value])
        
        # Get the unique parameters.
        unique_params = dict()
        for key, value in params.items():
            unique_params[key] = self.unique(value)

        return params, unique_params

    def check_type(self, obj):
        """
        Checks if the object is a list of strings.
        Parameters:
            obj: Object to check (type: list)
        Returns: True if it's a list of strings, False otherwise.
        """
        return bool(obj) and all(isinstance(elem, str) for elem in obj)

    def unique(self, obj: list):
        """
        Returns the unique elements of a list.
        Parameters: 
            obj: List of elements (type: list)
        Returns: List of unique elements (type: list)
        """
        # Check if it's a list of strings.
        if self.check_type(obj):
            return list(Counter(obj).keys())
        else: 
            return np.unique(np.array(obj))

    def get_object_keys_as_list(self, obj):
        """
        Returns the keys of a dictionary as a list.
        Parameters:
            obj: Dictionary (type: dict)
            
        Returns: List of keys (type: list)
        """
        return list(obj.keys())

    def get_frame_times_ms(self, block):
        """
        Returns the frame times in ms for each epoch.
        Parameters:
            block: Block dictionary (type: dict)

        Returns: Dictionary with the frame times for each epoch (type: dict)
        """
        epoch_frames = dict()

        # Grab the epoch starts and convert to ms.
        if isinstance(block['epochStarts'], list):
            epoch_start = np.zeros((len(block['epochStarts']),))
            count = 0
            for start in block['epochStarts']:
                epoch_start[count] = round(start / self.sampleRate * 1000)
                count += 1
        else:
            epoch_start = list()
            epoch_start.append(block['epochStarts'])
        
        # Now, get the frame times in msec.
        count = 0
        for f in block['frameTimesMs']:
            if count < len(epoch_start):
                f = np.array(f).astype(float)
                # Add the epoch start to get the frame timing re to the MEA DAQ clock.
                f = f + epoch_start[count]
                # Add it to the dictionary.
                epoch_frames[count] = f
            count += 1
        
        return epoch_frames

    def search_data_file(self, 
                        protocol_id: str,
                        file_name: Union[str, list]):
        """
        Searches for a file or list of files associated with a particular stimulus protocol.

        Parameters:
            protocol_id: name of protocol to search for (e.g. 'gratings')
            file_name: name of file to search for (e.g. 'data003'); can be a list of files
        Returns:
            protocol: dictionary containing the protocol information
        """

        if type(file_name) is str:
            file_name = [file_name]

        protocol = self.searchProtocol(protocol_id)
        
        if protocol is None:
            print('Warning: Failed to find protocol ' + protocol_id)
            return protocol

        files_found = False
        pos_group = list()
        pos_block = list()
        for group_idx, g in enumerate(protocol['group']):
            for block_idx, b in enumerate(g['block']):
                if self.get_data_file_from_block(b) in file_name:
                    files_found = True
                    pos_group.append(group_idx)
                    pos_block.append(block_idx)
                
        if not files_found:
            protocol = None
            print('Warning: Failed to find specified data files associated with protocol ' + protocol_id)
            return protocol
        

        new_protocol = dict()
        groups = list()
        for group_idx, g in enumerate(protocol['group']):
            if (group_idx in pos_group):
                blocks = list()
                for block_idx, b in enumerate(g['block']):
                    if (group_idx in pos_group and block_idx in pos_block):
                        blocks.append(b)
                new_group = dict()
                new_group['label'] = g['label']
                new_group['block'] = blocks
                groups.append(new_group)
        new_protocol['label'] = protocol['label']
        new_protocol['group'] = groups

        new_protocol, _ = self.check_protocol_format(new_protocol)
        return new_protocol

    def searchProtocol(self, 
                        protocolStr: str, 
                        groupStr: str=None):
        """
        Searches for a protocol and group based on a search string.
        Parameters:
            protocolStr: String to search for in the protocol label (type: str)
            groupStr: String to search for in the group label (type: str)

        Returns: Dictionary with the protocol and group information (type: dict)
        """
        protocol = None

        # Search through the protocols and find protocol labels that match the search string.
        for p in self.metadata['protocol']:
            if (p['label'] != None and protocolStr.lower() in p['label'].lower()):
                if protocol is None:
                    protocol = p
                else:
                    protocol['group'].extend(p['group'])

        # The code is expecting a list of groups, so if there's only one group, you need to turn it into a list...
        if isinstance(protocol['group'], dict):
            g = list()
            g.append(protocol['group'])
            protocol['group'] = g

        # Search within the groups for the group string.
        if groupStr != None:
            group = list()
            for g in protocol['group']:
                if (g['label'] != None and groupStr.lower() in g['label'].lower()):
                    group.append(g)
            protocol['group'] = group

        # Check the data structure.
        self.check_protocol_format(protocol)

        return protocol
    
    def searchGroup(self, group, protocolStr):
        return None

    def check_protocol_format(self, protocol):
        """
        Checks the format of the protocol dictionary and makes sure it's in the expected format.
        Parameters:
            protocol: Dictionary containing the protocol information (type: dict)
        
        Returns: Dictionary containing the protocol information (type: dict)
        """
        data_files = list()

        # Make sure you have lists of dicts and not just dicts.
        if isinstance(protocol['group'], dict):
            g = list()
            g.append(protocol['group'])
            protocol['group'] = g

        g = list()
        for group in protocol['group']:
            if ('block' in group):
                if isinstance(group['block'], dict):
                    b = list()
                    b.append(group['block'])
                    group['block'] = b

                # Grab the epoch blocks
                b = list()
                for block in group['block']:
                    dataFile = block['dataFile'] #self.get_data_file_from_block(block)
                    if len(dataFile) > 0:
                        data_files.append(dataFile)
                        if isinstance(block['epoch'], dict):
                            e = list()
                            e.append(block['epoch'])
                            block['epoch'] = e
                        b.append(block)
                if len(b) > 0:
                    group['block'] = b
                    g.append(group)
        protocol['group'] = g

        return protocol, data_files


class Dataset(object):

    def __init__(self, 
            experimentName: str,
            sort_algorithm: str='kilosort2',
            stim_time_name: str='stimTime'):
        """
        Constructor method
        Parameters:
            datapath: dataset name
            wn_datarun: white noise datarun
            estim_datarun: estim data run
            movie_xml_str: white noise movie XML
        """

        self.vcd = None
        self.cell_ids = None
        self.experimentName = experimentName
        self.sort_algorithm = sort_algorithm
        self.stim_time_name = stim_time_name
        self.M = Metadata(experimentName)
  
        # Make sure the metadata is loaded.
        if self.M.metadata is None:
            raise Exception('Metadata not loaded')

    def parameters_to_dict(self, parameters: dict, out_dict: dict=None):
        """
        Converts a dictionary of parameter to a dictionary that can be saved to an analysis format (.mat, .pkl).
        Parameters:
            parameters: dictionary of parameters
        Returns:
            out_dict: dictionary of parameters
        """
        if out_dict is None:
            out_dict = dict()
        for param_name in parameters.keys():
            if isinstance(parameters[param_name],list) and len(parameters[param_name]) > 0:
                if isinstance(parameters[param_name][0], str):
                    out_dict[param_name] = parameters[param_name] #np.array(parameters[param_name], dtype='U')
                else:
                    out_dict[param_name] = np.array(parameters[param_name])
            else:
                out_dict[param_name] = parameters[param_name]
        return out_dict
    
    def get_rf_parameters(self, file_name: Union[str, list]=None, chunk_name: str=None, scale_type: str='microns') -> np.ndarray:
        """ Pull the receptive field parameters from the vision table for the specified chunk. 
        Parameters:
            :param file_name: Name of the data file (type: str or list of strings)
            :param chunk_name: Name of the chunk (type: str)
            :param scale_type: Type of scale to use for the receptive field parameters (type: str, default: 'microns')
        Returns: 
            :returns: Array containing the receptive field parameters [x0,y0,sigma_x,sigma_y,theta] (type: np.ndarray)
            :rtype: np.ndarray
        """

        if file_name is None or chunk_name is None:
            raise Exception('Must specify both a file_name and chunk_name')

        if type(file_name) is str:
            file_name = [file_name]

        # Load the vision table for the chunk.
        filepath = os.path.join(self.M.sortpath, self.experimentName, chunk_name, self.sort_algorithm)
        if not os.path.exists(filepath):
            print('Warning: Could not find chunk ' + chunk_name + ' for file ' + file_name + ' in ' + filepath + '. Trying to load from analysis directory.')
            filepath = os.path.abspath(os.path.join(self.M.sortpath,'../../analysis/', self.experimentName, chunk_name, self.sort_algorithm))
            if not os.path.exists(filepath):
                print('Warning: Could not find chunk ' + chunk_name + ' for file ' + file_name + ' in ' + filepath + '.')
                return None
        try:
            vcd = vl.load_vision_data(filepath, self.sort_algorithm, include_params=True, include_runtimemovie_params=True)
            
            # Get the scale factors.
            if scale_type == 'microns':
                scale_factor = vcd.runtimemovie_params.micronsPerStixelX
            else:
                scale_factor = float(vcd.runtimemovie_params.pixelsPerStixelX)
            
            # Get the cell Ids.
            chunk_ids = sorted(vcd.get_cell_ids()) # sort by convention!
            rf_fits = [vcd.get_stafit_for_cell(cell_id) for cell_id in chunk_ids]
            rf_fits = np.array(rf_fits)

            # Multiply the receptive field parameters by the scale factor.
            rf_fits[:,:4] *= scale_factor

            # Get the unique cell ids from the target data files.
            ids = list()
            for f in file_name:
                filepath = os.path.join(self.M.sortpath, self.experimentName, f, self.sort_algorithm)
                vcd = vl.load_vision_data(filepath, f, include_neurons=True)
                ids.extend(vcd.get_cell_ids())

            chunk_ids = np.unique(np.array(chunk_ids))
            ids = np.unique(np.array(ids))

            # Get the parameters in the correct order.
            rf_params = np.zeros((len(ids), 5))
            for i, id in enumerate(ids):
                idx = np.where(chunk_ids == id)[0]
                if len(idx) > 0:
                    rf_params[i,:] = rf_fits[idx,:]
        except Exception as error:
            print('An error occurred while loading the STA fits:', type(error).__name__, '-', error)
            return None
        
        return rf_params
    
    def get_cell_classifications(self, chunk_name: str) -> dict:
        """ Pull the cell classifications from the vision table for the specified chunk. 
        Arguments:
            :param chunk_name: Name of the chunk (type: str)
            :type chunk_name: str
        Returns:
            :return: Dictionary containing the cell classifications (type: dict)
            :rtype: dict"""
        # Load the vision table for the chunk.
        filepath = os.path.join(self.M.sortpath, self.experimentName, chunk_name, self.sort_algorithm)
        vcd = vl.load_vision_data(filepath, self.sort_algorithm, include_params=True)

        # cell_ids = sorted(vcd.get_cell_ids())
        # cell_types = [vcd.get_cell_type_for_cell(cell_id) for cell_id in cell_ids]
        cell_types = vcd.get_all_present_cell_types()
        # Create an ouput dictionary to hold cell id by type.
        cell_ids_by_type = dict()
        for cell_type in cell_types:
            try:
                cell_ids_by_type[cell_type] = vcd.get_all_cells_of_type(cell_type)
            except Exception as error:
                print('An error occurred while loading the cell types:', type(error).__name__, '-', error)
        return cell_ids_by_type

    def get_screen_size(self, file_name: Union[str, list]=None, chunk_name: str=None, scale_type: str='microns') -> np.ndarray:
        if file_name is None or chunk_name is None:
            raise Exception('Must specify both a file_name and chunk_name')

        if type(file_name) is str:
            file_name = [file_name]

        # Load the vision table for the chunk.
        filepath = os.path.join(self.M.sortpath, self.experimentName, chunk_name, self.sort_algorithm)
        vcd = vl.load_vision_data(filepath, self.sort_algorithm, include_params=True, include_runtimemovie_params=True)
        # Get the scale factors.
        if scale_type == 'microns':
            scale_factor = vcd.runtimemovie_params.micronsPerStixelX
        else:
            scale_factor = float(vcd.runtimemovie_params.pixelsPerStixelX)
        
        screen_x = vcd.runtimemovie_params.width * scale_factor
        screen_y = vcd.runtimemovie_params.height * scale_factor
        return screen_x, screen_y

    def get_interspike_interval(self, 
                    protocol_id: str=None,
                    group_str: str=None,
                    sort_algorithm: str='kilosort2',
                    file_name: Union[str, list]=None,
                    bin_edges: np.array=np.linspace(0,300,601)):
        """
        Compute the interspike interval histogram for a given protocol.
        Parameters:
            protocol_id: Name of stimulus protocol (type: str)
            group_str: Search string for the group label (type: str); default=None
            sort_algorithm: spike sorting algorithm (type: str); default='yass'
            file_name: file name(s) to search for (type: str or list); default=None
            bin_edges: bin edges for the histogram (type: np.array); default=np.linspace(0,300,601)
        """

        if file_name is not None:
            protocol = self.M.search_data_file(protocol_id, file_name=file_name)
        else: 
            protocol = self.M.searchProtocol(protocol_id, group_str)
        
        acf_dict = dict() # Autocorrelation dictionary

        # Get the number of groups.
        nGroups = len(protocol['group'])
        # Loop through the EpochGroups
        for groupCount, group in enumerate(protocol['group']):

            nBlocks = len(group['block'])
            # Grab the epoch blocks
            for blockCount, block in enumerate(group['block']):
                print('Processing group ' + str(groupCount+1) + ' of ' + str(nGroups) + ' and block ' + str(blockCount+1) + ' of ' + str(nBlocks))

                dataFile = self.M.get_data_file_from_block(block)

                # Try to load a data file.
                filepath = os.path.join(self.M.sortpath, self.experimentName, dataFile, sort_algorithm)

                self.load_vision_table(filepath, dataFile)

                # Compute the autocorrelation function.
                acf_tmp = self.get_interspike_intervals(bin_edges=bin_edges)

                for cell in self.cell_ids:
                    if cell in acf_dict:
                        acf_dict[cell] += acf_tmp[cell].astype(int)
                    else:
                        acf_dict[cell] = acf_tmp[cell].astype(int)

        # Sort the ACF dictionary by cluster Id.
        acf_dict = dict(sorted(acf_dict.items()))
        cluster_id = list(acf_dict.keys())

        isi = np.zeros((len(cluster_id), acf_dict[cluster_id[0]].shape[0]))
        acf = np.zeros((len(cluster_id), acf_dict[cluster_id[0]].shape[0]))
        for count, cell in enumerate(acf_dict):
            isi[count,:] = acf_dict[cell]
            if np.sum(isi[count,:]) > 1:
                acf[count,:] = isi[count,:] / np.sum(isi[count,:])
            else:
                acf[count,:] = isi[count,:]

        return acf, isi, np.array(cluster_id)

    def get_frame_times(self, 
                    protocolStr: str=None,
                    groupStr: str=None,
                    file_name: Union[str, list]=None) -> list:
        if protocolStr == None:
            print('Protocol ID is empty. Cannot continue...')
            return None

        if file_name is not None:
            protocol = self.M.search_data_file(protocolStr, file_name=file_name)
        else: 
            protocol = self.M.searchProtocol(protocolStr, groupStr)

        # Get the spike count dictionary.
        frame_times = list()
        for group in protocol['group']:
            for block in group['block']:
                epoch_frames = self.M.get_frame_times_ms( block )
                for count, epoch in tqdm(enumerate(block['epoch']), desc="Pulling frame times", total=len(block['epoch'])):
                    if count < len(epoch_frames):
                        frame_times.append(epoch_frames[count])
        return frame_times

    def build_cluster_quality_dict(self, filepath):
        cluster_quality_by_id = {}
        with open(filepath, 'r') as cluster_quality_file:
            cluster_quality_file.readline()

            remaining_lines = cluster_quality_file.readlines()

            for line in remaining_lines:
                data_list = line.strip('\n').split('\t')
                cluster_quality_by_id[int(data_list[0])+1] = data_list[1]

        return cluster_quality_by_id

    def get_spike_times_yass(self, filepath):
        spikes = np.load(filepath + '.npy')
        spike_times_vector = spikes[:,0]
        spike_identity_vector = spikes[:,1]
        n_spikes = spike_times_vector.shape[0]
        cell_ids = []
        spikes_by_cell_id = {}
        for i in range(n_spikes):
            spike_time = spike_times_vector[i]
            spike_id = spike_identity_vector[i]  + 1
            # we add 1 because Vision/MATLAB requires that real cells start at index 1
            # (MATLAB does 1-based indexing)
            if spike_id not in spikes_by_cell_id:
                spikes_by_cell_id[spike_id] = []
                cell_ids.append(spike_id)
            spikes_by_cell_id[spike_id].append(spike_time)
        cell_ids = np.unique(cell_ids)
        return spikes_by_cell_id, cell_ids
    
    def get_spike_times(self, filepath):
        spike_times_filepath = os.path.join(filepath, SPIKE_TIMES_FILENAME)
        spike_identity_filepath = os.path.join(filepath, SPIKE_IDENTITY_FILENAME)
        cluster_quality_dict_filepath = os.path.join(filepath, CLUSTER_QUALITY_FILENAME)

        spike_times_vector = np.load(spike_times_filepath)
        spike_identity_vector = np.load(spike_identity_filepath)
        quality_by_cluster_id = self.build_cluster_quality_dict(cluster_quality_dict_filepath)

        n_spikes = spike_times_vector.shape[0]

        goodClusters = []
        spikes_by_cell_id = {}
        for i in range(n_spikes):
            spike_time = spike_times_vector[i,0]
            spike_id = spike_identity_vector[i,0]  + 1
            # we add 1 because Vision/MATLAB requires that real cells start at index 1
            # (MATLAB does 1-based indexing)

            if quality_by_cluster_id[spike_id] == 'good':
                goodClusters.append(spike_id)
                if spike_id not in spikes_by_cell_id:
                    spikes_by_cell_id[spike_id] = []

                spikes_by_cell_id[spike_id].append(spike_time)
        
        cell_ids = np.unique(goodClusters)
        return spikes_by_cell_id, cell_ids

    def get_spike_rate_and_parameters(self, 
                    protocolStr: str=None,
                    groupStr: str=None,
                    param_names: list=[],
                    sort_algorithm: str='kilosort2',
                    file_name: Union[str, list]=None,
                    bin_rate: float=1000.0,
                    sample_rate: float=20000.0):
        """
        Get the binned spike count and desired epoch parameters
        Parameters:
            protocolStr: Substring within the protocol ID you want to search for. (Required)
            groupStr: Substring with the epoch group label to search for (optional).
            param_names: List of parameter names to search for.
            sort_algorith: Name of algorithm used for spike detection ('yass' or 'kilosort2')
            bin_rate: Bin rate for the spike count in Hz.
            sample_rate: Sample rate for the data in Hz (default=20000.0)
        """

        spike_dict, cluster_id, params, unique_params, pre_pts, stim_pts, tail_pts = self.get_spike_count_and_parameters(
            protocolStr, groupStr, param_names, sort_algorithm, file_name, bin_rate, sample_rate)

        # Convert from spike count to spike rate.
        for cell in spike_dict:
            spike_dict[cell] = psth(spike_dict[cell]*bin_rate, bin_rate)        

        return spike_dict, np.array(cluster_id), params, unique_params, pre_pts, stim_pts, tail_pts
    
    def get_spike_count_and_parameters(self, 
                    protocolStr: str=None,
                    groupStr: str=None,
                    param_names: list=[],
                    sort_algorithm: str='kilosort2',
                    file_name: Union[str, list]=None,
                    bin_rate: float=1000.0,
                    sample_rate: float=20000.0):
        """
        Get the binned spike count and desired epoch parameters
        Parameters:
            protocolStr: Substring within the protocol ID you want to search for. (Required)
            groupStr: Substring with the epoch group label to search for (optional).
            param_names: List of parameter names to search for.
            sort_algorith: Name of algorithm used for spike detection ('yass' or 'kilosort2')
            bin_rate: Bin rate for the spike count in Hz.
            sample_rate: Sample rate for the data in Hz (default=20000.0)
        """

        if protocolStr == None:
            print('Protocol ID is empty. Cannot continue...')
            return None

        if file_name is not None:
            protocol = self.M.search_data_file(protocolStr, file_name=file_name)
        else: 
            protocol = self.M.searchProtocol(protocolStr, groupStr)

        spike_dict = dict()

        # Take a quick pass through the data files.
        if not param_names:
            params = dict()
            unique_params = dict()
        else:
            params, unique_params = self.M.get_stimulus_parameters(protocol, param_names)

        # Get the spike count dictionary.
        spike_dict, cluster_id, pre_pts, stim_pts, tail_pts = self.get_spike_count_dict(protocol, 
                                        sort_algorithm=sort_algorithm,
                                        bin_rate=bin_rate,
                                        sample_rate=sample_rate)
        
        return spike_dict, np.array(cluster_id), params, unique_params, pre_pts, stim_pts, tail_pts

    def get_count_at_frame_multiple(self, 
                    protocolStr: str=None,
                    groupStr: str=None,
                    param_names: list=[],
                    sort_algorithm: str='kilosort2',
                    file_name: Union[str, list]=None,
                    frame_rate: float=60.31807657,
                    stride: int=1):
        """
        Get the binned spike count and desired epoch parameters
        Parameters:
            protocolStr: Substring within the protocol ID you want to search for. (Required)
            groupStr: Substring with the epoch group label to search for (optional).
            param_names: List of parameter names to search for.
            sort_algorith: Name of algorithm used for spike detection ('yass' or 'kilosort2')
            file_name: Name of the file to search for. (string or list of strings)
            frame_rate: Frame rate for the stimulus in Hz.
            stride: Number of bins per frame (default=1)
        Returns:
            spike_dict: Dictionary of spike counts for each cell.
            cluster_id: List of cluster IDs for each cell.
            params: Dictionary of stimulus parameters.
            unique_params: Dictionary of unique stimulus parameters.
            pre_pts: Number of pre-stimulus points.
            stim_pts: Number of stimulus points.
            tail_pts: Number of post-stimulus points.
            mean_frame_rate: Mean frame rate for the stimulus.
        """

        if protocolStr == None:
            print('Protocol ID is empty. Cannot continue...')
            return None

        if file_name is not None:
            protocol = self.M.search_data_file(protocolStr, file_name=file_name)
        else: 
            protocol = self.M.searchProtocol(protocolStr, groupStr)

        bin_rate = frame_rate*float(stride)

        spike_dict = dict()
        num_epochs, time_bins = self.get_data_size(protocol, bin_rate)
        
        # Create an empty numpy array to hold the data.
        pre_pts = np.zeros((num_epochs,))
        stim_pts = np.zeros((num_epochs,))
        tail_pts = np.zeros((num_epochs,))

        # Take a quick pass through the data files.
        if not param_names:
            params = dict()
            unique_params = dict()
        else:
            params, unique_params = self.M.get_stimulus_parameters(protocol, param_names)

        # Get the spike count dictionary.
        total_epoch_count = 0
        stimulus = Stimulus()
        frame_times = list()
        diff_frames = list()
        for group in protocol['group']:
            for block in group['block']:
                dataFile = self.M.get_data_file_from_block(block)

                # Try to load a data file.
                filepath = os.path.join(self.M.sortpath, self.experimentName, dataFile, sort_algorithm)
                
                self.load_vision_table(filepath, dataFile)
                epoch_frames = self.M.get_frame_times_ms( block )
                for count, epoch in tqdm(enumerate(block['epoch']), desc="Processing epochs in block file {}".format(dataFile), total=len(block['epoch'])):
                    if count < len(epoch_frames):
                        parameters = epoch['parameters']
                        frame_times = epoch_frames[count]

                        pre_frames = round(parameters['preTime']*1e-3 * frame_rate)
                        pre_pts[total_epoch_count] = round(parameters['preTime']*1e-3 * bin_rate)
                        stim_pts[total_epoch_count] = round(parameters[self.stim_time_name]*1e-3 * bin_rate)
                        tail_pts[total_epoch_count] = round(parameters['tailTime']*1e-3 * bin_rate)
                        
                        # Check for frame drops.
                        frame_times, transition_frames = stimulus.check_frame_times(frame_times)
                        # Compute the time between frames.
                        diff_frames.append(np.diff(frame_times))

                        # Get the binned spike count.
                        binned_spikes = self.get_binned_spikes(frame_times, stride)
                        # print("t: " + str(time_bins) + "; sp:" + str(binned_spikes.shape[0]) + " " + str(binned_spikes.shape[1]))
                        for cell_count, cell in enumerate(self.cell_ids):
                            if cell in spike_dict:
                                if binned_spikes.shape[0] <= time_bins:
                                    spike_dict[cell][total_epoch_count,:binned_spikes.shape[0]] = binned_spikes[:,cell_count]
                                else:
                                    spike_dict[cell][total_epoch_count,:time_bins] = binned_spikes[:time_bins,cell_count]
                            else:
                                spike_dict[cell] = np.zeros((num_epochs, time_bins)).astype(float) #empty_array
                                if binned_spikes.shape[0] <= time_bins:
                                    spike_dict[cell][total_epoch_count,:binned_spikes.shape[0]] = binned_spikes[:,cell_count]
                                else:
                                    spike_dict[cell][total_epoch_count,:time_bins] = binned_spikes[:time_bins,cell_count]
                    # Iterate the count.
                    total_epoch_count += 1
        # Compute the mean frame rate.
        diff_frames = np.concatenate(diff_frames)
        diff_frames = diff_frames[diff_frames < 20.0]
        mean_frame_rate = 1000.0 / np.mean(diff_frames)
        # Sort the STA dictionary by cluster.
        spike_dict = dict(sorted(spike_dict.items()))

        cluster_id = np.array(list(spike_dict.keys()))

        return spike_dict, cluster_id, params, unique_params, pre_pts, stim_pts, tail_pts, mean_frame_rate
    
    def get_spike_times_and_parameters(self, 
                    protocolStr: str=None,
                    groupStr: str=None,
                    param_names: list=[],
                    sort_algorithm: str='kilosort2',
                    file_name: Union[str, list]=None,
                    bin_rate: float=1000.0,
                    sample_rate: float=20000.0):
        """
        Get the spike times (in msec) and desired epoch parameters
        Parameters:
            protocolStr: Substring within the protocol ID you want to search for. (Required)
            groupStr: Substring with the epoch group label to search for (optional).
            param_names: List of parameter names to search for.
            sort_algorith: Name of algorithm used for spike detection ('yass' or 'kilosort2')
            bin_rate: Bin rate for the spike count in Hz.
            sample_rate: Sample rate for the data in Hz (default=20000.0)
        """

        if protocolStr == None:
            print('Protocol ID is empty. Cannot continue...')
            return None

        if file_name is not None:
            protocol = self.M.search_data_file(protocolStr, file_name=file_name)
        else: 
            protocol = self.M.searchProtocol(protocolStr, groupStr)
            
        spike_dict = dict()

        # Take a quick pass through the data files.
        if not param_names:
            params = dict()
            unique_params = dict()
        else:
            params, unique_params = self.M.get_stimulus_parameters(protocol, param_names)

        # Get the spike count dictionary.
        spike_dict, cluster_id, pre_pts, stim_pts, tail_pts = self.get_spike_times_dict(protocol, 
                                        sort_algorithm=sort_algorithm,
                                        bin_rate=bin_rate,
                                        sample_rate=sample_rate)
        
        return spike_dict, np.array(cluster_id), params, unique_params, pre_pts, stim_pts, tail_pts

    def get_data_size(self, protocol: dict,
                            bin_rate: float=100.0):
        """
        Get the total number of epochs and the number of time bins.
        Parameters:
            protocol: Protocol dictionary.
            bin_rate: Bin rate for the spike count in Hz.

        Returns:
            total_epoch_count: Total number of epochs.
            time_bins: Number of time bins.
        """
        total_epoch_count = 0 # Keep track of the total number of epochs.
        time_bins = 0

        # Loop through the EpochGroups
        for group in protocol['group']:
            # Grab the epoch blocks
            for block in group['block']:
                # Loop through each epoch.
                for epoch in block['epoch']:

                    parameters = epoch['parameters']

                    # Iterate the count.
                    total_epoch_count += 1

                    # Get the total number of time bins at the bin rate.
                    time_pts = np.ceil((parameters['preTime']+parameters[self.stim_time_name]+parameters['tailTime'])*1e-3 * bin_rate).astype(int)
                    if time_pts > time_bins:
                        time_bins = time_pts
                    
        return total_epoch_count, time_bins

    def get_spike_times_dict(self, protocol: dict,
                            sort_algorithm: str='yass',
                            bin_rate: float=1000.0,
                            sample_rate: float=20000.0):
        """
        Get the spike times (in msec) in a dictionary.
        Parameters:
            protocolStr: Substring within the protocol ID you want to search for. (Required)
            groupStr: Substring with the epoch group label to search for (optional).
            param_names: List of parameter names to search for.
            sort_algorith: Name of algorithm used for spike detection ('yass' or 'kilosort2')
            bin_rate: Bin rate for the spike count in Hz.
            sample_rate: Sample rate for the data in Hz (default=20000.0)
        """
        
        total_epoch_count = 0 # Keep track of the total number of epochs.

        num_epochs, _ = self.get_data_size(protocol, bin_rate)
        pre_pts = np.zeros((num_epochs,))
        stim_pts = np.zeros((num_epochs,))
        tail_pts = np.zeros((num_epochs,))

        sp_dict = dict()

        # Get the number of groups.
        nGroups = len(protocol['group'])
        # Loop through the EpochGroups
        for groupCount, group in enumerate(protocol['group']):
            nBlocks = len(group['block'])
            # Grab the epoch blocks
            for blockCount, block in enumerate(group['block']):
                print('Processing group ' + str(groupCount+1) + ' of ' + str(nGroups) + ' and block ' + str(blockCount+1) + ' of ' + str(nBlocks))

                dataFile = self.M.get_data_file_from_block(block)

                # Try to load a data file.
                filepath = os.path.join(self.M.sortpath, self.experimentName, dataFile, sort_algorithm)

                self.load_vision_table(filepath, dataFile)

                epoch_frames = self.M.get_frame_times_ms( block )

                # Loop through each epoch.
                for count, epoch in tqdm(enumerate(block['epoch']), desc="Processing epochs in block file {}".format(dataFile), total=len(block['epoch'])):
                    if count < len(epoch_frames):
                        parameters = epoch['parameters']

                        pre_pts[total_epoch_count] = round(parameters['preTime']*1e-3 * bin_rate)
                        stim_pts[total_epoch_count] = round(parameters[self.stim_time_name]*1e-3 * bin_rate)
                        tail_pts[total_epoch_count] = round(parameters['tailTime']*1e-3 * bin_rate)

                        duration=(parameters['preTime']+parameters[self.stim_time_name]+parameters['tailTime'])*1e-3 * bin_rate

                        # Get the spike times.
                        for cell in self.cell_ids:
                            spike_times = self.vcd.get_spike_times_for_cell(cell) / sample_rate * bin_rate # ms
                                
                            # Take the spikes the occurred between the beginning and end of the epoch.
                            spike_times = spike_times - epoch_frames[count][0] / 1000 * bin_rate
                            spike_times = spike_times[np.where((spike_times >= 0) * (spike_times < duration))[0]]
                            if cell in sp_dict:
                                sp_dict[cell][total_epoch_count] = spike_times
                            else:
                                sp_dict[cell] = np.zeros((num_epochs,)).astype(object)
                                sp_dict[cell][total_epoch_count] = spike_times
                        # Iterate the count.
                        total_epoch_count += 1
        
        # Sort the spike dictionary by the cluster id.
        sp_dict = dict(sorted(sp_dict.items()))
        cluster_id = np.array(list(sp_dict.keys()))
        spike_dict = np.zeros((len(sp_dict),num_epochs), dtype=object)
        for i, key in enumerate(sp_dict.keys()):
            spike_dict[i,:] = sp_dict[key]

        return spike_dict, cluster_id, pre_pts, stim_pts, tail_pts

    def get_spike_count_dict(self, protocol: dict,
                            sort_algorithm: str='yass',
                            bin_rate: float=100.0,
                            sample_rate: float=20000.0):
        """
        Get the spike counts in a dictionary.
        Parameters:
            protocol: Substring within the protocol ID you want to search for. (Required)
            sort_algorith: Name of algorithm used for spike detection ('yass' or 'kilosort2')
            bin_rate: Bin rate for the spike count in Hz.
            sample_rate: Sample rate for the data in Hz (default=20000.0)
        """
        spike_dict = dict()
        
        total_epoch_count = 0 # Keep track of the total number of epochs.

        num_epochs, time_bins = self.get_data_size(protocol, bin_rate)
        # Create an empty numpy array to hold the data.
        pre_pts = np.zeros((num_epochs,))
        stim_pts = np.zeros((num_epochs,))
        tail_pts = np.zeros((num_epochs,))

        # Get the number of groups.
        nGroups = len(protocol['group'])
        # Loop through the EpochGroups
        for groupCount, group in enumerate(protocol['group']):
            nBlocks = len(group['block'])
            # Grab the epoch blocks
            for blockCount, block in enumerate(group['block']):
                print('Processing group ' + str(groupCount+1) + ' of ' + str(nGroups) + ' and block ' + str(blockCount+1) + ' of ' + str(nBlocks))

                dataFile = self.M.get_data_file_from_block(block)

                # Try to load a data file.
                filepath = os.path.join(self.M.sortpath, self.experimentName, dataFile, sort_algorithm)

                self.load_vision_table(filepath, dataFile)

                epoch_frames = self.M.get_frame_times_ms( block )

                # Loop through each epoch.
                for count, epoch in tqdm(enumerate(block['epoch']), desc="Processing epochs in block file {}".format(dataFile), total=len(block['epoch'])):
                    if count < len(epoch_frames):
                        parameters = epoch['parameters']

                        pre_pts[total_epoch_count] = round(parameters['preTime']*1e-3 * bin_rate)
                        stim_pts[total_epoch_count] = round(parameters[self.stim_time_name]*1e-3 * bin_rate)
                        tail_pts[total_epoch_count] = round(parameters['tailTime']*1e-3 * bin_rate)

                        # Get the binned spike count.
                        # if bin_rate == 60.0:
                        #     binned_spikes = self.get_binned_spikes(frame_times=epoch_frames, sample_rate=sample_rate)
                        # else:
                        binned_spikes = self.get_binned_count(epoch_frames[count][0], 
                            duration=parameters['preTime']+parameters[self.stim_time_name]+parameters['tailTime'], 
                            bin_rate=bin_rate,
                            sample_rate=sample_rate)

                        for cell_count, cell in enumerate(self.cell_ids):
                            if cell in spike_dict:
                                spike_dict[cell][total_epoch_count,:binned_spikes.shape[0]] = binned_spikes[:,cell_count]
                            else:
                                spike_dict[cell] = np.zeros((num_epochs, time_bins)).astype(float) #empty_array
                                spike_dict[cell][total_epoch_count,:binned_spikes.shape[0]] = binned_spikes[:,cell_count]
                        # Iterate the count.
                        total_epoch_count += 1

        # Sort the STA dictionary by cluster.
        spike_dict = dict(sorted(spike_dict.items()))

        cluster_id = list(spike_dict.keys())

        return spike_dict, cluster_id, pre_pts, stim_pts, tail_pts

    def get_interspike_intervals(self, bin_edges=np.linspace(0,300,301)):
        isi = dict()
        for cell in self.cell_ids:
            spike_times = self.vcd.get_spike_times_for_cell(cell) / 20000 * 1000 # ms
            
            # Compute the interspike interval
            if len(spike_times) > 1:
                isi_tmp = np.diff(spike_times)
                isi[cell] = np.histogram(isi_tmp,bins=bin_edges)[0]
            else:
                isi[cell] = np.zeros((len(bin_edges)-1,)).astype(int)
        
        return isi

    def merge_neurons_files(self, filepath: str, filenames: list, dataset_name: str):
        """
        Merge the neurons files from multiple data files.
        Parameters:
            filepath: Path to the data files (string).
            filenames: List of data file names (list)
            dataset_name: Name of the dataset (string)
        """
        spikes_by_cell_id_np = {}
        n_samples = 0
        for dataFile in filenames:
            fpath = os.path.join(filepath, dataFile, 'yass') 
            # fpath = os.path.join(filepath, dataFile, 'kilosort2') #os.path.join(self.M.sortpath, self.experimentName, dataFile, 'kilosort2')
            self.load_vision_table(fpath, dataFile)

            for cell in self.cell_ids:
                spike_times = self.vcd.get_spike_times_for_cell(cell) + n_samples
                if cell in spikes_by_cell_id_np:
                    np.append(spikes_by_cell_id_np[cell], np.array(spike_times))
                else:
                    spikes_by_cell_id_np[cell] = np.array(spike_times)
        
            if n_samples == 0:
                ttl_times = self.vcd.get_ttl_times()
            else: 
                np.append(ttl_times,self.vcd.get_ttl_times()+n_samples)

            n_samples += self.vcd.get_n_samples()

        with vw.NeuronsFileWriter(filepath, dataset_name) as nfw:
            nfw.write_neuron_file(spikes_by_cell_id_np, ttl_times, n_samples)

    def get_binned_count(self, start_sample, # Start sample in milliseconds
                duration: float, # Duration in milliseconds
                bin_rate: float=100.0, 
                sample_rate: float=20000.0):
        """
        Get the binned spike count for a given time window.
        Parameters:
            start_sample: Start sample in milliseconds (int).
            duration: Duration in milliseconds (float).
            bin_rate: Bin rate in Hz (float).
            sample_rate: Sample rate in Hz (float).

        Returns:
            binned_spikes: Binned spike count (numpy array).
        """
        # Bin width in ms.
        bin_width = 1000 / bin_rate

        # Number of time bins.
        num_bins = np.floor(duration / bin_width + 1).astype(int)

        # Bin times (in samples)
        time_bins = np.linspace(start_sample, start_sample+duration-1, num_bins).astype(int)

        # Loop through cells and bin spike train.
        binned_spikes = []

        for cell in self.cell_ids:
            spike_times = self.vcd.get_spike_times_for_cell(cell) / sample_rate * 1000 # ms
            binned_spikes.append(np.histogram(spike_times,bins=time_bins)[0])
        
        return np.asarray(binned_spikes).T
    
    
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

    def load_vision_table(self, filepath: str, dataFile: str):
        """
        Load the vision data table.
        Parameters:
            filepath: Path to the data files (string).
            dataFile: Data file name (string)
        """
        self.vcd = vl.load_vision_data(filepath, dataFile, include_neurons=True)
        # Get the cell Ids.
        self.cell_ids = sorted(self.vcd.get_cell_ids()) # sort by convention!

    def map_wn_to_ns_cellids(ns_vcd,ns_cellids,wn_vcd,
                            wn_cellids,celltypes_dict,
                            corr_dict = {'on parasol': .95,'off parasol': .95,
                                        'on midget': .95,'off midget': .95,
                                        'a1': .95},
                            mask=False,n_sigmas=None):
        """
        Maps WN to NS EIs according to a threshold value of correlation. Computes
        EI power over space, both with and without masking (user choice). Does a
        pass over the NS cellids and finds the corresponding WN cell. If none is
        found the cell doesn't get mapped (does not appear in the dictionary).
        Also writes a field for the normalized x,y locations for each RF.
        Parameters:
            ns_vcd: natural scenes vision data object
            ns_cellids: natural scenes cellids to map
            wn_cellids: white noise cellids to map
            celltypes_dict: dictionary mapping white noise cell ids to celltype.
        """

        channel_noise = wn_vcd.channel_noise

        # Initialize a dictionary and loop over the cells.
        cellids_dict = dict()

        for key in ['wn_to_ns','ns_to_wn','celltypes']:
            cellids_dict[key] = dict()

        for wn_cell in wn_cellids:

            # Get the cell type and write as well.
            celltype = celltypes_dict[wn_cell].lower()

            # Hardcode these for now TODO FIXME
            if "on" in celltype and "parasol" in celltype:
                celltype = 'on parasol'
            elif "off" in celltype and "parasol" in celltype:
                celltype = "off parasol"
            elif "on" in celltype and "midget" in celltype:
                celltype = 'on midget'
            elif "off" in celltype and "midget" in celltype:
                celltype = 'off midget'
            elif "a1" in celltype:
                celltype = 'a1'
            else:
                continue

            # If masking, only look at the significant indices. 
            wn_cell_ei = wn_vcd.get_ei_for_cell(wn_cell).ei

            if mask and n_sigmas is not None:
                sig_inds = np.argwhere(np.abs(np.amin(wn_cell_ei,axis=1))
                                    > n_sigmas * channel_noise).flatten()
                wn_cell_ei_power = np.zeros(wn_cell_ei.shape[0])
                wn_cell_ei_power[sig_inds] = np.sum(wn_cell_ei[sig_inds,:]**2,
                                                    axis=1)
            else:
                wn_cell_ei_power = np.sum(wn_cell_ei**2,axis=1)

            corrs = []

            for ns_cell in ns_cellids:
                ns_cell_ei = ns_vcd.get_ei_for_cell(ns_cell).ei

                if mask and n_sigmas is not None:
                    sig_inds = np.argwhere(np.abs(np.amin(ns_cell_ei,axis=1))
                                        > n_sigmas * channel_noise).flatten()
                    ns_cell_ei_power = np.zeros(ns_cell_ei.shape[0])
                    ns_cell_ei_power[sig_inds] = np.sum(ns_cell_ei[sig_inds,:]**2,
                                                    axis=1)
                else:
                    ns_cell_ei_power = np.sum(ns_cell_ei**2,axis=1)

                corr = np.corrcoef(wn_cell_ei_power,ns_cell_ei_power)[0,1]
                corrs.append(corr)

            # Take the cell with the largest correlation.
            if np.max(corrs) < corr_dict[celltype]:
                continue

            max_ind = np.argmax(np.asarray(corrs))
            cellids_dict['wn_to_ns'][wn_cell] = ns_cellids[max_ind]
            cellids_dict['ns_to_wn'][ns_cellids[max_ind]] = wn_cell
            cellids_dict['celltypes'][wn_cell] = celltype

            # Once the cell has been mapped, remove it (hack) FIXME.
            ns_cellids.remove(ns_cellids[max_ind]) 

class Stimulus(object):
    def __init__(self):
        """
        Constructor method.
        """
        self.stimulus = None

    def fix_frames(self, frames: np.ndarray, transition_frames: Union[np.ndarray, list]=None):
        """
        Find and fix dropped frames.

        Parameters:
            frames: 4D array of frames (n_frames, x, y, n_colors).
            transition_frames: 1D array of number of frames between transitions.
        
        Returns:
            new_frames: 4D array of frames with dropped frames fixed.
        """
        new_frames = frames
        frame_idx = 0
        for idx in range(len(transition_frames)):
            for idx2 in range(transition_frames[idx]):
                new_frames[frame_idx,:,:,:] = frames[idx,:,:,:]
        return new_frames
    
    def check_frame_times(self, frame_times: np.ndarray, frame_rate: float=60.0): 
        """
        Check the frame times for dropped frames.

        Parameters:
            frame_times: 1D array of frame times.
            frame_rate: frame rate of the stimulus.

        Returns:
            frame_times: 1D array of frame times with dropped frames fixed.
        """
        # Get the frame durations in milliseconds.
        frame_interval = 1000/frame_rate
        d_frames = np.diff(frame_times)
        # Get the number of frames between transitions/check for drops.
        transition_frames = np.round( frame_interval / d_frames ).astype(np.int32)
        # Check for frame drops.
        if np.amax(transition_frames) > 1:
            n_frames = np.sum(transition_frames)+1
            f_times = np.zeros((n_frames,), dtype=np.float64)
            frame_count = 0
            for idx in range(len(frame_times)-1):
                if transition_frames[idx] > 1:
                    this_frame = frame_times[idx]
                    next_frame = frame_times[idx+1]
                    new_times = np.linspace(this_frame, next_frame, transition_frames[idx], endpoint=False)
                    for new_t in new_times:
                        f_times[frame_count] = new_t
                        frame_count += 1
                else:
                    f_times[frame_count] = frame_times[idx]
                    frame_count += 1
                # Add in the last frame time.
                f_times[-1] = frame_times[-1]
            return f_times, transition_frames
        else: 
            return frame_times, transition_frames
    
    # Get the frame sequence for the FastNoiseStimulus
    def getFastNoiseFrames(self, numXStixels: int, 
                            numYStixels: int, 
                            numXChecks: int, 
                            numYChecks: int, 
                            chromaticClass: str, 
                            numFrames: int, 
                            stepsPerStixel: int, 
                            seed: int, 
                            frameDwell: int=1,
                            gaussianFilter: bool=False,
                            filterSdStixels: float=1.0):
        """
        Get the frame sequence for the FastNoiseStimulus.
        Parameters:
            numXStixels: number of stixels in the x direction.
            numYStixels: number of stixels in the y direction.
            numXChecks: number of checks in the x direction.
            numYChecks: number of checks in the y direction.
            chromaticClass: chromatic class of the stimulus.
            numFrames: number of frames in the stimulus.
            stepsPerStixel: number of steps per stixel.
            seed: seed for the random number generator.
            frameDwell: number of frames to dwell on each frame.

        Returns:
        frames: 4D array of frames (n_frames, x, y, n_colors).
        """
        # Seed the random number generator.
        np.random.seed( seed )

        # First, generate the larger grid of stixels.
        if (chromaticClass == 'BY'):
            tfactor = 2
        elif (chromaticClass == 'RGB'):
            tfactor = 3
        else: # Black/white checks
            tfactor = 1

        # Get the size of the time dimension; expands for RGB, etc.
        tsize = np.ceil(numFrames*tfactor/frameDwell).astype(int)
        
        if (tfactor == 2 and (tsize % 2) != 0):
            tsize += 1
        
        # Generate the random grid of stixels.
        gridValues = np.random.rand(tsize, numXStixels*numYStixels)
        gridValues = np.reshape(gridValues, (tsize, numXStixels, numYStixels))
        gridValues = np.transpose(gridValues, (0, 2, 1))
        gridValues = np.round(gridValues)
        gridValues = (2*gridValues-1).astype(np.float32) # Convert to contrast

        # Filter the stixels if indicated.
        if gaussianFilter:
            for i in range(tsize):
                # frame_tmp = gaussian_filter(gridValues[i,:,:], sigma=filterSdStixels, order=0, mode='wrap', radius=np.ceil(2*filterSdStixels).astype(int))
                frame_tmp = gaussian_filter(gridValues[i,:,:], sigma=filterSdStixels, order=0, mode='wrap')
                gridValues[i,:,:] = 0.5*frame_tmp/np.std(frame_tmp)
            gridValues[gridValues > 1.0] = 1.0
            gridValues[gridValues < -1.0] = -1.0

        # Translate to the full grid
        fullGrid = np.zeros((tsize,numYStixels*stepsPerStixel,numXStixels*stepsPerStixel), dtype=np.float32)
        # fullGrid = np.zeros((tsize,numYStixels*stepsPerStixel,numXStixels*stepsPerStixel), dtype=np.uint8)

        for k in range(numYStixels*stepsPerStixel):
            yindex = math.floor(k/stepsPerStixel)
            for m in range(numXStixels*stepsPerStixel):
                xindex = math.floor(m/stepsPerStixel)
                fullGrid[:, k, m] = gridValues[:, yindex, xindex]

        # Generate the motion trajectory of the larger stixels.
        np.random.seed( seed ) # Re-seed the number generator

        # steps = np.round( (stepsPerStixel-1) * np.random.rand(tsize, 2) )
        steps = np.round( (stepsPerStixel-1) * np.random.rand(tsize, 2) )
        steps[:,0] = (stepsPerStixel-1) - steps[:,0]
        # steps = (stepsPerStixel-1) - np.round( (stepsPerStixel-1) * np.random.rand(tsize, 2) )
        # Get the frame values for the finer grid.
        # frameValues = np.zeros((tsize,numYChecks,numXChecks),dtype=np.uint8)
        frameValues = np.zeros((tsize,numYChecks,numXChecks),dtype=np.float32)
        for k in range(tsize):
            x_offset = steps[math.floor(k/tfactor), 0].astype(int)
            y_offset = steps[math.floor(k/tfactor), 1].astype(int)
            frameValues[k,:,:] = fullGrid[k, y_offset : numYChecks+y_offset, x_offset : numXChecks+x_offset]

        # Create your output stimulus. (t, y, x, color)
        stimulus = np.zeros((np.ceil(tsize/tfactor).astype(int),numYChecks,numXChecks,3), dtype=np.float32)

        # Get the pixel values into the proper color channels
        if (chromaticClass == 'BY'):
            stimulus[:,:,:,0] = frameValues[0::2,:,:]
            stimulus[:,:,:,1] = frameValues[0::2,:,:]
            stimulus[:,:,:,2] = frameValues[1::2,:,:]
        elif (chromaticClass == 'RGB'):
            stimulus[:,:,:,0] = frameValues[0::3,:,:]
            stimulus[:,:,:,1] = frameValues[1::3,:,:]
            stimulus[:,:,:,2] = frameValues[2::3,:,:]
        else: # Black/white checks
            stimulus[:,:,:,0] = frameValues
            stimulus[:,:,:,1] = frameValues
            stimulus[:,:,:,2] = frameValues
        # return stimulus

        # Deal with the frame dwell.
        if frameDwell > 1:
            stim = np.zeros((numFrames,numYChecks,numXChecks,3), dtype=np.float32)
            for k in range(numFrames):
                idx = np.floor(k / frameDwell).astype(int)
                stim[k,:,:,:] = stimulus[idx,:,:,:]
            return stim
        else:
            return stimulus

    def get_spatial_noise_frames(self, numXStixels: int, 
                            numYStixels: int, 
                            numXChecks: int, 
                            numYChecks: int, 
                            chromaticClass: str, 
                            unique_frames: int, 
                            repeat_frames: int,
                            stepsPerStixel: int, 
                            seed: int, 
                            frameDwell: int=1,
                            gaussianFilter: bool=False,
                            filterSdStixels: float=1.0):
        """
        Get the frame sequence for the FastNoiseStimulus.
        Parameters:
            numXStixels: number of stixels in the x direction.
            numYStixels: number of stixels in the y direction.
            numXChecks: number of checks in the x direction.
            numYChecks: number of checks in the y direction.
            chromaticClass: chromatic class of the stimulus.
            numFrames: number of frames in the stimulus.
            stepsPerStixel: number of steps per stixel.
            seed: seed for the random number generator.
            frameDwell: number of frames to dwell on each frame.

        Returns:
        frames: 4D array of frames (n_frames, x, y, n_colors).
        """
        # Seed the random number generator.
        np.random.seed( seed )

        # First, generate the larger grid of stixels.
        if (chromaticClass == 'BY'):
            tfactor = 2
        elif (chromaticClass == 'RGB'):
            tfactor = 3
        else: # Black/white checks
            tfactor = 1

        # Get the size of the time dimension; expands for RGB, etc.
        numFrames = unique_frames + repeat_frames
        tsize = np.ceil(numFrames*tfactor/frameDwell).astype(int)
        usize = np.ceil(unique_frames*tfactor/frameDwell).astype(int)
        rsize = np.ceil(repeat_frames*tfactor/frameDwell).astype(int)
        
        if (tfactor == 2 and (tsize % 2) != 0):
            tsize += 1
            rsize += 1
        
        # Generate the random grid of stixels.
        gridValues = np.zeros((tsize, numXStixels*numYStixels), dtype=np.float32)
        gridValues[:usize,:] = np.random.rand(usize, numXStixels*numYStixels)
        # Repeating sequence.
        if repeat_frames > 0:
            # Reseed the generator.
            np.random.seed( 1 )
            gridValues[usize:,:] = np.random.rand(rsize, numXStixels*numYStixels)
        gridValues = np.reshape(gridValues, (tsize, numXStixels, numYStixels))
        gridValues = np.transpose(gridValues, (0, 2, 1))
        gridValues = np.round(gridValues)
        gridValues = (2*gridValues-1).astype(np.float32) # Convert to contrast

        # Filter the stixels if indicated.
        if gaussianFilter:
            for i in range(tsize):
                # frame_tmp = gaussian_filter(gridValues[i,:,:], sigma=filterSdStixels, order=0, mode='wrap', radius=np.ceil(2*filterSdStixels).astype(int))
                frame_tmp = gaussian_filter(gridValues[i,:,:], sigma=filterSdStixels, order=0, mode='wrap')
                gridValues[i,:,:] = 0.5*frame_tmp/np.std(frame_tmp)
            gridValues[gridValues > 1.0] = 1.0
            gridValues[gridValues < -1.0] = -1.0

        # Translate to the full grid
        fullGrid = np.zeros((tsize,numYStixels*stepsPerStixel,numXStixels*stepsPerStixel), dtype=np.float32)
        # fullGrid = np.zeros((tsize,numYStixels*stepsPerStixel,numXStixels*stepsPerStixel), dtype=np.uint8)

        for k in range(numYStixels*stepsPerStixel):
            yindex = math.floor(k/stepsPerStixel)
            for m in range(numXStixels*stepsPerStixel):
                xindex = math.floor(m/stepsPerStixel)
                fullGrid[:, k, m] = gridValues[:, yindex, xindex]

        # Generate the motion trajectory of the larger stixels.
        np.random.seed( seed ) # Re-seed the number generator

        # steps = np.round( (stepsPerStixel-1) * np.random.rand(tsize, 2) )
        steps = np.round( (stepsPerStixel-1) * np.random.rand(tsize, 2) )
        steps[:,0] = (stepsPerStixel-1) - steps[:,0]
        # steps = (stepsPerStixel-1) - np.round( (stepsPerStixel-1) * np.random.rand(tsize, 2) )
        # Get the frame values for the finer grid.
        # frameValues = np.zeros((tsize,numYChecks,numXChecks),dtype=np.uint8)
        frameValues = np.zeros((tsize,numYChecks,numXChecks),dtype=np.float32)
        for k in range(tsize):
            x_offset = steps[math.floor(k/tfactor), 0].astype(int)
            y_offset = steps[math.floor(k/tfactor), 1].astype(int)
            frameValues[k,:,:] = fullGrid[k, y_offset : numYChecks+y_offset, x_offset : numXChecks+x_offset]

        # Create your output stimulus. (t, y, x, color)
        stimulus = np.zeros((np.ceil(tsize/tfactor).astype(int),numYChecks,numXChecks,3), dtype=np.float32)

        # Get the pixel values into the proper color channels
        if (chromaticClass == 'BY'):
            stimulus[:,:,:,0] = frameValues[0::2,:,:]
            stimulus[:,:,:,1] = frameValues[0::2,:,:]
            stimulus[:,:,:,2] = frameValues[1::2,:,:]
        elif (chromaticClass == 'RGB'):
            stimulus[:,:,:,0] = frameValues[0::3,:,:]
            stimulus[:,:,:,1] = frameValues[1::3,:,:]
            stimulus[:,:,:,2] = frameValues[2::3,:,:]
        else: # Black/white checks
            stimulus[:,:,:,0] = frameValues
            stimulus[:,:,:,1] = frameValues
            stimulus[:,:,:,2] = frameValues
        # return stimulus

        # Deal with the frame dwell.
        if frameDwell > 1:
            stim = np.zeros((numFrames,numYChecks,numXChecks,3), dtype=np.float32)
            for k in range(numFrames):
                idx = np.floor(k / frameDwell).astype(int)
                stim[k,:,:,:] = stimulus[idx,:,:,:]
            return stim
        else:
            return stimulus

    def get_sparse_noise_frames(self, 
                            numXStixels: int, 
                            numYStixels: int, 
                            numXChecks: int, 
                            numYChecks: int, 
                            chromaticClass: str, 
                            numFrames: int, 
                            stepsPerStixel: int, 
                            seed: int, 
                            frameDwell: int=18,
                            pixelDensity: float=0.01):
        """
        Generate sparse noise frames.

        Parameters:
            numXStixels: number of stixels in the x direction.
            numYStixels: number of stixels in the y direction.
            numXChecks: number of checks in the x direction.
            numYChecks: number of checks in the y direction.
            chromaticClass: 'BY', 'RGB', or 'BW'.
            numFrames: number of frames to generate.
            stepsPerStixel: number of steps per stixel.
            seed: random number generator seed.
            frameDwell: number of frames to dwell on each frame.
            pixelDensity: fraction of pixels to be active.
        
        Returns:
            stimulus: a numpy array of size (numFrames, numYChecks, numXChecks, 3)
        """
        # Seed the random number generator.
        np.random.seed( seed )

        # First, generate the larger grid of stixels.
        if (chromaticClass == 'BY'):
            tfactor = 2
        elif (chromaticClass == 'RGB'):
            tfactor = 3
        else: # Black/white checks
            tfactor = 1

        # Get the size of the time dimension; expands for RGB, etc.
        tsize = np.ceil(numFrames*tfactor/frameDwell).astype(int)
        
        if (tfactor == 2 and (tsize % 2) != 0):
            tsize += 1

        gridValues = np.random.rand(tsize, numXStixels*numYStixels)
        gridValues = np.reshape(gridValues, (tsize, numXStixels, numYStixels))
        gridValues = np.transpose(gridValues, (0, 2, 1))

        # Get the active pixels based on the pixelDensity.
        if (chromaticClass == 'BY'):
            gridValues[(gridValues < 1-pixelDensity/2) & (gridValues > pixelDensity/2)] = 0.5
            gridValues[gridValues > 0.5] = 1
            gridValues[gridValues < 0.5] = 0
        elif (chromaticClass == 'RGB'):
            gridValues[(gridValues < 1-pixelDensity/2) & (gridValues > pixelDensity/2)] = 0.5
            gridValues[gridValues > 0.5] = 1
            gridValues[gridValues < 0.5] = 0
        else: # Black/white checks
            gridValues[(gridValues < 1-pixelDensity/2) & (gridValues > pixelDensity/2)] = 0.5
            gridValues[gridValues > 0.5] = 1
            gridValues[gridValues < 0.5] = 0

        gridValues = (2*gridValues-1).astype(np.float32) # Convert to contrast

        # Translate to the full grid
        fullGrid = np.zeros((tsize,numYStixels*stepsPerStixel,numXStixels*stepsPerStixel), dtype=np.float32)

        for k in range(numYStixels*stepsPerStixel):
            yindex = math.floor(k/stepsPerStixel)
            for m in range(numXStixels*stepsPerStixel):
                xindex = math.floor(m/stepsPerStixel)
                fullGrid[:, k, m] = gridValues[:, yindex, xindex]

        # Generate the motion trajectory of the larger stixels.
        np.random.seed( seed ) # Re-seed the number generator

        steps = np.round( (stepsPerStixel-1) * np.random.rand(tsize, 2) )
        steps[:,0] = (stepsPerStixel-1) - steps[:,0]

        # steps = (stepsPerStixel-1) - np.round( (stepsPerStixel-1) * np.random.rand(tsize, 2) )
        # Get the frame values for the finer grid.
        frameValues = np.zeros((tsize,numYChecks,numXChecks),dtype=np.float32)
        for k in range(tsize):
            x_offset = steps[math.floor(k/tfactor), 0].astype(int)
            y_offset = steps[math.floor(k/tfactor), 1].astype(int)
            frameValues[k,:,:] = fullGrid[k, y_offset : numYChecks+y_offset, x_offset : numXChecks+x_offset]

        # Create your output stimulus. (t, y, x, color)
        stimulus = np.zeros((np.ceil(tsize/tfactor).astype(int),numYChecks,numXChecks,3), dtype=np.float32)

        # Get the pixel values into the proper color channels
        if (chromaticClass == 'BY'):
            stimulus[:,:,:,0] = frameValues[0::2,:,:]
            stimulus[:,:,:,1] = frameValues[0::2,:,:]
            stimulus[:,:,:,2] = frameValues[1::2,:,:]
        elif (chromaticClass == 'RGB'):
            stimulus[:,:,:,0] = frameValues[0::3,:,:]
            stimulus[:,:,:,1] = frameValues[1::3,:,:]
            stimulus[:,:,:,2] = frameValues[2::3,:,:]
        else: # Black/white checks
            stimulus[:,:,:,0] = frameValues
            stimulus[:,:,:,1] = frameValues

        # Deal with the frame dwell.
        if frameDwell > 1:
            stim = np.zeros((numFrames,numYChecks,numXChecks,3), dtype=np.float32)
            for k in range(numFrames):
                idx = np.floor(k / frameDwell).astype(int)
                stim[k,:,:,:] = stimulus[idx,:,:,:]
            return stim
        else:
            return stimulus

    # Get the frame sequence for the PinkNoise Stimulus
    def getPinkNoiseFrames(self, 
                            numXChecks: int, 
                            numYChecks: int, 
                            numFrames: int, 
                            noiseContrast: float=0.35, 
                            spatialAmplitude: float=1.0,
                            temporalAmplitude: float=0.1,
                            chromaticClass: str='BY',
                            seed: int=1):

        # Seed the random number generator.
        np.random.seed( seed )

        # Generate the spatial frequencies.
        x = np.concatenate((np.arange(np.floor(float(numXChecks)/2.0)+1.0),np.arange(-np.ceil(float(numXChecks)/2.0)+1.0,0.0))) / float(numXChecks)
        x = np.abs(x)
        # Reproduce the frequencies along every row.
        x = np.tile(np.reshape(x,(numXChecks,1)),(1,numYChecks))

        # Y-axis frequencies.
        y = np.concatenate((np.arange(np.floor(float(numYChecks)/2.0)+1.0), np.arange(-np.ceil(float(numYChecks)/2.0)+1.0,0.0))) / float(numXChecks)
        y = np.abs(y)
        y = np.tile(np.reshape(y,(1,numYChecks)),(numXChecks,1))

        # Spatial frequency amplitude spectrum.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            sf = (x**2 + y**2) ** -(spatialAmplitude/2.0)
            sf[sf == np.inf] = 0.0
            sf *= 0.5
            sf = sf.T

        # Temporal frequencies.
        if temporalAmplitude > 0.0:
            t = np.concatenate((np.arange(np.floor(float(numFrames)/2.0)+1.0), np.arange(-np.ceil(float(numFrames)/2.0)+1.0,0.0))) / float(numFrames)
            t = np.abs(t)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                t_freq = t ** -temporalAmplitude
                t_freq[t_freq == np.inf] = 0.0
                t_freq /= np.max(t_freq)
                t_freq[0] = 1.0

            # Generate the spatiotemporal spectrum.
            st_freq = np.outer(sf.ravel(), t_freq)
            st_freq = np.reshape(st_freq,(numYChecks,numXChecks,numFrames,1))

        # Generate the random phases.
        if chromaticClass == 'RGB':
            phi = np.random.rand(3,numFrames,numXChecks,numYChecks)
            st_freq = np.tile(st_freq,(1,1,1,3))
        elif chromaticClass == 'BY':
            phi = np.random.rand(2,numFrames,numXChecks,numYChecks)
            st_freq = np.tile(st_freq,(1,1,1,2))
        else: # Grayscale
            phi = np.random.rand(1,numFrames,numXChecks,numYChecks)

        ## insertions
        phi *= 2*np.pi

        prod = np.cos(phi) + 1j*np.sin(phi)
        prod *= np.transpose(st_freq,(3,2,1,0))
        ft = np.fft.ifftn(prod).real

        if chromaticClass == 'RGB':
            stimulus = ft
        elif chromaticClass == 'BY':
            stimulus = np.zeros((3,numFrames,numXChecks,numYChecks))
            stimulus[0,:,:,:] = ft[0,:,:,:].copy()
            stimulus[1,:,:,:] = ft[0,:,:,:].copy()
            stimulus[2,:,:,:] = ft[1,:,:,:].copy()
        else: # Grayscale
            stimulus = np.tile( ft, (3,1,1,1) )

        # Normalize the stimulus.
        stimulus /= np.std(stimulus)
        stimulus *= noiseContrast
        stimulus[stimulus > 1.0] = 1.0
        stimulus[stimulus < -1.0] = -1.0

        # Create your output stimulus. (t, y, x, color)
        stimulus = np.transpose(stimulus, axes=(1,3,2,0))

        return stimulus
    
    def upsample_frames(self, frames, upsample_factor: int):
        nt = frames.shape[0]*upsample_factor-1
        stim = np.zeros((nt,frames.shape[1],frames.shape[2],frames.shape[3]), dtype=np.float32)
        for k in range(nt):
            idx = np.floor(k / upsample_factor).astype(int)
            stim[k,:,:,:] = frames[idx,:,:,:]
        return stim


class Analysis(object):
    def __init__(self):
        self.chicken = None

    def get_stimulus_parameters(self, protocol, param_names):
        orientations = list()
        contrasts = list()
        for group in protocol['group']:
            for block in group['block']:
                dataFile = self.M.get_data_file_from_block(block)

                # Try to load a data file.
                filepath = os.path.join(self.M.sortpath, self.experimentName, dataFile, self.sort_algorithm)

                self.D.load_vision_table(filepath, dataFile)

                for count, epoch in enumerate(block['epoch']):

                    parameters = epoch['parameters']

                    # Get the orientation and contrast.
                    orientations.append(parameters['orientation'])
                    contrasts.append(parameters['contrast'])
        return orientations, contrasts
    
    def compute_DSI(self, orientations, # Orientation of stimuli in degrees.
                y):

        # If the y-values are all zero, don't bother running.
        if np.all(y == 0):
            return 0.0, 0.0
        
        # Get the orientation angles in radians.
        angles = orientations.astype('float') / 180.0 * np.pi
        w = y / np.max(y)

        # Get the vector sum of responses across the angles.
        r = np.sum(w * np.exp(1j*angles)) / np.sum(w)

        # DSI is the absolute value of the vector sum.
        dsi_global = np.abs(r)

        # The preferred angle is the complex phase of the vector sum.
        pref_angle = np.angle(r) / np.pi * 180 # Convert to degrees.
        return dsi_global, pref_angle

    # Orientation Selectivity Index
    def compute_OSI(self, orientations, # Orientation of stimuli in degrees.
                y):

        # If the y-values are all zero, don't bother running.
        if np.all(y == 0):
            return 0.0, 0.0
        
        # Get the orientation angles in radians.
        angles = orientations.astype('float') / 180.0 * np.pi
        w = y / np.max(y)

        # Get the vector sum of responses across the angles.
        r = np.sum(w * np.exp(2j*angles)) / np.sum(w)

        # OSI is the absolute value of the vector sum.
        osi_global = np.abs(r)

        # The preferred angle is the half of the complex phase of the vector sum.
        pref_angle = np.angle(r) / 2.0 / np.pi * 180 # Convert to degrees.
        return osi_global, pref_angle

    # Save data to a MATLAB compatible file.
    def save_mat(self, filepath, mdic):
        # Decide which format to use based on the size of the data dictionary.
        if getsizeof(mdic) < 4e9:
            savemat(filepath, mdic)
        else:
            hdf5storage.savemat( filepath, mdic, format=7.3, matlab_compatible=True, compress=False )
        return None

    def load_mat(self, filepath):
        mdic = hdf5storage.loadmat(filepath)
        return mdic

    def save_pickle(self, filepath, d):
        pickle.dump(d, open(filepath, 'wb'))
        return None

    def load_sta_pickle(self, filepath):
        d = pickle.load( open(filepath, 'rb') )
        return d


