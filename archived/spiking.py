"""
Module of functions for analysis of spiking properties: EIs, ACFs, etc.
"""
import numpy as np
import scipy as sp
from scipy import linalg
import src.config as cfg
from scipy import stats, signal
from itertools import combinations
import electrode_map as elcmp
import pdb

def compute_acf(vcd,cell,delays=cfg.ACF_DELAYS):
    """
    Computes the autocorrelation function of the spike train with lags given
    in delays. Omits the zero bin. L2-normalizes the histogram.
    Parameters:
        vcd: vision cell data table object
        cell: cell id
        delays: delay (in ms)
    Returns:
        tuple of raw and L2 normalized spike train autocorrelation function.
    """
    def bin_spikes():
        """
        Helper function for binning spike trains with 1 ms precision. to be 
        used ONLY for correlations, and not for model fitting.

        Returns:
            binned spike train with 1 ms precision.
        """
        n_sec = vcd.n_samples / cfg.FS
        bin_edges = np.linspace(0, n_sec* 1000,
                         (int(n_sec/cfg.ACF_BIN_SIZE) * 1000)+1)
        spiketimes = (vcd.get_spike_times_for_cell(cell) / cfg.FS) * 1000

        return np.histogram(spiketimes,bin_edges)[0]
        
    # Get the binned spikes.
    binned_spikes = bin_spikes()

    # Find nonzero indices (spike times) and initialize histogram
    spiketimes1 = np.argwhere(binned_spikes).flatten()
    spiketimes2 = np.argwhere(binned_spikes).flatten()
    acf_out = np.zeros(delays.shape)

    # Loop through spike train and get frequency of each delay.
    for spiketime in spiketimes1:
        diff = spiketimes2 - spiketime
        inds = np.where((diff >= -cfg.ACF_DELAY) & (diff <= cfg.ACF_DELAY))[0]
        
        for ind in inds:
            
            if np.argwhere(delays == diff[ind]).flatten().shape[0] == 0:
                 continue
               
            delay_ind = np.argwhere(delays == diff[ind]).flatten()[0]
            acf_out[delay_ind] += 1

    return acf_out,acf_out / linalg.norm(acf_out)

'''
The following  functions are for computing axon conduction velocity from the 
EI. The original code was written by Moosa Zaidi, with some refractoring and
cleanup done by me.
'''
def get_axonal_conduction_velocity(ei, upsample_factor = 10, threshold_factor=60, 
min_num_axonal_elecs = 8, num_pairs=10):
    """ 
    Estimates the axonal conduciton velocity of the input EI 
  
    see detailed decription in "Distinguishing midget and parasol cells using axon conduction velocity" 
    portion of the Methods sections of https://doi.org/10.1101/2020.08.23.263608
    Parameters: 
    ei (num_electrodes x num_time_points numpy array): the EI whose conduction velocity needs to be computed
    upsample_factor (int): factor by which to upsample EI in time.
    threshold_factor (int): number of robust standard deviations of voltage to use as threshold for electrode significance
    min_num_axonal_elecs (int): minimum number of axonal electrodes required to return condcution velocity estimate.
                                if number of axonal electrodes less than this number, return numpy.nan
    num_pairs (int): number of axonal electrode pairs to use in final estimate (e.g. top 10 weighted electrode pairs)
  
    Returns: 
    numpy.float64: estimated axonal conduction velocity of input ei. Returns numpy.nan if number of axonal electrodes
                   is less than min_num_axonal_elecs
  
    """

    # Based on the array, get the coords and the adj mat.
    if ei.shape[0] == 519:
        array_id = 1502
    elif ei.shape[0] == 512:
        array_id = 502
    else:
        raise Exception("Input EI of unknown shape.")

    coords = elcmp.get_litke_array_coordinates_by_array_id(array_id)
    adj_mat = get_adj_mat(coords)

    #upsample ei in time
    ei_upsampled = upsample_ei(ei, upsample_factor)
    
    #compute electrode significance threshold  as a multiple of the robust standard deviation of voltage values  
    #(robust std of ei voltages represents noise level)

    # Make sure it's robust to the zero channels on 519.
    non_zero_inds = np.argwhere(np.amin(ei,axis=1) != 0).flatten()
    std = stats.median_absolute_deviation(ei[non_zero_inds,:],axis=None)
    threshold = threshold_factor*std
    # determine the type of each electrode ("low_amp", "soma", "axon", "dendrite", or "mixed")
    electrode_types = filter_ei_for_electrode_types(ei_upsampled,thr=threshold)
   
    axonal_elecinds = get_axonal_elec_ind(electrode_types, adj_mat)
    if len(axonal_elecinds) < min_num_axonal_elecs:
        return np.nan

    min_ts = np.argmin(ei_upsampled,axis=1) #time of minimum volage per electrode (in time index of upsampled ei)
    #Across all pairs of axonal electrodes, compute for each pair: the conduction velocity estimated for that pair and a
    #weight for that pair.  The weight of a pair of electrodes is the product of their maximum voltage amplitudes. Exclude 
    #pairs whose difference in time at minimum are less than or equal to the original sampling period.
    velocities = []
    weights = []
    min_ts = np.argmin(ei_upsampled,axis=1) #time of minimum volage per electrode (in time index of upsampled ei)

    # Get the largest electrode ind.
    max_elec_ind = np.argmin(np.amin(ei,axis=1))

    for e1, e2 in combinations(axonal_elecinds, 2):
        dt = 1.0 * np.abs(min_ts[e1]-min_ts[e2]) / upsample_factor
        if dt <= 1:
            continue

        dist = np.linalg.norm(coords[e1]-coords[e2])
        velocities.append(dist/dt)
        weights.append(np.max(np.abs(ei_upsampled[e1]))*np.max(np.abs(ei_upsampled[e2])))

    #return weighted average of conduction velocity estimates using the top num_pairs axonal electrode pairs
    tups = []
    for idx in range(len(velocities)):
        tups.append((weights[idx],velocities[idx]))
    tups = sorted(tups)
    tups= tups[-num_pairs:]
    tups = np.array(tups)

    if tups.shape[0] == 0:
        return np.nan

    return np.average(tups[:,1],weights=tups[:,0])

def get_adj_mat(coords):
    """ 
    Return adjacency matrix for 512/519 electrode array
  
    Returns: 
    dict: key: index of query electrodes, value: array of indices of all electrodes that are within 70 um
                                                 of the query electrode (including the query electrode)
    """

    # Figure out spacing.
    if coords.shape[0] == 519:
        spacing = 35
    elif coords.shape[0] == 512:
        spacing = 70 
    else:
        raise Exception("Coordindates of unexpected shape.")

    adj_mat = {}

    for idx in range(coords.shape[0]):
        adj_mat[idx]= np.where(np.linalg.norm(coords-coords[idx],axis=1)<=spacing)[0]

    return adj_mat

def upsample_ei(ei, upsample_factor):
    """ 
    upsample ei in time according by factor of upsample_factor using Fourier interpolation
    """
    (num_electrodes, num_timepoints) = np.shape(ei)
    ei_upsampled = np.zeros((num_electrodes, num_timepoints*upsample_factor))
    for idx in range(num_electrodes):
        ei_upsampled[idx,:] = signal.resample(ei[idx,:], num_timepoints*upsample_factor)
    return ei_upsampled

# determine the type of each electrode ("low_amp", "soma", "axon", "dendrite", or "mixed")
#adapted from Eric functions in eilib
def filter_ei_for_electrode_types(ei, thr = 5):
    electrode_types = []
    for ei_per_elec in ei:
        if (np.min(ei_per_elec) > -np.abs(thr)):
            electrode_types.append("low_amp")
        else:
            #electrode_types.append(el.axonorsomaRatio(ei_per_elec))
            electrode_types.append(axon_or_soma_ratio(ei_per_elec))
    return np.array(electrode_types)

def get_axonal_elec_ind(electrode_types, adj_mat):
    """ 
    Returns the indices of the axonal electrodes. Any neighbors of soma electrodes are excluded
    Parameters: 
    electrode_types (length num_electrodes array of strings): each element is the type of that electrode
  
    Returns: 
    array of ints: indices of axonal electrodes
  
    """
    axonal_elecinds = np.where(electrode_types=='axon')[0]
    axonal_elecinds = set(axonal_elecinds)

    #exclude soma electrodes and their neighbors
    soma_elecinds = np.where(electrode_types=='soma')[0]
    soma_neighbors = set(soma_elecinds)

    for elec in soma_elecinds:
        soma_neighbors = soma_neighbors.union(set(adj_mat[elec]))
    axonal_elecinds= axonal_elecinds.difference(soma_neighbors)

    # Enforce that the axonal electrodes are only neighbors with axons.
    axonal_elecinds_pruned = []

    for axonal_elecind in axonal_elecinds:
        neighbor_types = electrode_types[adj_mat[axonal_elecind]]

        if "soma" in neighbor_types or "dendrite" in neighbor_types or\
            "mixed" in neighbor_types:
                continue

        axonal_elecinds_pruned.append(axonal_elecind)

    #return axonal_elecinds
    return axonal_elecinds_pruned

def axon_or_soma_ratio(wave,uppBound=1.6,lowBound=0.05):
    """
    This function was written by Sasi Madugula and is part of the 
    heuristic-based algorithm for determining electrode-compartment.
    """
    try:
        #Get index and value of (negative) min {only one}
        minind = np.argmin(wave)
        minval = np.min(wave)

        #Get max vals to either side of min
        maxvalLeft = np.max(wave[0:minind])
        maxvalRight = np.max(wave[minind:])

        if np.abs(minval) < max(maxvalLeft,maxvalRight):
            rCompt = 'dendrite'
        else:
            if maxvalRight == 0:
                ratio = 0
            else:
                ratio = maxvalLeft/maxvalRight
            if ratio > uppBound:
                rCompt = 'axon'
            elif ratio < lowBound: #FUDGED
                rCompt = 'soma'
            else:
                rCompt = 'mixed'
    except ValueError:
        rCompt = 'error' #wave is abnormally shaped (usually has min at leftmost or rightmost point)

    return rCompt