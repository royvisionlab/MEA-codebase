import numpy as np
import os
import torch
from kilosort import run_kilosort, io

import time
import kilosort
from kilosort import (
    preprocessing,
    datashift,
    template_matching,
    io,
    spikedetect,
    clustering_qr,
    CCG,
    PROBE_DIR
)
from kilosort.parameters import DEFAULT_SETTINGS

# import clustering_qr

from kilosort.run_kilosort import set_files, initialize_ops, compute_preprocessing, compute_drift_correction, save_sorting, detect_spikes, cluster_spikes


def detect_spikes(ops, device, bfile, tic0=np.nan, progress_bar=None):
    """Run spike sorting algorithm and save intermediate results to `ops`.
    
    Parameters
    ----------
    ops : dict
        Dictionary storing settings and results for all algorithmic steps.
    device : torch.device
        Indicates whether `pytorch` operations should be run on cpu or gpu.
    bfile : kilosort.io.BinaryFiltered
        Wrapped file object for handling data.
    tic0 : float; default=np.nan.
        Start time of `run_kilosort`.
    progress_bar : TODO; optional.
        Informs `tqdm` package how to report progress, type unclear.

    Returns
    -------
    st : np.ndarray
        1D vector of spike times for all clusters.
    clu : np.ndarray
        1D vector of cluster ids indicating which spike came from which cluster,
        same shape as `st`.
    tF : np.ndarray
        TODO
    Wall : np.ndarray
        TODO

    """
    tic = time.time()
    st0, tF, ops = spikedetect.run(ops, bfile, device=device, progress_bar=progress_bar)
    tF = torch.from_numpy(tF)
    if len(st0) == 0:
        raise ValueError('No spikes detected, cannot continue sorting.')

    tic = time.time()
    clu, Wall = clustering_qr.run(ops, st0, tF, mode='spikes', device=device,
                                  progress_bar=progress_bar)
    Wall3 = template_matching.postprocess_templates(Wall, ops, clu, st0, device=device) 
    tic = time.time()
    st, tF, ops = template_matching.extract(ops, bfile, Wall3, device=device,
                                                 progress_bar=progress_bar)

    return st, tF, Wall, clu

def cluster_spikes(st, tF, ops, device, bfile, tic0=np.nan, progress_bar=None):
    clu, Wall = clustering_qr.run(ops, st, tF,  mode = 'template', device=device,
                                  progress_bar=progress_bar)
    Wall, clu, is_ref = template_matching.merging_function(ops, Wall, clu, st[:,0],
                                                           device=device)
    clu = clu.astype('int32')
    bfile.close()
    return clu, Wall

def get_settings(
        num_channels: int=512, 
        electrode_spacing: float=60.0, 
        use_drift_correction: bool=False,
        batch_size: int=40000):
    if use_drift_correction:
        nblocks = 1
    else:
        nblocks = 0
    settings = {'n_chan_bin': num_channels,
                'fs': 20000,
                'batch_size': batch_size,
                'nblocks': nblocks, # 0 turns off drift correction; default is 1
                'Th_universal': 8,
                'Th_learned': 6,
                'whitening_range': 32, # Number of nearby channels to include in the whitening
                'dmin': electrode_spacing, # 60 Check this...
                'dminx': electrode_spacing, # 60 Check this...
                'min_template_size': 10, # 10 Check this...
                'max_channel_distance': 10, # Templates farther away than this from their nearest channel will not be used.
                'nearest_chans': 7, # Number of nearest channels to consider when finding local maxima during spike detection.
                'n_pcs': 6, # default is 6, we used 3 before
                'Th_single_ch': 4.5, # Default is 6... we used 4.5 before
                'x_centers': 10, # Number of x-positions to use when determining center points
                'acg_threshold': 0.2, # Fraction of refractory period violations that are allowed in the ACG compared to baseline; used to assign "good" units.
                'ccg_threshold': 0.1 # Fraction of refractory period violations that are allowed in the CCG compared to baseline; used to perform splits and merges.
                } 
    return settings

experiment_name = '20240229C'
chunk_name = 'chunk1'
electrode_pitch = 60
device = 'cuda'
use_car = False
num_channels = 512

bin_path = '/data/data/sorted/20240229C/chunk1.bin'
results_path = '/data/data/sorted/20240229C/chunk1/kilosort4/'
csv_path = '/data/data/sorted/20240229C/chunk1.csv'

probe_path = os.path.abspath( os.path.join('/home/mike/Documents/git_repos/manookin-lab/MEA/src/utilities', '..','pipeline_utilities','kilosort', 'LITKE_512_ARRAY' + '.mat') )

probe = io.load_probe(probe_path)
settings = get_settings(num_channels=num_channels, electrode_spacing=electrode_pitch, use_drift_correction=False)
settings = {**DEFAULT_SETTINGS, **settings}

filename, data_dir, results_dir, probe = \
        set_files(settings, bin_path, probe, None, None, results_path)

data_dtype = 'int16'
device = torch.device('cuda:1')
invert_sign = True
do_CAR = use_car
tic0 = time.time()
save_extra_vars = True
ops = initialize_ops(settings, probe, data_dtype, do_CAR, invert_sign, device)

progress_bar=None
file_object=None
ops = compute_preprocessing(ops, device, tic0=tic0, file_object=file_object)
np.random.seed(1)
torch.cuda.manual_seed_all(1)
torch.random.manual_seed(1)
ops, bfile, st0 = compute_drift_correction(
    ops, device, tic0=tic0, progress_bar=progress_bar,
    file_object=file_object
    )

# Sort spikes and save results
st,tF, _, _ = detect_spikes(ops, device, bfile, tic0=tic0,
        progress_bar=progress_bar)

# Cluster spikes
clu, Wall = cluster_spikes(st, tF, ops, device, bfile, tic0=tic0,
        progress_bar=progress_bar)

# Save intermediate results.
# from pathlib import Path
# results_dir = Path('/data/data/sorted/20240229C/chunk1/')
# np.save(results_dir / 'st.npy', st)
# np.save(results_dir / 'tF.npy', tF.cpu().numpy())
# np.save(results_dir / 'Wall.npy', Wall.cpu().numpy())
# np.save(results_dir / 'clu.npy', clu)

ops, similar_templates, is_ref, est_contam_rate, kept_spikes = \
        save_sorting(ops, results_path, st, clu, tF, Wall, bfile.imin, tic0,
                     save_extra_vars=save_extra_vars)

# Post-process the results. Make sure to run this!!!



bfile = io.BinaryFiltered(
        ops['filename'], n_chan_bin, fs, NT, nt, twav_min, chan_map, 
        hp_filter=hp_filter, whiten_mat=whiten_mat, device=device, do_CAR=do_CAR,
        invert_sign=invert, dtype=dtype, tmin=tmin, tmax=tmax,
        artifact_threshold=artifact, shift=shift, scale=scale,
        file_object=file_object
        )

reader = io.BinaryRWFile(filename='/data/data/sorted/20240229C/chunk1.bin',
                         n_chan_bin=512,
                         fs=20000,
                         NT = 60000,
                         nt = 61,
                         nt0min = 20,
                         device=torch.device('cpu'),
                         dtype = 'int16',
                         write = False)

reader.n_batches
ibatch = 0
X, inds = reader.padded_batch_to_torch(ibatch, return_inds=True)

for ibatch in tqdm(np.arange(bfile.n_batches), miniters=200 if progress_bar else None, 
                        mininterval=60 if progress_bar else None):
        X = bfile.padded_batch_to_torch(ibatch, ops)

## test region
from clustering_qr import xy_up, x_centers, y_centers, get_data_cpu, cluster
xy, iC = xy_up(ops)
>>> iclust_template = st0[:,5].astype('int32')
>>> xcup, ycup = ops['xcup'], ops['ycup']
dmin = ops['dmin']
dminx = ops['dminx']
nskip = ops['settings']['cluster_downsampling']
ycent = y_centers(ops)
xcent = x_centers(ops)
nsp = st.shape[0]

# Get positions of all grouping centers
ycent_pos, xcent_pos = np.meshgrid(ycent, xcent)
ycent_pos = torch.from_numpy(ycent_pos.flatten())
xcent_pos = torch.from_numpy(xcent_pos.flatten())
# Compute distances from templates
center_distance = (
    (xy[0,:] - xcent_pos.unsqueeze(-1))**2
    + (xy[1,:] - ycent_pos.unsqueeze(-1))**2
    )
# Add some randomness in case of ties
center_distance += 1e-20*torch.rand(center_distance.shape)
# Get flattened index of x-y center that is closest to template
minimum_distance = torch.min(center_distance, 0).indices

kk=0; jj=0;
>>> ii = kk + jj*ycent.size
>>> ii
0
>>> ix = (minimum_distance == ii)
tF = torch.from_numpy(tF)
Xd, ch_min, ch_max, igood  = get_data_cpu(
                ops, xy, iC, iclust_template, tF, ycent[kk], xcent[jj], dmin=dmin,
                dminx=dminx, ix=ix
                )

Xd.shape

from clustering_qr import kmeans_plusplus
nclust = 200
seed = 1
Xg = Xd.to(device)
iclust_init =  kmeans_plusplus(Xg, niter = nclust, seed = seed, device=device)

# ops, st, clu, tF, Wall, similar_templates, is_ref, est_contam_rate = \
#     run_kilosort(settings = settings, 
#                 probe = probe,
#                 filename = bin_path,
#                 results_dir = results_path,
#                 data_dtype = 'int16', # Check the data type
#                 device = torch.device('cuda'), # device = torch.device('cuda:1'), device = torch.device('cuda')
#                 invert_sign = True, # Check this; Invert the sign of the data as expected by kilosort4 (was False)
#                 do_CAR = use_car,
#                 save_extra_vars=True)


