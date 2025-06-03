import argparse
import csv
import numpy as np
import os
import torch
from kilosort import run_kilosort, io
import config as cfg
# import pathlib

os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "expandable_segments:True"

SPIKE_TIMES_FILENAME = 'spike_times.npy'
SPIKE_IDENTITY_FILENAME = 'spike_clusters.npy'
CLUSTER_QUALITY_FILENAME = 'cluster_KSLabel.tsv'

def build_cluster_quality_dict(filepath):
    cluster_quality_by_id = {}
    with open(filepath, 'r') as cluster_quality_file:
        cluster_quality_file.readline()
        remaining_lines = cluster_quality_file.readlines()
        for line in remaining_lines:
            data_list = line.strip('\n').split('\t')
            cluster_quality_by_id[int(data_list[0])] = data_list[1]
    return cluster_quality_by_id

def save_cluster_quality_dict(cluster_quality_dict, filepath):
    stype = 'KSLabel'
    with open(filepath, 'w') as f:
        f.write(f'cluster_id\t{stype}\n')
        for i,p in cluster_quality_dict.items():
            f.write(f'{i}\t{p}\n')

def get_settings(
        num_channels: int=512, 
        sample_rate: int=20000,
        batch_size: int=60000,
        use_drift_correction: bool=False,
        threshold_universal: int=8,
        threshold_learned: int=6,
        electrode_spacing: float=60.0, 
        min_template_size: int=10,
        max_channel_distance: int=10,
        max_peels: int=100,
        nearest_chans: int=19, # First ring:7, second ring:19, third ring: 37
        num_pcs: int=6,
        threshold_single_ch: float=3.0,
        x_centers: int=10,
        acg_threshold: float=0.2,
        ccg_threshold: float=0.2,
        whitening_range: int=37):
    
    if use_drift_correction:
        nblocks = 1
    else:
        nblocks = 0
    settings = {'n_chan_bin': num_channels,
                'fs': sample_rate,
                'batch_size': batch_size,
                'nblocks': nblocks, # 0 turns off drift correction; default is 1
                'Th_universal': threshold_universal,
                'Th_learned': threshold_learned,
                'nskip': 25, # Batch stride for computing whitening matrix.
                'whitening_range': whitening_range, # Number of nearby channels to include in the whitening (default: 32)
                'highpass_cutoff': 300, # Critical frequency for highpass Butterworth filter applied to data.
                'sig_interp': 20, # Approximate spatial smoothness scale in units of microns.
                'dmin': electrode_spacing, # 60 Check this...
                'dminx': electrode_spacing, # 60 Check this...
                'min_template_size': min_template_size, # 10 Check this...
                'max_channel_distance': max_channel_distance, # Templates farther away than this from their nearest channel will not be used.
                'max_peels': max_peels, # Number of iterations to do over each batch of data in the matching pursuit step. More iterations may detect more overlapping spikes.
                'nearest_chans': nearest_chans, # Number of nearest channels to consider when finding local maxima during spike detection.
                'templates_from_data': True, # Indicates whether spike shapes used in universal templates should be estimated from the data or loaded from the predefined templates.
                'n_templates': 6, # Number of single-channel templates to use for the universal templates (only used if templates_from_data is True).
                'n_pcs': num_pcs, # default is 6, we used 3 before; Number of single-channel PCs to use for extracting spike features (only used if templates_from_data is True).
                'Th_single_ch': threshold_single_ch, # Default is 6... we used 4.5 before
                'x_centers': x_centers, # Number of x-positions to use when determining center points
                'acg_threshold': acg_threshold, # Fraction of refractory period violations that are allowed in the ACG compared to baseline; used to assign "good" units.
                'ccg_threshold': ccg_threshold, # Fraction of refractory period violations that are allowed in the CCG compared to baseline; used to perform splits and merges.
                'cluster_downsampling': 20, # Inverse fraction of nodes used as landmarks during clustering (can be 1, but that slows down the optimization).
                'duplicate_spike_ms': 0.25 # Time in ms for which subsequent spikes from the same cluster are assumed to be artifacts. A value of 0 disables this step.
                } 
    return settings

def postprocess(experiment_name, sorter_name, results_path, csv_path, SAVE_PATH):
    file_names = list()
    file_sizes = list()
    with open(csv_path, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in reader:
            s = row[0].split('/')
            file_names.append(s[-1])
            file_sizes.append(row[1])
    # Get the spike times and cluster identities
    spike_times_filepath = os.path.join(results_path, SPIKE_TIMES_FILENAME)
    spike_identity_filepath = os.path.join(results_path, SPIKE_IDENTITY_FILENAME)
    cluster_quality_dict_filepath = os.path.join(results_path, CLUSTER_QUALITY_FILENAME)
    spike_times_vector = np.load(spike_times_filepath)
    spike_identity_vector = np.load(spike_identity_filepath)
    quality_by_cluster_id = build_cluster_quality_dict(cluster_quality_dict_filepath)
    file_sizes = np.array(file_sizes).astype(int)
    count_offset = 0
    for file_count in range(len(file_names)):
        spike_idx = np.where((spike_times_vector >= count_offset) & (spike_times_vector < count_offset + file_sizes[file_count]))[0]
        if len(spike_idx) > 0:
            # Make the directory if it doesn't exist.
            if not os.path.exists(os.path.join(SAVE_PATH, experiment_name, file_names[file_count], sorter_name)):
                os.makedirs(os.path.join(SAVE_PATH, experiment_name, file_names[file_count], sorter_name))
            spike_times = spike_times_vector[spike_idx] - count_offset
            spike_identity = spike_identity_vector[spike_idx]
            np.save(os.path.join(SAVE_PATH, experiment_name, file_names[file_count], sorter_name, SPIKE_TIMES_FILENAME), spike_times)
            np.save(os.path.join(SAVE_PATH, experiment_name, file_names[file_count], sorter_name, SPIKE_IDENTITY_FILENAME), spike_identity)
            save_cluster_quality_dict(quality_by_cluster_id, os.path.join(SAVE_PATH, experiment_name, file_names[file_count], sorter_name, CLUSTER_QUALITY_FILENAME))
        count_offset += file_sizes[file_count]



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Run kilosort4 and process the output.')
    parser.add_argument('experiment_name', type=str, help='folder containing Kilosort outputs, must already exist')
    parser.add_argument('chunk_name', type=str, help='name of data file/folder (e.g. data020)')
    parser.add_argument('-a','--algorithm', default='kilosort4', type=str, help='Sorting algorithm used (yass or kilosort2)')
    parser.add_argument('-e', '--electrode_pitch', default=60, type=int, help='Electrode array pitch')
    parser.add_argument('-d','--device', default='cuda', type=str, help='The device to run the sorter on (default is cuda; e.g. --device "cuda:1" for torch.device("cuda:1") or --device "cpu" for torch.device("cpu")')
    parser.add_argument('-c','--use_car', action='store_true', help='Whether to get spike times or spike counts')
    parser.add_argument('-n', '--nearest_chans', default=19, type=int, help='Nearest channels to consider when finding local maxima during spike detection')
    parser.add_argument('-b', '--batch_size', default=60000, type=int, help='Batch size for processing data')
    parser.add_argument('-at', '--acg_threshold', default=0.2, type=float, help='Threshold for ACG')
    parser.add_argument('-ct', '--ccg_threshold', default=0.2, type=float, help='Threshold for CCG')
    parser.add_argument('-w', '--whitening_range', default=37, type=int, help='Range for whitening')
    parser.add_argument('-p', '--max_peels', default=200, type=int, help='Maximum number of peels')
    parser.add_argument('-cc', '--clear_cache', action='store_true', help='Clear the cache')

    args = parser.parse_args()

    sorter_name = 'kilosort4'

    if ('mmfs1' in os.getcwd()) or ('gscratch' in os.getcwd()): # Hyak
        SORT_PATH = '/gscratch/scrubbed/retina/data/sorted/'
    else:
        SORT_PATH, _ = cfg.get_data_paths()
    SAVE_PATH = cfg.get_save_path()
    experiment_name = args.experiment_name
    chunk_name = args.chunk_name
    electrode_pitch = args.electrode_pitch
    device = args.device
    use_car = args.use_car
    nearest_chans = args.nearest_chans
    batch_size = args.batch_size
    acg_threshold = args.acg_threshold
    ccg_threshold = args.ccg_threshold
    whitening_range = args.whitening_range
    max_peels = args.max_peels
    clear_cache = args.clear_cache

    print('Common Average Referencing settting is ' + str(use_car))

    # probe_path = '/home/mike/Documents/git_repos/manookin-lab/MEA/src/pipeline_utilities/kilosort/LITKE_512_ARRAY.mat'
    bin_path = os.path.join(SORT_PATH, experiment_name, chunk_name + '.bin') #'/data/data/sorted/20240401C/chunk1.bin'
    results_path = os.path.join(SAVE_PATH, experiment_name, chunk_name, sorter_name) #'/data/data/sorted/20240401C/chunk1/kilosort4/'
    csv_path = os.path.join(SAVE_PATH, experiment_name, chunk_name + '.csv') #'/data/data/sorted/20240401C/chunk1.csv'

    # pathlib.Path(__file__).parent.resolve()
    if electrode_pitch == 30:
        num_channels = 519
        probe_name = 'LITKE_519_ARRAY_30UM'
    elif electrode_pitch == 120:
        num_channels = 519
        probe_name = 'LITKE_519_ARRAY_120UM'
    elif electrode_pitch == 60:
        num_channels = 512
        probe_name = 'LITKE_512_ARRAY'
    elif electrode_pitch == 200:
        num_channels = 60
        probe_name = 'MCS_60_ARRAY_200UM'
    else:
        num_channels = 512
        probe_name = 'LITKE_512_ARRAY'
    dir_path = os.path.dirname( os.path.realpath(__file__) )
    probe_path = os.path.abspath( os.path.join(dir_path, '..','pipeline_utilities','kilosort', probe_name + '.mat') )

    probe = io.load_probe(probe_path)
    settings = get_settings(num_channels=num_channels, 
            electrode_spacing=electrode_pitch, 
            max_peels=max_peels,
            nearest_chans=nearest_chans,
            use_drift_correction=False,
            batch_size=batch_size,
            acg_threshold=acg_threshold,
            ccg_threshold=ccg_threshold,
            whitening_range=whitening_range
            )

    ops, st, clu, tF, Wall, similar_templates, is_ref, est_contam_rate, kept_spikes = \
        run_kilosort(settings = settings, 
                    probe = probe,
                    filename = bin_path,
                    results_dir = results_path,
                    data_dtype = 'int16', # Check the data type
                    device = torch.device( device ), # device = torch.device('cuda:1'), device = torch.device('cuda')
                    invert_sign = True, # Check this; Invert the sign of the data as expected by kilosort4 (was False)
                    do_CAR = use_car, # Common average referencing
                    save_extra_vars=False,
                    clear_cache=clear_cache)
    
    # imin = 0
    # run_kilosort.save_sorting(ops, results_path, st, clu, tF, Wall, imin, save_extra_vars=False)

    # Post-process the results.
    postprocess(experiment_name, sorter_name, results_path, csv_path, SAVE_PATH)


# Follow this through the stack. This is how to use predefined templates.
# parameters.templates_from_data = True
# python -m pip install kilosort[gui] --ignore-installed --force-reinstall
# import inspect
# import re

# ret2 = inspect.getsourcelines(run_kilosort)
# ret_pattern = r'return\s*(.*)\n*$'

# combined_ret = "".join(ret2[0][-2:]) # Combine the last two lines of the function

# len(combined_ret.split(','))
# out2 = re.findall(ret_pattern, combined_ret)

# pip install -e ".[gui]"

