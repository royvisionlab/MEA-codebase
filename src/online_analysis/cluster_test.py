import sys
sys.path.append('../analysis/')

import numpy as np
# import matplotlib.pyplot as plt
import os
import platform
from sorting_hat import cluster_compact, cluster_full, get_cluster_labels, write_cluster_labels
from scipy.io import loadmat

mdic = loadmat('/Users/michaelmanookin/Documents/Data/rawdata/MEA_Data/analysis/20220712C/kilosort2/small_params.mat')

timecourse_matrix = mdic['timecourse_matrix']
acf = mdic['acf']
total_spikes = mdic['total_spikes']
cluster_id = mdic['cluster_id'].flatten()
rf_area = mdic['rf_area'].flatten()
hull_area = mdic['hull_area'].flatten()
significance_area = mdic['significance_area'].flatten()

if cluster_id.ndim == 2:
    cluster_id = cluster_id.flatten()

cell_clusters = cluster_compact(timecourse_matrix, acf, rf_area, total_spikes, variance_or_components=0.85,
        refractory_threshold=0.2, snr_threshold=2.0, isi_binning=1.0, count_threshold=300)
category_list = get_cluster_labels(cell_clusters, cluster_id)
write_cluster_labels('/Users/michaelmanookin/Documents/Data/rawdata/MEA_Data/analysis/20220712C/kilosort2/noise/online_clusters.txt',category_list)
print("Wrote automated clustering")

cell_clusters = cluster_full(timecourse_matrix, acf, hull_area, significance_area, total_spikes, variance_or_components=0.85,
        refractory_threshold=0.2, snr_threshold=2.0, isi_binning=1.0, count_threshold=300)
category_list = get_cluster_labels(cell_clusters, cluster_id)
write_cluster_labels('/Users/michaelmanookin/Documents/Data/rawdata/MEA_Data/analysis/20220712C/kilosort2/noise/full_clusters.txt',category_list)
print("Wrote full clustering")