import os, shutil
import numpy as np
import diptest
from typing import Tuple
from tqdm import tqdm
from pathlib import Path

from kilosort import CCG

SPIKE_TIMES_FILENAME = 'spike_times.npy'
SPIKE_IDENTITY_FILENAME = 'spike_clusters.npy'
CLUSTER_QUALITY_FILENAME = 'cluster_KSLabel.tsv'

MAX_SAMPLE_SIZE = 72000
MIN_SPIKE_COUNT = 300 # Minimum number of spikes required to process a cluster.

def build_cluster_quality_dict(file_path: str):
    cluster_quality_by_id = {}
    with open(file_path, 'r') as cluster_quality_file:
        cluster_quality_file.readline()
        remaining_lines = cluster_quality_file.readlines()
        for line in remaining_lines:
            data_list = line.strip('\n').split('\t')
            cluster_quality_by_id[int(data_list[0])] = data_list[1]
    return cluster_quality_by_id

def save_cluster_quality_dict(ops: dict, spike_clusters: np.ndarray, spike_times: np.ndarray, template_amplitudes: np.ndarray, results_dir: str):
    # contamination ratio
    acg_threshold = ops['settings']['acg_threshold']
    ccg_threshold = ops['settings']['ccg_threshold']
    is_ref, est_contam_rate = CCG.refract(spike_clusters, spike_times / ops['fs'],
                                          acg_threshold=acg_threshold,
                                          ccg_threshold=ccg_threshold)

    # write properties to *.tsv
    stypes = ['ContamPct', 'KSLabel'] #stypes = ['ContamPct', 'Amplitude', 'KSLabel']
    ks_labels = [['mua', 'good'][int(r)] for r in is_ref]
    props = [est_contam_rate*100, ks_labels] #props = [est_contam_rate*100, template_amplitudes, ks_labels]
    for stype, prop in zip(stypes, props):
        with open((results_dir / f'cluster_{stype}.tsv'), 'w') as f:
            f.write(f'cluster_id\t{stype}\n')
            for i,p in enumerate(prop):
                if stype != 'KSLabel':
                    f.write(f'{i}\t{p:.1f}\n')
                else:
                    f.write(f'{i}\t{p}\n')
        if stype == 'KSLabel':
            shutil.copyfile((results_dir / f'cluster_{stype}.tsv'), 
                            (results_dir / f'cluster_group.tsv'))

def kmeans_plus_plus(data: np.ndarray, k: int=2) -> np.ndarray:
    # Choose the first center uniformly at random.
    centers = [data[np.random.choice(data.shape[0]),:]]
    for i in range(1,k):
        # Compute the distance from each data point to the nearest center.
        distances = np.array([min([np.linalg.norm(x-c) for c in centers]) for x in data])
        # Choose the next center from the data points with probability proportional to the square of the distance to the nearest center.
        centers.append(data[np.random.choice(data.shape[0],p=distances**2/np.sum(distances**2),replace=False)])
    return np.array(centers)

def kmeans_cluster(data: np.ndarray, k: int=2) -> Tuple[np.ndarray, np.ndarray]:
    centers = kmeans_plus_plus(data, k)
    while True:
        # Assign each data point to the nearest center.
        labels = np.argmin(np.array([np.linalg.norm(data-c,axis=1) for c in centers]),axis=0)
        # Recompute the cluster centers as the mean of the data points assigned to the cluster.
        new_centers = np.array([data[labels == i].mean(axis=0) for i in range(k)])
        if np.all(centers == new_centers):
            break
        centers = new_centers
    return centers, labels

def evaluate_cluster(Xs: np.ndarray):
    try:
        if np.ndim(Xs) == 3:
            # Combine the last two dimensions of the numpy array.
            Xs = Xs.reshape(Xs.shape[:-2] + (-1,))
        # K-means clustering.
        _, labels = kmeans_cluster(data=Xs, k=2)
        # Labels need to be {-1,1}
        labels = 2.0 * labels.astype(float) - 1.0
        # Get the label weights.
        w = np.ones((Xs.shape[0],1))
        w[labels>0] = np.mean(labels<0)
        w[labels<0] = np.mean(labels>0)
        # Compute the weighted covariance matrix.
        CC = Xs.T @ (Xs * w)
        # Ridge regression.
        CC = CC + 0.01 * np.eye( CC.shape[0] )
        b = np.linalg.solve( CC, labels @ (Xs * w) )
        xproj = Xs @ b
        # Compute the dip statistic and p-value.
        if xproj.shape[0] > MAX_SAMPLE_SIZE:
            idx = np.random.choice(xproj.shape[0],MAX_SAMPLE_SIZE,replace=False)
            xproj = xproj[idx]
        _, pval = diptest.diptest(xproj)
    except Exception as error:
        print(error)
        return 1.0, np.zeros(Xs.shape[0])
    return pval, labels

def split_cluster(Xs: np.ndarray, spike_idx: np.ndarray, out: list):
    pval, labels = evaluate_cluster(Xs[spike_idx])
    if pval < 0.05:
        u_labels = np.unique(labels)
        for ii in range(2):
            idx = np.where(labels == u_labels[ii])[0]
            split_cluster(Xs, spike_idx[idx], out)
    else:
        u_labels = np.unique(labels)
        for ii, label in enumerate(u_labels):
            idx = np.where(labels == label)[0]
            out.append(spike_idx[idx])

def clean_clusters(spike_clusters: np.ndarray, 
                   pc_features: np.ndarray,
                   is_cleaner: bool=True) -> np.ndarray:
    unique_clusters = np.unique(spike_clusters) 
    split_clusters = dict()
    for target_cluster in tqdm(unique_clusters, desc='Evaluating clusters', unit='clusters', leave=False):
        try:
            target_spikes = np.where(spike_clusters == target_cluster)[0]
            if target_spikes.shape[0] < MIN_SPIKE_COUNT:
                continue
            Xs = pc_features[target_spikes,...]
            if np.all(Xs == 0.0):
                continue
            if is_cleaner:
                Xs = np.transpose(Xs, (0,2,1))
            # Combine the last two dimensions of the numpy array.
            Xs = Xs.reshape(Xs.shape[:-2] + (-1,))
            # Evaluate the cluster.
            # score, dip, pval, labels = evaluate_single_cluster(Xs)
            pval, labels = evaluate_cluster(Xs)
            if pval < 0.05:
                # print(f'Cluster {target_cluster} is bimodal.')
                # Split the cluster.
                # spike_indices = spike_clusters[target_spikes]
                u_labels = np.unique(labels)
                split_list = np.zeros(2, dtype=object)
                for ii, label in enumerate(u_labels):
                    idx = np.where(labels == label)[0]
                    split_list[ii] = target_spikes[idx]
                split_clusters[target_cluster] = split_list
        except Exception as error:
            print(error)
        clean_clusters = np.copy(spike_clusters)
        if len(split_clusters) > 0:
            cluster_count = np.max(spike_clusters) + 1
            print(f'Split {len(split_clusters)} clusters.')
            for key, value in split_clusters.items():
                print(f'Cluster {key} has been split into {len(value)} clusters.')
                for ii in range(1,len(value)):
                    clean_clusters[value[ii]] = cluster_count.astype(np.int32)
                    cluster_count += 1
    return clean_clusters


class KilosortCleaner(object):
    def __init__(self):
        self.spike_times = None
        self.spike_clusters = None
        self.templates = None
        self.spike_templates = None
        self.pc_features = None
        self.pc_feature_ind = None
        self.ops = None

    def load_npy(self, file_path):
        self.spike_times = np.load(os.path.join(file_path, 'spike_times.npy'))
        self.spike_clusters = np.load(os.path.join(file_path, 'spike_clusters.npy'))
        self.templates = np.load(os.path.join(file_path, 'templates.npy'))
        self.spike_templates = np.load(os.path.join(file_path, 'spike_templates.npy'))
        self.pc_features = np.load(os.path.join(file_path, 'pc_features.npy'))
        self.pc_feature_ind = np.load(os.path.join(file_path, 'pc_feature_ind.npy'))
        self.ops = np.load(os.path.join(file_path, 'ops.npy'), allow_pickle=True).item()


results_dir = '/data/data/sorted/20240229C/chunk1/kilosort4_clean'
results_dir = Path(results_dir)
np.save((results_dir / 'spike_times.npy'), spike_times)
np.save((results_dir / 'spike_clusters.npy'), spike_clusters)

# Post-process the results. Make sure to run this!!!
