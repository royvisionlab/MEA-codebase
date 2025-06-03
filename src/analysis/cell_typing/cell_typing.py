
import numpy as np
import hdf5storage
import os
import argparse
from sklearn.cluster import KMeans
from sklearn.decomposition import TruncatedSVD, PCA
from scipy.ndimage import gaussian_filter
from typing import Tuple

# from numpy import mean
# from numpy import std
# from sklearn.model_selection import train_test_split
# from sklearn.model_selection import cross_val_score
# from sklearn.model_selection import RepeatedStratifiedKFold
# from sklearn.pipeline import Pipeline
# from sklearn.linear_model import LogisticRegression
# from matplotlib import pyplot
 
# get a list of models to evaluate
# def get_models():
# 	models = dict()
# 	for i in range(1,20):
# 		steps = [('svd', TruncatedSVD(n_components=i)), ('m', LogisticRegression())]
# 		models[str(i)] = Pipeline(steps=steps)
# 	return models
 
# # evaluate a give model using cross-validation
# def evaluate_model(model, X, y):
# 	cv = RepeatedStratifiedKFold(n_splits=10, n_repeats=3, random_state=1)
# 	scores = cross_val_score(model, X, y, scoring='accuracy', cv=cv, n_jobs=-1, error_score='raise')
# 	return scores

# def combine_dims(a, start=0, count=2):
#     """ Reshapes numpy array a by combining count dimensions, 
#         starting at dimension index start """
#     s = a.shape
#     return np.reshape(a, s[:start] + (-1,) + s[start+count:])

# Return the peak (positive or negative) value of the first max of the time filter.
def find_timecourse_peak(timecourse):
    min_val = np.min(timecourse)
    max_val = np.max(timecourse)
    # Determine which came first.
    min_idx = np.where(timecourse == min_val)[0]
    max_idx = np.where(timecourse == max_val)[0]
    if max_idx[0] < min_idx[0]:
        return max_val
    else:
        return min_val

def initial_polarity_split(timecourse_matrix):
    yellow_pc = PCA(n_components=2).fit_transform(timecourse_matrix[:,:,0])
    blue_pc = PCA(n_components=2).fit_transform(timecourse_matrix[:,:,2])
    idx_blue = np.where((blue_pc[:,0] > 0.3) & (yellow_pc[:,1] < -0.5))[0]
    idx_off = np.where(yellow_pc[:,0] < 0.0)[0]
    idx_on = np.where(yellow_pc[:,0] >= 0.0)[0]
    return idx_blue, idx_off, idx_on

def get_weighted_pca(X, variance_or_components=0.9):
    my_model = PCA(n_components=variance_or_components)
    pc = my_model.fit_transform(X)
    # Weight by explained variance.
    # w = my_model.explained_variance_ratio_ / np.sum(my_model.explained_variance_ratio_)
    # Normalize.
    pc /= np.nanstd(pc)
    # for i in range(pc.shape[0]):
    #     pc[i,:] *= w
    return pc

def get_components_and_weights(X, variance_or_components=0.9)->Tuple[np.ndarray, np.ndarray]:
    """Run PCA on input and get normalized PC and wts

    Args:
        X (np.ndarray): input matrix
        variance_or_components (float, optional): PCA n_components parameter, can specify # components or % variance explained. Defaults to 0.9.

    Returns:
        pc (np.ndarray): normalized principle components matrix
        weights (np.ndarray): normalized weights of principle components
    """
    my_model = PCA(n_components=variance_or_components)
    pc = my_model.fit_transform(X - X.mean(axis=0))
    # Normalize.
    pc /= np.std(pc)
    # Weight by explained variance.
    total_variance_explained = np.sum(my_model.explained_variance_ratio_)
    weights = my_model.explained_variance_ratio_ / total_variance_explained
    return pc, weights

def cluster_cells(hull_area, gauss_params, significance_area, timecourse_matrix, acf, num_clusters):
    # Create the design matrix.
    X0 = np.zeros((hull_area.shape[0], 5))
    X0[:,0] = hull_area / np.std(hull_area)
    X0[:,1] = np.mean(gauss_params[:,(4,5)], axis=1) / np.std(gauss_params[:,(4,5)])
    
    # TODO look into why this is often all 0's and leads to divide by 0 error
    X0[:,2] = significance_area / np.std(significance_area)
    yellow_t = get_weighted_pca(timecourse_matrix[:,:,0]) #my_model.fit_transform(timecourse_matrix[:,:,0])
    blue_t = get_weighted_pca(timecourse_matrix[:,:,2])
    acf_t = get_weighted_pca(acf)
    reduced = np.append(yellow_t, X0, axis=1)
    reduced = np.append(reduced, blue_t, axis=1)
    reduced = np.append(reduced, acf_t, axis=1)
    
    # Replace NaNs with 0.
    reduced = np.nan_to_num(reduced, nan=0.0, posinf=0.0, neginf=0.0)
    
    # Determine the optimal number of clusters using the elbow method.
    num_clusters = get_num_clusters(reduced)

    kmeans = KMeans(n_clusters=num_clusters, init='k-means++', max_iter=30000, n_init=1, verbose=False).fit(reduced)
    labels = kmeans.labels_
    return labels

def compute_rank1(sta):
    nt, ny, nx = sta.shape[1:4]
    r2 = np.zeros((sta.shape[0],))
    for i in range(sta.shape[0]):
        s = sta[i,:,:,:,0]
        # Pass the STA through a Gaussian spatial filter.
        s = gaussian_filter(s, sigma=0.75)
        s_mat = s.reshape((nt, ny*nx))
        # Extract the spatiotemporal features using SVD.
        space, singular_values, time = np.linalg.svd(s_mat.T, full_matrices=False)
        # Get the rank-1 approximation
        r1 = np.outer(time[0,:], space[:,0]) * singular_values[0]
        # Get the r^2 value
        coeff = np.corrcoef(s_mat.ravel(), r1.ravel())
        r2[i] = coeff[0,1]
    return r2

def get_num_clusters(pc: np.ndarray, cluster_threshold: float=0.9)->int:
    """! Get optimal number of KMeans clusters k using elbow method
    1. Test a range of clusters to
    2. Generate WSS (within cluster sum of squares, aka inertia) distribution
    3. Make normalized cumulative WSS distribution
    4. Select optimal k as the lowest k whose normalized cumulative WSS > threshold

    Args:
        @param pc (np.ndarray): Input matrix of principal components [n_samples x n_components]
        @param cluster_threshold (float, optional): Defaults to 0.9.

    Returns:
        @return int: optimal number of KMeans clusters
    """
    distortions = []
    # Test number of clusters from range (1, min(n_samples, n_components)) as n_components must be <= n_samples
    K = range(1, np.min([pc.shape[1], pc.shape[0]]) + 1)

    for k in K:
        # print(k)
        mdl = KMeans(n_clusters=k, n_init=10)
        mdl.fit(pc)
        distortions.append(mdl.inertia_)

    # Normalize the distortions.
    distortions -= np.min(distortions)
    distortions /= np.sum(distortions)
    c_dist = np.cumsum(distortions)
    # Get the optimal k value (first above threshold)
    return np.argwhere(c_dist > cluster_threshold)[0][0] + 1 # Adding 1 to get value of k

def cluster_experiment(output_path, 
            sta, 
            timecourse_matrix, 
            acf, 
            significance_maps, 
            hull_area,
            spike_count, 
            gauss_params, 
            cluster_id, 
            refractory_threshold=0.1,
            snr_threshold=2.0):
    """
    Computes STA by convolving the spiketrain with the stimulus. This 
    implementation uses straight up numpy.

    Parameters:
        stimulus: stimulus tensor of size t, y, x, color (8 bit ints).
        binned_spikes: spiketrain matrix of size t,cells
        stride: the stride of convolution (interval)
        depth: number of lags for convolution
        norm: boolean whether to normalize STA for nice plotting
    Returns:
        STA tensor, indexed by cell,t,y,x,color channel
    """

    # Calculate the goodness of fit for the rank1 RF approximation.
    r2 = compute_rank1(sta)

    # Get the space-time variant cells; low R^2
    idx_var = np.where(r2 < 0.0)[0]

    # Calculate SNR.
    np.seterr(invalid='ignore')
    snr = np.divide(np.std(timecourse_matrix[:,1:31,(0,2)],axis=1), np.std(timecourse_matrix[:,31:,(0,2)],axis=1))
    # snr = np.std(timecourse_matrix[:,1:31,(0,2)],axis=1) / np.std(timecourse_matrix[:,31:,(0,2)],axis=1)
    
    # Replace NaNs with 0 as there may be all 0 timecourses.
    snr = np.nan_to_num(snr, nan=0.0, posinf=0.0, neginf=0.0) 
    snr = np.nanmax(snr,axis=1)

    # acf has 1ms bins, n_bin_max is index of bins to count spikes upto. 2 means first two bins will be included.
    n_bin_max = 4
    # Calculate the percentage of spikes that violate refractoriness (<2ms isi)
    pct_refractory = np.sum(acf[:,:n_bin_max], axis=1) * 100

    # May have to play with the cutoff
    idx_bad = np.argwhere((snr < snr_threshold) | (pct_refractory > refractory_threshold))[:,0]

    # Get the area of significance in the maps.
    significance_area = np.sum(significance_maps, axis=(1,2))

    idx_blue, idx_off, idx_on = initial_polarity_split(timecourse_matrix)

    # Remove any potentially invariant cells from the other arrays.
    bool_array = np.where(np.in1d(idx_var, idx_bad))[0]
    idx_var = np.delete(idx_var, bool_array)
    bool_array = np.where(np.in1d(idx_blue, idx_var))[0]
    idx_blue = np.delete(idx_blue, bool_array)
    bool_array = np.where(np.in1d(idx_off, idx_var))[0]
    idx_off = np.delete(idx_off, bool_array)
    bool_array = np.where(np.in1d(idx_on, idx_var))[0]
    idx_on = np.delete(idx_on, bool_array)

    # Remove any blue cells from the other arrays.
    bool_array = np.where(np.in1d(idx_off, idx_blue))[0]
    idx_off = np.delete(idx_off, bool_array)
    bool_array = np.where(np.in1d(idx_on, idx_blue))[0]
    idx_on = np.delete(idx_on, bool_array)

    # Remove the bad cells.
    bool_array = np.where(np.in1d(idx_blue, idx_bad))[0]
    idx_blue = np.delete(idx_blue, bool_array)
    bool_array = np.where(np.in1d(idx_off, idx_bad))[0]
    idx_off = np.delete(idx_off, bool_array)
    bool_array = np.where(np.in1d(idx_on, idx_bad))[0]
    idx_on = np.delete(idx_on, bool_array)

    # Create a category list.
    category_list = list()

    if len(idx_blue) > 2:
        blue_labels = cluster_cells(hull_area[idx_blue],gauss_params[idx_blue,...], significance_area[idx_blue], timecourse_matrix[idx_blue,...], acf[idx_blue,...],2)
        for i, label in enumerate(blue_labels):
            idx = idx_blue[i]
            category_list.append(str(cluster_id[idx]) + '  All/Blue/nc' + str(label))

        cat_offset = np.max(blue_labels) + 1
    else:
        cat_offset = 1

    # Off cells
    off_labels = cluster_cells(hull_area[idx_off],gauss_params[idx_off,...], significance_area[idx_off], timecourse_matrix[idx_off,...], acf[idx_off,...],10)
    for i, label in enumerate(off_labels):
        idx = idx_off[i]
        category_list.append(str(cluster_id[idx]) + '  All/Off/nc' + str(label + cat_offset))

    cat_offset += np.max(off_labels) + 1

    on_labels = cluster_cells(hull_area[idx_on],gauss_params[idx_on,...], significance_area[idx_on], timecourse_matrix[idx_on,...], acf[idx_on,...],10)
    for i, label in enumerate(on_labels):
        idx = idx_on[i]
        category_list.append(str(cluster_id[idx]) + '  All/On/nc' + str(label + cat_offset))

    for idx in idx_var:
        category_list.append(str(cluster_id[idx]) + '  All/stvar')

    # Bad cells.
    count_threshold = 10000
    idx_bad_low = idx_bad[np.where(spike_count[idx_bad] <= count_threshold)[0]]
    idx_bad_high = idx_bad[np.where(spike_count[idx_bad] > count_threshold)[0]]
    for idx in cluster_id[idx_bad_high]:
        category_list.append(str(idx) + '  All/lowsnr')

    for idx in cluster_id[idx_bad_low]:
        category_list.append(str(idx) + '  All/crap')

    # Kick out the clusters to a classifications file.
    with open(output_path, 'w') as f:
        for line in category_list:
            f.write(line)
            f.write('\n')
    return category_list


def compute_snr(timecourse_matrix):
    half_time = np.floor(timecourse_matrix.shape[1]/2).astype(int) + 1
    np.seterr(all='ignore')
    snr = np.divide(np.std(timecourse_matrix[:,1:half_time,(0,2)], axis=1), np.std(timecourse_matrix[:,half_time:,(0,2)], axis=1))
    snr[np.isnan(snr)] = 0.0 # Set NaN's to zero.
    snr = np.nanmax(snr, axis=1)
    return snr

def compute_refractory_pct(isi, isi_binning=0.5):
    """
    Calculate the percentage of spikes that violate refractoriness (<2ms isi) 
    Parameters:
        acf: interspike interval distribution (sum=1)
        acf_binning: bin width of isi distribution in msec (default: 0.5)
    """
    num_bins = np.floor(2.0 / isi_binning).astype(int)
    return np.sum(isi[:,:num_bins], axis=1) * 100

def find_bad_cells(timecourse_matrix, isi, total_spikes, refractory_threshold=0.1, 
            snr_threshold=2.0, isi_binning=0.5, count_threshold=300)->np.ndarray:
    """Get index of cells that violate refractoriness and below snr_threshold

    Args:
        timecourse_matrix (_type_): _description_
        isi (bool): _description_
        total_spikes (_type_): _description_
        refractory_threshold (float, optional): _description_. Defaults to 0.1.
        snr_threshold (float, optional): _description_. Defaults to 2.0.
        isi_binning (float, optional): _description_. Defaults to 0.5.
        count_threshold (int, optional): _description_. Defaults to 300.

    Returns:
        np.ndarray: index of bad cells
    """
    # Calculate SNR.
    snr = compute_snr(timecourse_matrix)
    # Calculate the percentage of spikes that violate refractoriness (<2ms isi)
    pct_refractory = compute_refractory_pct(isi, isi_binning)
    # May have to play with the cutoff
    idx_bad = np.argwhere((snr < snr_threshold) | (pct_refractory > refractory_threshold) | (total_spikes < count_threshold))[0]
    return idx_bad

def get_cluster_labels(cell_clusters, cluster_id):
    cluster_id = cluster_id.flatten()
    category_list = list()
    for count, value in enumerate(cell_clusters):
        if value == -2:
            category_list.append(str(cluster_id[count]) + '  All/Crap')
        elif value == -1:
            category_list.append(str(cluster_id[count]) + '  All/good/lowSNR')
        elif value < 100:
            category_list.append(str(cluster_id[count]) + '  All/good/RGC/nc' + str(value))
        else:
            category_list.append(str(cluster_id[count]) + '  All/good/Amacrine/nc' + str(value))
    return category_list

def write_cluster_labels(output_path, category_list):
    with open(output_path, 'w') as f:
        for line in category_list:
            f.write(line)
            f.write('\n')

def cluster_compact(timecourse_matrix, isi, rf_area, total_spikes, variance_or_components=0.85,
            refractory_threshold=0.2, snr_threshold=2.0, isi_binning=0.5, count_threshold=300):
    cell_clusters = np.zeros((timecourse_matrix.shape[0],), dtype=np.int32)
    # Find the bad cells
    idx_bad = find_bad_cells(timecourse_matrix, isi, total_spikes, refractory_threshold=refractory_threshold, snr_threshold=snr_threshold, isi_binning=isi_binning, count_threshold=count_threshold)
    cell_clusters[idx_bad] = -1
    # Find the good cells.
    idx_good = np.argwhere(cell_clusters == 0)[:,0]
    # Initialize your design matrix
    pc = rf_area[idx_good] / np.nanstd(rf_area[idx_good])
    if pc.ndim == 1:
        pc = pc[:,np.newaxis]
    # sample_weight = np.ones((1,))
    # Run SVD on the time courses and ISI distributions.
    pc_red = get_weighted_pca(timecourse_matrix[idx_good,:,0], variance_or_components=variance_or_components)
    pc_blue = get_weighted_pca(timecourse_matrix[idx_good,:,2], variance_or_components=variance_or_components)
    pc_isi = get_weighted_pca(isi[idx_good,:], variance_or_components=variance_or_components)
    # pc_red, w_red = get_components_and_weights(timecourse_matrix[idx_good,:,0], variance_or_components=variance_or_components)
    # pc_blue, w_blue = get_components_and_weights(timecourse_matrix[idx_good,:,2], variance_or_components=variance_or_components)
    # pc_isi, w_isi = get_components_and_weights(isi[idx_good,:], variance_or_components=variance_or_components)
    pc = np.append(pc, pc_red, axis=1)
    pc = np.append(pc, pc_blue, axis=1)
    pc = np.append(pc, pc_isi, axis=1)
    # sample_weight = np.append(sample_weight, w_red, axis=0)
    # sample_weight = np.append(sample_weight, w_blue, axis=0)
    # sample_weight = np.append(sample_weight, w_isi, axis=0)
    # Use the elbow method to find the optimal number of clusters for KMeans.
    num_clusters = get_num_clusters(pc)
    # Run KMeans
    kmeans = KMeans(n_clusters=num_clusters, init='k-means++', max_iter=3000, n_init=1, verbose=False).fit(pc)
    labels = kmeans.labels_
    cell_clusters[idx_good] = labels
    return cell_clusters

def cluster_full(timecourse_matrix, isi, hull_area, significance_area, total_spikes, variance_or_components=0.85,
            refractory_threshold=0.2, snr_threshold=2.0, isi_binning=0.5, count_threshold=300):
    cell_clusters = np.zeros((timecourse_matrix.shape[0],), dtype=np.int32)
    # Find the bad cells
    idx_bad = find_bad_cells(timecourse_matrix, isi, total_spikes, refractory_threshold=refractory_threshold, snr_threshold=snr_threshold, isi_binning=isi_binning, count_threshold=count_threshold)
    cell_clusters[idx_bad] = -1
    # Find the good cells.
    idx_good = np.argwhere(cell_clusters == 0)[:,0]
    # Initialize your design matrix
    pc = hull_area[idx_good] / np.nanstd(hull_area[idx_good])
    if pc.ndim == 1:
        pc = pc[:,np.newaxis]
    # Get the ratio of hull area to significance area (should be different than one for spotty RFs).
    rf_ratio = significance_area[idx_good] / hull_area[idx_good]
    rf_ratio /= np.nanstd(rf_ratio)
    if rf_ratio.ndim == 1:
        rf_ratio = rf_ratio[:,np.newaxis]
    pc = np.append(pc, rf_ratio, axis=1)
    # Run SVD on the time courses and ISI distributions.
    pc_red = get_weighted_pca(timecourse_matrix[idx_good,:,0], variance_or_components=variance_or_components)
    pc_blue = get_weighted_pca(timecourse_matrix[idx_good,:,2], variance_or_components=variance_or_components)
    pc_isi = get_weighted_pca(isi[idx_good,:], variance_or_components=variance_or_components)
    pc = np.append(pc, pc_red, axis=1)
    pc = np.append(pc, pc_blue, axis=1)
    pc = np.append(pc, pc_isi, axis=1)
    # Use the elbow method to find the optimal number of clusters for KMeans.
    num_clusters = get_num_clusters(pc)
    # Run KMeans
    kmeans = KMeans(n_clusters=num_clusters, init='k-means++', max_iter=3000, n_init=1, verbose=False).fit(pc)
    labels = kmeans.labels_
    cell_clusters[idx_good] = labels
    return cell_clusters

def elbow_cluster(timecourse_matrix, acf, hull_area, significance_area, spike_count, cluster_id,
            refractory_threshold=0.1,
            snr_threshold=2.0):
    X = np.zeros((timecourse_matrix.shape[0],timecourse_matrix.shape[1]*2+acf.shape[1]))
    X[:,:timecourse_matrix.shape[1]] = timecourse_matrix[:,:,0]
    X[:,timecourse_matrix.shape[1]:2*timecourse_matrix.shape[1]] = timecourse_matrix[:,:,2]
    X[:,2*timecourse_matrix.shape[1]:] = acf

    my_model = PCA(0.9)
    pc = my_model.fit_transform(X)

    X0 = np.zeros((timecourse_matrix.shape[0],2))
    X0[:,0] = hull_area / np.std(hull_area)
    X0[:,1] = significance_area / np.std(significance_area)
    pc = np.append(pc, X0, axis=1)

    distortions = []
    K = range(1, pc.shape[1])
    for k in K:
        mdl = KMeans(n_clusters=k)
        mdl.fit(pc)
        distortions.append(mdl.inertia_)
    # Normalize the distortions.
    distortions -= np.min(distortions)
    distortions /= np.sum(distortions)
    c_dist = np.cumsum(distortions)
    # Get the optimal k value (first above threshold)
    kmeans = KMeans(n_clusters=np.argwhere(c_dist > 0.90)[0][0], init='k-means++', max_iter=30000, n_init=1, verbose=False).fit(pc)
    labels = kmeans.labels_

    # Calculate SNR.
    snr = np.std(timecourse_matrix[:,1:31,(0,2)],axis=1) / np.std(timecourse_matrix[:,31:,(0,2)],axis=1)
    snr = np.max(snr,axis=1)

    # Calculate the percentage of spikes that violate refractoriness (<2ms isi)
    pct_refractory = np.sum(acf[:,:4], axis=1) * 100

    # May have to play with the cutoff
    idx_bad = np.argwhere((snr < snr_threshold) | (pct_refractory > refractory_threshold))[:,0]
    
    category_list = list()

    # Bad cells.
    count_threshold = np.median(spike_count)*0.25
    idx_bad_low = idx_bad[np.where(spike_count[idx_bad] <= count_threshold)[0]]
    idx_bad_high = idx_bad[np.where(spike_count[idx_bad] > count_threshold)[0]]
    for idx in cluster_id[idx_bad_high]:
        category_list.append(str(idx) + '  All/lowsnr')

    for idx in cluster_id[idx_bad_low]:
        category_list.append(str(idx) + '  All/crap')

    # Remove the bad clusters.
    for i, label in enumerate(labels):
        if not np.in1d(i, idx_bad)[0]:
            category_list.append(str(cluster_id[i]) + '  All/nc' + str(label))

    # Kick out the clusters to a classifications file.
    with open(output_path, 'w') as f:
        for line in category_list:
            f.write(line)
            f.write('\n')

def cluster_data(timecourse_matrix, isi, hull_area, significance_area, total_spikes, variance_or_components=0.85,
            refractory_threshold=0.2, snr_threshold=2.0, isi_binning=0.5, count_threshold=300, max_clusters=25):
    cell_clusters = np.zeros((timecourse_matrix.shape[0],), dtype=np.int32)
    # Find the bad and good cells
    
    # Calculate the percentage of spikes that violate refractoriness (<2ms isi)
    pct_refractory = compute_refractory_pct(isi, isi_binning)

    # Calculate SNR.
    snr = compute_snr(timecourse_matrix)

    # May have to play with the cutoff
    idx_bad = np.argwhere((pct_refractory > refractory_threshold) | (total_spikes < count_threshold))[:,0]
    idx_good = np.argwhere((snr >= snr_threshold) & (pct_refractory <= refractory_threshold) & (total_spikes >= count_threshold))[:,0]
    idx_low_snr = np.argwhere((snr < snr_threshold))[:,0]

    cell_clusters[idx_bad] = -2
    cell_clusters[idx_low_snr] = -1
    # Initialize your design matrix
    # pc = np.divide(hull_area[idx_good], np.nanstd(hull_area[idx_good])
    pc = hull_area[idx_good] / np.nanstd(hull_area[idx_good])
    
    if pc.ndim == 1:
        pc = pc[:,np.newaxis]
    # Get the ratio of hull area to significance area (should be different than one for spotty RFs).
    rf_ratio = significance_area[idx_good] / hull_area[idx_good]
    rf_ratio /= np.nanstd(rf_ratio)
    if rf_ratio.ndim == 1:
        rf_ratio = rf_ratio[:,np.newaxis]
    pc = np.append(pc, rf_ratio, axis=1)

    # Run SVD on the time courses and ISI distributions.
    pc_red = get_weighted_pca(timecourse_matrix[idx_good,:,0], variance_or_components=variance_or_components)
    pc_blue = get_weighted_pca(timecourse_matrix[idx_good,:,2], variance_or_components=variance_or_components)
    pc_isi = get_weighted_pca(isi[idx_good,:], variance_or_components=variance_or_components)
    pc = np.append(pc, pc_red, axis=1)
    pc = np.append(pc, pc_blue, axis=1)
    pc = np.append(pc, pc_isi, axis=1)
    # Use the elbow method to find the optimal number of clusters for KMeans.
    
    pc = np.nan_to_num(pc)
    num_clusters = get_num_clusters(pc)
    if num_clusters > max_clusters:
        num_clusters = max_clusters
    
    # Run KMeans
    kmeans = KMeans(n_clusters=num_clusters, init='k-means++', max_iter=3000, n_init=10, verbose=False).fit(pc)
    labels = kmeans.labels_
    cell_clusters[idx_good] = labels + 1
    return cell_clusters



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Convert Kilosort output to Vision .neurons format')
    parser.add_argument('filepath', type=str, help='folder containing Kilosort outputs, must already exist')
    parser.add_argument('output_path', type=str, help='folder containing Kilosort outputs, must already exist')
    
    args = parser.parse_args()

    print('Loading .mat')
    mdic = hdf5storage.loadmat(args.filepath)
    print('Done loading, clustering now')

    cluster_experiment(args.output_path, 
            mdic['sta'], 
            mdic['timecourse_matrix'], 
            mdic['acf'], 
            mdic['significance_maps'], 
            mdic['hull_area'],
            np.sum(mdic['spike_count'][:,(0,2)], axis=1), 
            mdic['gauss_params'], 
            mdic['cluster_id'], 
            refractory_threshold=0.1,
            snr_threshold=2.0)

# filepath = '/Users/michaelmanookin/Documents/Data/rawdata/MEA_Data/Analysis/20220712C/yass/20220712C_sta.mat'
# output_path = '/Users/michaelmanookin/Documents/Data/rawdata/MEA_Data/Analysis/20220712C/yass/noise/noise.classification.txt'
# # filepath = '/Users/michaelmanookin/Documents/Data/rawdata/MEA_Data/Analysis/20220712C/kilosort2/20220712C_sta.mat'
# # output_path = '/Users/michaelmanookin/Documents/Data/rawdata/MEA_Data/Analysis/20220712C/kilosort2/noise/noise.classification.txt'
# mdic = hdf5storage.loadmat(args.filepath)

# lobe_pts = range(4,13)
# refractory_threshold = 0.1 # No more than 0.1% of spikes can occur in first 2 msec of ISI.

# # Grab the properties to test.
# acf = mdic['acf']
# gauss_params = mdic['gauss_params']
# gabor_params = mdic['gabor_params']
# timecourse_matrix = mdic['timecourse_matrix']
# hull_area = mdic['hull_area']
# significance_maps = mdic['significance_maps']
# cluster_id = mdic['cluster_id']
# sta = mdic['sta']
# spike_count = np.sum(mdic['spike_count'][:,(0,2)], axis=1)
# spatial_maps = mdic['spatial_maps']
# # spike_count = np.sum(mdic['isi'],axis=1)

# # foo = np.reshape(timecourse_matrix[:,:,(0,2)],(timecourse_matrix.shape[0],timecourse_matrix.shape[1]*2))
# # my_model = PCA(n_components=20)


# # blue_pc = my_model.fit_transform(timecourse_matrix[:,:,2])
# # t_pc = my_model.fit_transform(np.reshape(timecourse_matrix[:,:,(0,2)],(timecourse_matrix.shape[0],timecourse_matrix.shape[1]*2)))

# # Weighted k-means!! Weight the components by their explained variance.
# my_model = PCA(0.9)
# t_pc = my_model.fit_transform(np.reshape(timecourse_matrix[:,:,(0,2)],(timecourse_matrix.shape[0],timecourse_matrix.shape[1]*2)))




# yellow_pc = PCA(n_components=2).fit_transform(timecourse_matrix[:,:,0])

# # Compute the STA over the lobe points and compare to the spatial maps to pick out 
# spatial_peaks = np.mean(sta[:,lobe_pts,:,:,:],axis=1)


# # Test SVD for space-time invariance.
# s = sta[171,:,:,:,0] # ST-variant
# s = sta[35,:,:,:,0] # ST-invariant

# nt, ny, nx = sta.shape[1:4]

# r2 = np.zeros((sta.shape[0],))
# for i in range(sta.shape[0]):
#     s = sta[i,:,:,:,0]
#     # Pass the STA through a Gaussian spatial filter.
#     s = gaussian_filter(s, sigma=0.75)
#     s_mat = s.reshape((nt, ny*nx))
#     space, singular_values, time = np.linalg.svd(s_mat.T, full_matrices=False)
#     r1 = np.outer(time[0,:], space[:,0]) * singular_values[0]
#     # Get the r^2 value
#     coeff = np.corrcoef(s_mat.ravel(), r1.ravel())
#     r2[i] = coeff[0,1]





# # s_mat = s.reshape((nt*ny, nx))
# s_mat = s.reshape((nt, ny*nx))
# # Subtract the mean.
# # s_mat -= s_mat.mean(axis=0, keepdims=True)

# # Temporal components
# t_pc = TruncatedSVD(n_components=5).fit(s_mat.T - s_mat.mean(axis=1, keepdims=True).T)
# # Spatial Components
# s_pc = TruncatedSVD(n_components=5).fit(s_mat - s_mat.mean(axis=0, keepdims=True))

# space, singular_values, time = np.linalg.svd(s_mat.T, full_matrices=False)

# r1 = np.outer(time[0,:], space[:,0]) * singular_values[0]

# coeff = np.corrcoef(s_mat.ravel(), r1.ravel())

# s_components = s_pc.components_

# c = t_pc.components_
# v = t_pc.explained_variance_ratio_
# v /= np.sum(v)

# filter_std = np.std(c, axis=1)



# # Calculate the S/T invariant and peak maps for the yellow channel.
# inv_map = np.mean(spatial_maps[:,:,:,(0,1)],axis=3)
# pk_map = np.mean(spatial_peaks[:,:,:,(0,1)],axis=3)

# # Normalize.
# n_inv = np.linalg.norm(inv_map,axis=(1,2))
# n_pk = np.linalg.norm(pk_map,axis=(1,2))
# for i in range(inv_map.shape[0]):
#     inv_map[i,...] /= n_inv[i]
#     pk_map[i,...] /= n_pk[i]

# # Take the dot product to get the correlation
# rf_corr = np.zeros((inv_map.shape[0],))
# for i in range(inv_map.shape[0]):
#     rf_corr[i] = np.sum(inv_map[i,...] * pk_map[i,...])

# # Possibly space-time invariant
# # rf_corr = 1 - np.abs(rf_corr)
# # idx_var = np.where(rf_corr > 0.9)[0]
# idx_var = np.where(rf_corr < 0.0)[0]

# # Compu

# # inv_map /= np.tile(n[:,np.newaxis,np.newaxis],(inv_map.shape[1],inv_map.shape[2]))
# # pk_map /= np.tile(n[:,np.newaxis,np.newaxis],(pk_map.shape[1],pk_map.shape[2]))

# # Get the area of significance in the maps.
# significance_area = np.sum(significance_maps, axis=(1,2))

# # Calculate the percentage of spikes that violate refractoriness (<2ms isi)
# pct_refractory = np.sum(acf[:,:1], axis=1) * 100

# # The first cut is based on significance. Non-significant should get cut first...
# # n_sig = np.zeros((sta.shape[0],))
# # for idx, cell_sta in enumerate(sta):
# #     sig_stixels = np.argwhere(lnp.calculate_sig_stixels(cell_sta, alpha=1/10000)[0])
# #     n_sig[idx] = sig_stixels.shape[0]

# # pca = PCA(n_components=2)
# # pc = pca.fit_transform(timecourse_matrix[:,:,0])
# # # Compute the loadings.
# # X_train_pca = pca.transform(X_train)

# # Calculate SNR.
# snr = np.std(timecourse_matrix[:,1:31,(0,2)],axis=1) / np.std(timecourse_matrix[:,31:,(0,2)],axis=1)
# snr = np.max(snr,axis=1)

# # May have to play with the cutoff
# # idx_bad = np.argwhere(n_sig < 1)
# idx_bad = np.argwhere((snr < 2.0) | (pct_refractory > refractory_threshold))[:,0]

# # Calculate the magnitude of yellow/blue temporal components.
# # yellow_t = np.median(timecourse_matrix[:,lobe_pts,0], axis=1)
# # blue_t = np.median(timecourse_matrix[:,lobe_pts,2], axis=1)
# # # Better way. 
# # for i in range(timecourse_matrix.shape[0]):
# #     yellow_t[i] = find_timecourse_peak(timecourse_matrix[i,:,0])
# #     blue_t[i] = find_timecourse_peak(timecourse_matrix[i,:,2])

# # # Normalize the amplitudes.
# # mval = np.max(np.abs(np.append(yellow_t,blue_t)))
# # yellow_t /= mval
# # blue_t /= mval

# # # Try the initial sorting by contrast polarity.
# # idx_blue = np.where((blue_t > 0.2) & (yellow_t < -0.2))[0]
# # idx_off = np.where(yellow_t < 0.0)[0]
# # idx_on = np.where(yellow_t >= 0.0)[0]
# # idx_off = np.where(yellow_t < 0.0 & ((blue_t <= 0.1) & (yellow_t >= -0.1)))

# idx_blue, idx_off, idx_on = initial_polarity_split(timecourse_matrix)

# # Remove any potentially invariant cells from the other arrays.
# bool_array = np.where(np.in1d(idx_var, idx_bad))[0]
# idx_var = np.delete(idx_var, bool_array)
# bool_array = np.where(np.in1d(idx_blue, idx_var))[0]
# idx_blue = np.delete(idx_blue, bool_array)
# bool_array = np.where(np.in1d(idx_off, idx_var))[0]
# idx_off = np.delete(idx_off, bool_array)
# bool_array = np.where(np.in1d(idx_on, idx_var))[0]
# idx_on = np.delete(idx_on, bool_array)

# # Remove any blue cells from the other arrays.
# bool_array = np.where(np.in1d(idx_off, idx_blue))[0]
# idx_off = np.delete(idx_off, bool_array)
# bool_array = np.where(np.in1d(idx_on, idx_blue))[0]
# idx_on = np.delete(idx_on, bool_array)


# # Remove the bad cells.
# bool_array = np.where(np.in1d(idx_blue, idx_bad))[0]
# idx_blue = np.delete(idx_blue, bool_array)
# bool_array = np.where(np.in1d(idx_off, idx_bad))[0]
# idx_off = np.delete(idx_off, bool_array)
# bool_array = np.where(np.in1d(idx_on, idx_bad))[0]
# idx_on = np.delete(idx_on, bool_array)

# # Create a category list.
# category_list = list()

# blue_labels = cluster_cells(hull_area[idx_blue],gauss_params[idx_blue,...], significance_area[idx_blue], timecourse_matrix[idx_blue,...], acf[idx_blue,...],2)
# for i, label in enumerate(blue_labels):
#     idx = idx_blue[i]
#     category_list.append(str(cluster_id[idx]) + '  All/Blue/nc' + str(label))

# cat_offset = np.max(blue_labels) + 1

# # Off cells
# off_labels = cluster_cells(hull_area[idx_off],gauss_params[idx_off,...], significance_area[idx_off], timecourse_matrix[idx_off,...], acf[idx_off,...],10)
# for i, label in enumerate(off_labels):
#     idx = idx_off[i]
#     category_list.append(str(cluster_id[idx]) + '  All/Off/nc' + str(label + cat_offset))

# cat_offset += np.max(off_labels) + 1

# on_labels = cluster_cells(hull_area[idx_on],gauss_params[idx_on,...], significance_area[idx_on], timecourse_matrix[idx_on,...], acf[idx_on,...],10)
# for i, label in enumerate(on_labels):
#     idx = idx_on[i]
#     category_list.append(str(cluster_id[idx]) + '  All/On/nc' + str(label + cat_offset))

# for idx in idx_var:
#     category_list.append(str(cluster_id[idx]) + '  All/stvar')

# # Bad cells.
# count_threshold = 20000
# # idx_bad_low = cluster_id[np.where(spike_count[idx_bad] <= count_threshold)[0]]
# # idx_bad_high = cluster_id[np.where(spike_count[idx_bad] > count_threshold)[0]]
# idx_bad_low = idx_bad[np.where(spike_count[idx_bad] <= count_threshold)[0]]
# idx_bad_high = idx_bad[np.where(spike_count[idx_bad] > count_threshold)[0]]
# for idx in cluster_id[idx_bad_high]:
#     category_list.append(str(idx) + '  All/lowsnr')

# for idx in cluster_id[idx_bad_low]:
#     category_list.append(str(idx) + '  All/crap')

# # Kick out the clusters to a classifications file.
# with open(output_path, 'w') as f:
#     for line in category_list:
#         f.write(line)
#         f.write('\n')



# model = Pipeline([('svd', TruncatedSVD(n_components=12)), ('m', LogisticRegression())])

# # get the models to evaluate
# models = get_models()
# # evaluate the models and store results
# results, names = list(), list()
# for name, model in models.items():
# 	scores = evaluate_model(model, X, y)
# 	results.append(scores)
# 	names.append(name)
# 	print('>%s %.3f (%.3f)' % (name, mean(scores), std(scores)))

