"""
Functions related to fitting LNP model to data.
"""
import sys, os
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import numpy.ma as ma
import scipy as sp
from scipy import stats
# import config as cfg
import pdb
from tqdm import tqdm
# from progressbar import *
from scipy.linalg import blas
from scipy import linalg
import torch
from scipy.spatial import ConvexHull
import cv2
import gc
from scipy.optimize import curve_fit
import analysis_utils
import warnings
from typing import Tuple
from sklearn.decomposition import PCA
from scipy import ndimage
from scipy.signal import savgol_filter

from scipy.ndimage import gaussian_filter
# from analysis_utils import Analysis

# Set number of threads for torch.
# torch.set_num_threads(cfg.TORCH_N_THREADS)

#tfun = @(p,t)(p(1) .* (((t./p(2)).^p(6))./(1+((t./p(2)).^p(6)))) .* exp(-((t./p(3)))) .* cos(((2.*pi.*t)./p(4))+(2*pi*p(5)/360)));
# tparams = [62.562745558862289   0.061433021205849   0.016320868351090   0.560874309169834  52.432435563319949 4.112884082863220];
def linear_filter_function(t, *params):
    """ Linear filter function. 
    Arguments:
        t (np.ndarray): Time vector.
        params (list): List of parameters.
    Returns:
        f (np.ndarray): Filter output.
    """
    f = params[0] * (((t/params[1]) ** params[5]) / 
        (1+((t/params[1])**params[5]))) * np.exp(-((t/params[2]))) * np.cos(((2*np.pi*t) / params[3])+(2*np.pi*params[4]/360))
    return f

# Fit the linear filter functions
def fit_linear_filters(time_course, bin_rate=120.0):
    t = np.arange(0, time_course.shape[0]) / bin_rate
    time_fits = np.zeros((time_course.shape))
    # Define the boundaries.
    non_bds = ([-1000, 0.01, 0.005, 0.01, -360, 4.0], 
           [1000, 0.3, 0.5, 1.0, 360.0, 10.0])
    t_params = np.zeros((time_course.shape[1], 6))
    for j in range(0,time_course.shape[1]):
        # Quick check of the filter sign.
        lobe_pts = np.sum(time_course[np.round(0.03*bin_rate).astype(int)-1 : np.round(0.1*bin_rate).astype(int),j])
        if lobe_pts > 0.0:
            sgn = 1
            filter_peak = np.argmax(time_course[:,j], axis=0)
        else:
            sgn = -1
            filter_peak = np.argmin(time_course[:,j], axis=0)
        # Initial parameter estimates.
        p0 = [sgn*7, filter_peak/bin_rate, 0.06, 0.6, 90.0, 5.0]
        # if lobe_pts > 0.0:
        #     p0 = [62.562745558862289,  0.061433021205849,  0.016320868351090,  0.560874309169834, 52.432435563319949, 4.112884082863220]
        # else:
        #     p0 = [-62.562745558862289,  0.061433021205849,  0.016320868351090,  0.560874309169834, 52.432435563319949, 4.112884082863220]
        popt, _ = curve_fit(f=linear_filter_function, xdata=t, ydata=time_course[:,j].T, p0=p0, bounds=non_bds, max_nfev=1e05)
        time_fits[:,j] = linear_filter_function(t, *popt).T
        t_params[j,:]  = popt
    return time_fits, t_params

def zero_sta_low_freq(cell_sta: np.ndarray, threshold_frequency: float=2.0, bin_rate: float=120.0, high_threshold_frequency: float=None) -> np.ndarray:
    """ Zero out the low spatial frequency components of the STA.
    Parameters:
        cell_sta: The STA of the cell (size: t,y,x,c).
        low_freq: The low frequency cutoff.
    Returns:
        sta: The STA with the low frequency components zeroed out.
    """
    # Compute the FFT.
    fft_sta = np.fft.fft(cell_sta, axis=0)
    # Get the corresponding frequencies
    frequencies = np.fft.fftfreq(cell_sta.shape[0], d=1.0/bin_rate)
    # Create a mask to zero out low frequencies
    if high_threshold_frequency is None:
        mask = np.abs(frequencies) > threshold_frequency
    else:
        mask = np.logical_and(np.abs(frequencies) > threshold_frequency, np.abs(frequencies) < high_threshold_frequency)
    # Apply the mask.
    fft_sta *=  mask[:,np.newaxis,np.newaxis,np.newaxis]
    # Reconstruct the signal using the inverse FFT
    filtered_sta = np.fft.ifft(fft_sta, axis=0)
    return filtered_sta.real

def calculate_map_and_time_course(cell_sta: np.ndarray, threshold: float=3.0):
    sta_depth,_,_,n_phosphors = cell_sta.shape
    timecourse_matrix = np.zeros((sta_depth,n_phosphors))
    sig_stixels = calculate_sig_stixels_variance(cell_sta, threshold=threshold)
    if sig_stixels.shape[0] == 0:
        return timecourse_matrix, None
    for stixel in sig_stixels:
        timecourse_matrix += cell_sta[:,stixel[0],stixel[1],:]
    # Divide by the max L2 norm.
    l2 = np.max(np.linalg.norm(timecourse_matrix, axis=0))
    timecourse_matrix /= l2
    return timecourse_matrix, sig_stixels

def compute_spatial_map(cell_sta: np.ndarray, timecourse_matrix: np.ndarray) -> np.ndarray:
    """
    Compute the spatial RF map for a cell.

    Parameters: 
        cell_sta: STA of a cell (size: t,y,x,c) with zero mean
        timecourse_matrix: Time course matrix for the cell (size: t,c)

    Returns:
        spatial_map: Spatial RF map (size: y,x)
    """
    _,x,y,n_phosphors = cell_sta.shape
    spatial_map = np.zeros((x,y,n_phosphors))
    for j in range(timecourse_matrix.shape[1]):
        spatial_map[:,:,j] = np.dot(cell_sta[:,:,:,j].T, timecourse_matrix[:,j]).T
    return spatial_map

def compute_hull(significance_map: np.ndarray):
    """
    Fit a convex hull to the significance map.

    Parameters:
        significance_map: Map of significant stixels (y,x)
    
    Returns: Tuple containing hull vertices and hull area
    """
    x, y = significance_map.shape
    scale = np.floor(300/x).astype(int)
    if scale > 1:
        res = cv2.resize(significance_map, dsize=(y*scale, x*scale), interpolation=cv2.INTER_CUBIC)
    else:
        scale = 1
        res = significance_map
    res /= np.max(res)
    # Get the significant points.
    points = np.argwhere(res > 0.1)
    # Fit a convex hull to the significant stixels.
    if points.shape[0] < 4:
        return np.zeros((1,2)), 0.0
    hull = ConvexHull(points)
    hull_area = hull.volume / np.power(scale, 2)
    hull_vertices = points[hull.vertices,:] / scale
    hull_vertices = hull_vertices[:,::-1]
    return hull_vertices, hull_area

def in_hull(p, hull):
    """
    Test if points in `p` are in `hull`

    `p` should be a `NxK` coordinates of `N` points in `K` dimensions
    `hull` is either a scipy.spatial.Delaunay object or the `MxK` array of the 
    coordinates of `M` points in `K`dimensions for which Delaunay triangulation
    will be computed
    """
    if hull.shape[0] < 4:
        return None
    
    from scipy.spatial import Delaunay
    if not isinstance(hull,Delaunay):
        hull = Delaunay(hull)

    return hull.find_simplex(p)>=0

def map_from_points(points: np.ndarray, x_size: int, y_size: int) -> np.ndarray:
    """
    Create a significane map from array of significant stixels

    Parameters:
        points: x/y locations of significant stixels (size: n,2)
        x_size: number of stixels along the x-axis
        y_size: number of stixels along the y-axis
    
    Returns: significance map (size: x,y)
    """
    my_map = np.zeros((x_size, y_size))
    if points.shape[0] > 0:
        for stixel in points:
            my_map[stixel[0],stixel[1]] = 1
    return my_map

def compute_mean_in_hull(h_vertices, space_map: np.ndarray):
    """
    Compute the average intensity inside of a convex hull.

    Parameters:
        h_vertices: Vertices of the convex hull
        space_map: Spatial RF map (y,x)
    
    Returns: mean pixel value
    """
    try:
        x_size, y_size = space_map.shape
        X, Y = np.meshgrid(np.linspace(1, x_size, x_size), 
                np.linspace(1, y_size, y_size), indexing="ij")
        xy = np.zeros((len(X.flatten()),2))
        xy[:,0] = X.flatten()
        xy[:,1] = Y.flatten()
        foo = in_hull(xy, h_vertices)
        pts = np.argwhere(foo)[:,0]
        if len(pts) < 3:
            return 0.0
        hull_points = xy[pts,::-1] - 1
        hull_points = hull_points.astype(int)
        my_mean = 0.0
        if hull_points.shape[0] > 0:
            for pt in hull_points:
                my_mean += space_map[pt[0], pt[1]]
            my_mean /= float(hull_points.shape[0])
        return my_mean
    except Exception as e:
        return 0.0

def compute_timecourse_in_hull(cell_sta: np.ndarray, h_vertices: np.ndarray) -> np.ndarray:
    """
    Compute the significant stixels based on variance of the time course.

    Parameters:
        cell_sta: STA of a cell size t,y,x,c with zero mean
        h_vertices: Vertices of the convex hull
    
    Returns: timecourse_matrix (t,c)
    """
    sta_depth,x_size,y_size,n_phosphors = cell_sta.shape
    timecourse_matrix = np.zeros((sta_depth,n_phosphors))
    x_size, y_size = cell_sta.shape[1:3]
    X, Y = np.meshgrid(np.linspace(1, x_size, x_size), 
            np.linspace(1, y_size, y_size), indexing="ij")
    xy = np.zeros((len(X.flatten()),2))
    xy[:,0] = X.flatten()
    xy[:,1] = Y.flatten()
    foo = in_hull(xy, h_vertices)
    hull_points = xy[np.argwhere(foo)[:,0],::-1] - 1
    hull_points = hull_points.astype(int)
    if hull_points.shape[0] == 0:
        return None
    for stixel in hull_points:
        timecourse_matrix += cell_sta[:,stixel[0],stixel[1],:]
    # Divide by the max L2 norm.
    l2 = np.max(np.linalg.norm(timecourse_matrix, axis=0))
    timecourse_matrix /= l2
    return timecourse_matrix

def compute_timecourse_from_index(cell_sta: np.ndarray, sig_stixels: np.ndarray) -> np.ndarray:
    """
    Compute the significant stixels based on variance of the time course.

    Parameters:
        cell_sta: STA of a cell size t,y,x,c with zero mean
        sig_stixels: x/y indices of statistically signficant stixels
    
    Returns: timecourse_matrix (t,c)
    """
    sta_depth,_,_,n_phosphors = cell_sta.shape
    timecourse_matrix = np.zeros((sta_depth,n_phosphors))
    if sig_stixels.shape[0] == 0:
        return timecourse_matrix
    for stixel in sig_stixels:
        timecourse_matrix += cell_sta[:,stixel[0],stixel[1],:]
    # Divide by the max L2 norm.
    l2 = np.max(np.linalg.norm(timecourse_matrix, axis=0))
    timecourse_matrix /= l2
    return timecourse_matrix

def find_timecourse_sign(time_course: np.ndarray) -> float:
    """
    Find the sign of the time course based on the polarity of the filter.
    Parameters:
        time_course: time course of a cell (size: t)
    Returns: filter polarity
    """
    extrema = np.abs([np.min(time_course), np.max(time_course)])
    positions = [np.argmin(time_course), np.argmax(time_course)]
    if positions[0] < positions[1]:
        filter_polarity = -1
    else:
        filter_polarity = 1
    if extrema[1] > 2 * extrema[0]:
        filter_polarity = 1
    elif extrema[0] > 2 * extrema[1]:
        filter_polarity = -1
    if filter_polarity > 0:
        peak_position = positions[1]
    else:
        peak_position = positions[0]
    return filter_polarity, peak_position

def filter_sta(cell_sta: np.ndarray, filter_sigma: float=0.75) -> np.ndarray:
    """
    Filter the STA with a Gaussian filter.

    Parameters:
        cell_sta: STA of a cell size t,y,x,c with zero mean
        filter_sigma: Sigma of the Gaussian filter.

    Returns:
        filtered_sta: Filtered STA
    """
    sta_depth,sta_height,sta_width,n_phosphors = cell_sta.shape
    # Filter the STA.
    filtered_sta = np.zeros((sta_depth,sta_height,sta_width,n_phosphors))
    for phosphor in range(n_phosphors):
        filtered_sta[:,:,:,phosphor] = gaussian_filter(cell_sta[:,:,:,phosphor], sigma=filter_sigma)
    return filtered_sta

def calculate_time_course(cell_sta: np.ndarray, threshold: float=3.0, max_stixels: int=3) -> np.ndarray:
    """
    Compute the significant stixels based on variance of the time course.

    Parameters:
        cell_sta: STA of a cell size t,y,x,c with zero mean
        threshold: Threshold for significant stixels in stdev above the mean.

    Returns:
        time_course: Computed time course (t,c)
    """
    sta_depth,_,_,n_phosphors = cell_sta.shape
    time_course = np.zeros((sta_depth,n_phosphors))
    # Deal with the first half of the STA only.
    half_pts = np.ceil(sta_depth/2).astype(int)
    # Pass the STA through a Gaussian spatial filter.
    f_sta = filter_sta(cell_sta[1:half_pts,...])
    # Find the significant stixels.
    sig_stixels = calculate_sig_stixels_variance(f_sta, threshold=threshold)
    if sig_stixels.shape[0] == 0:
        return time_course
    t_mat = np.zeros((sig_stixels.shape[0], f_sta.shape[0], n_phosphors))
    on_idx = list()
    off_idx = list()
    for index, stixel in enumerate(sig_stixels):
        t_mat[index,:,:] = f_sta[:,stixel[0],stixel[1],:]
        # Separately gather the On and Off time courses for the green channel.
        sign_mean,_ = find_timecourse_sign(t_mat[index,:,1])
        if sign_mean > 0:
            on_idx.append(index)
        else:
            off_idx.append(index)
    # Compute the mean to get the On and Off filters.
    on_filter = np.mean(t_mat[on_idx,:],axis=0)
    off_filter = np.mean(t_mat[off_idx,:],axis=0)
    if np.max(on_filter) > -np.min(off_filter):
        proj = np.dot(t_mat[:,:,1], on_filter[:,1])
    else:
        proj = np.dot(t_mat[:,:,1], off_filter[:,1])
    n_vals = np.min([len(proj), max_stixels])
    ind = np.argpartition(proj, -n_vals)[-n_vals:]
    sig_stixels = sig_stixels[ind,:]
    for stixel in sig_stixels:
        time_course += cell_sta[:,stixel[0],stixel[1],:]
    # Divide by the max L2 norm.
    l2 = np.max(np.linalg.norm(time_course, axis=0))
    time_course /= l2
    return time_course

def calculate_time_course_svd(cell_sta: np.ndarray, threshold: float=2.0, max_stixels: int=3) -> np.ndarray:
    """
    Compute the significant stixels based on variance of the time course.

    Parameters:
        cell_sta: STA of a cell size t,y,x,c with zero mean
        threshold: Threshold for significant stixels in stdev above the mean.

    Returns:
        time_course: Computed time course (t,c)
    """
    try:
        sta_depth,_,_,n_phosphors = cell_sta.shape
        time_course = np.zeros((sta_depth,n_phosphors))
        sig_stixels = calculate_sig_stixels_variance(cell_sta, threshold=threshold)
        if sig_stixels.shape[0] == 0:
            return time_course
        elif sig_stixels.shape[0] <= max_stixels:
            for stixel in sig_stixels:
                time_course += cell_sta[:,stixel[0],stixel[1],:]
        else:
            half_pts = np.ceil(sta_depth/2).astype(int)
            st_mat = np.zeros((sig_stixels.shape[0], cell_sta.shape[0], n_phosphors))
            for index, stixel in enumerate(sig_stixels):
                st_mat[index,:,:] = cell_sta[:,stixel[0],stixel[1],:]
            
            t_mean = np.mean(st_mat[:,:half_pts,1], axis=0)
            sign_mean,_ = find_timecourse_sign(t_mean)
            t_comp = PCA(n_components=1).fit(st_mat[:,:half_pts,1]).components_[0]
            sign_t,_ = find_timecourse_sign(t_comp)
            if sign_mean != sign_t:
                t_comp *= -1.0
            proj = np.dot(st_mat[:,:half_pts,index], t_comp)
            n_vals = np.min([len(proj), max_stixels])
            ind = np.argpartition(proj, -n_vals)[-n_vals:]
            time_course = np.mean(st_mat[ind,:,:], axis=0)
    except:
        return calculate_time_course(cell_sta, threshold=threshold)
    
    # Divide by the max L2 norm.
    l2 = np.max(np.linalg.norm(time_course, axis=0))
    time_course /= l2
    return time_course

def calculate_time_course_count(cell_sta: np.ndarray, max_stixels: int=3) -> np.ndarray:
    """
    Compute the significant stixels based on variance of the time course.

    Parameters:
        cell_sta: STA of a cell size t,y,x,c with zero mean
        threshold: Threshold for significant stixels in stdev above the mean.

    Returns:
        time_course: Computed time course (t,c)
    """
    sta_depth,_,_,n_phosphors = cell_sta.shape
    time_course = np.zeros((sta_depth,n_phosphors))
    # Deal with the first half of the STA only.
    half_pts = np.ceil(sta_depth/2).astype(int)
    f_sta = cell_sta[1:half_pts,...].copy()
    # Find the significant stixels.
    sig_stixels = calculate_sig_stixels_count(f_sta, max_stixels=max_stixels)
    if sig_stixels.shape[0] == 0:
        return time_course
    for stixel in sig_stixels:
        time_course += cell_sta[:,stixel[0],stixel[1],:]
    # Divide by the max L2 norm.
    l2 = np.max(np.linalg.norm(time_course, axis=0))
    time_course /= l2
    return time_course

def threshold_space_map(space_map: np.ndarray, filter_sigma: float=1.0, threshold: float=2.0) -> np.ndarray:
    """
    Compute the significant stixels based on variance of the time course.

    Parameters:
        cell_sta: STA of a cell size t,y,x,c with zero mean
        threshold: Threshold for determining significant stixels in standard deviations above the mean.

    Returns: numpy array denoting the x/y values for the significant stixels.
    """
    space_map = gaussian_filter(space_map, sigma=filter_sigma)
    sig_stixels = np.argwhere((space_map / np.std(space_map)) >= threshold)
    return sig_stixels

def compute_filtered_map(cell_sta: np.ndarray, threshold: float=2.0) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute the significant stixels based on variance of the time course.

    Parameters:
        cell_sta: STA of a cell size t,y,x,c with zero mean
        threshold: Threshold for determining significant stixels in standard deviations above the mean.

    Returns: numpy array denoting the x/y values for the significant stixels.
    """
    _, x, y, _ = cell_sta.shape
    sig_stixels = calculate_sig_stixels_variance(cell_sta, threshold=threshold)
    # Make the significance map.
    significance_map = np.zeros((x,y))
    if sig_stixels.shape[0] > 0:
        tc = compute_timecourse_from_index(cell_sta, sig_stixels)
        space_map = compute_spatial_map(cell_sta, tc)
        sig_stixels = np.argwhere(space_map / np.std(space_map) >= threshold)
    if sig_stixels.shape[0] > 0:
        for stixel in sig_stixels:
            significance_map[stixel[0],stixel[1]] = 1  
    return significance_map, sig_stixels

def threshold_space_map(space_map: np.ndarray, filter_sigma: float=1.0, threshold: float=2.0) -> np.ndarray:
    """
    Compute the significant stixels based on variance of the time course.

    Parameters:
        cell_sta: STA of a cell size t,y,x,c with zero mean
        threshold: Threshold for determining significant stixels in standard deviations above the mean.

    Returns: numpy array denoting the x/y values for the significant stixels.
    """
    space_map = gaussian_filter(space_map, sigma=filter_sigma)
    sig_stixels = np.argwhere((space_map / np.std(space_map)) >= threshold)
    return sig_stixels

def find_best_map_threshold(cell_sta: np.ndarray, threshold: float=2.0, filter_sigma: float=1.0) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Find the best threshold and signficant stixels for an RF map.

    Parameters:
        cell_sta: STA of a cell size t,y,x,c with zero mean
        threshold: Threshold for determining significant stixels in standard deviations above the mean.
        filter_sigma: standard deviation for the Gaussian spatial filter (in stixels).

    Returns: Tuple containing the significant stixels, best threshold and the mean intensity within the hull for each threshold.
    """
    _, x, y, _ = cell_sta.shape
    sig_stixels = calculate_sig_stixels_variance(cell_sta, threshold=threshold)
    space_map = compute_spatial_map(cell_sta, compute_timecourse_from_index(cell_sta, sig_stixels))
    space_map = np.mean(space_map, axis=2)
    space_map -= np.mean(space_map)
    thresholds = np.arange(2.0, 6.0, 0.5)
    mean_intensity = np.zeros((len(thresholds),))
    for thresh_count, thresh in enumerate(thresholds):
        sig_stixels = threshold_space_map(space_map, filter_sigma=filter_sigma, threshold=thresh)
        if sig_stixels.shape[0] > 0:
            h_vert, _ = compute_hull(map_from_points(sig_stixels,x,y))
            mean_intensity[thresh_count] = compute_mean_in_hull(h_vert, space_map)
    if np.all(mean_intensity <= 0.0):
        best_threshold = threshold
        sig_stixels = threshold_space_map(space_map, filter_sigma=filter_sigma, threshold=best_threshold)
    else:
        # Compute the derivative of the mean values.
        mean_intensity /= np.max(mean_intensity)
        diff_mean = np.diff(mean_intensity)
        max_idx = np.argwhere(diff_mean == np.max(diff_mean))[0][0]+1
        if max_idx is not None:
            best_threshold = thresholds[max_idx]
            sig_stixels = threshold_space_map(space_map, filter_sigma=filter_sigma, threshold=best_threshold)
        else:
            best_threshold = threshold
            sig_stixels = threshold_space_map(space_map, filter_sigma=filter_sigma, threshold=best_threshold)
    return sig_stixels, best_threshold, mean_intensity

def compute_significance_map(cell_sta: np.ndarray, threshold: float=2.0, filter_sigma: float=1.0) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute the significant stixels based on variance of the time course.

    Parameters:
        cell_sta: STA of a cell size t,y,x,c with zero mean
        threshold: Threshold for determining significant stixels in standard deviations above the mean.
        filter_sigma: standard deviation for the Gaussian spatial filter (in stixels).

    Returns: Tuple containing the significance map and a numpy array denoting the x/y values for the significant stixels.
    """
    _, x, y, _ = cell_sta.shape
    try:
        sig_stixels, _, _ = find_best_map_threshold(cell_sta, threshold=threshold, filter_sigma=filter_sigma)
    except Exception as e:
        print(e)
        sig_stixels = calculate_sig_stixels_variance(cell_sta, threshold=threshold)
    # Make the significance map.
    significance_map = np.zeros((x,y))
    if sig_stixels.shape[0] > 0:
        significance_map = map_from_points(sig_stixels, x, y)
    return significance_map, sig_stixels

def calculate_sig_stixels_count(cell_sta: np.ndarray, max_stixels: int=3) -> np.ndarray:
    """
    Compute the significant stixels based on variance of the time course.

    Parameters:
        cell_sta: STA of a cell size t,y,x,c with zero mean
        threshold: Threshold for determining significant stixels in standard deviations above the mean.

    Returns: numpy array denoting the x/y values for the significant stixels.
    """
    # Compare the variance/stdev in the first vs second half of the STA
    if cell_sta.shape[0] == 3:
        v_base = np.var(cell_sta[(0,2),...], axis=0)
        v = np.var(cell_sta[(0,1),...], axis=0)
    else:
        half_pts = np.ceil(cell_sta.shape[0]/2).astype(int)
        # Compute the baseline variance in the second half of the sta.
        v_base = np.var(cell_sta[half_pts:,...], axis=0)
        v_base = np.mean(v_base, axis=2)
        # Compute the variance along the time axis.
        v = np.var(cell_sta[:half_pts,...], axis=0)
    # Average across the color channels.
    v = np.mean(v, axis=2)
    # Get the top n stixels.
    flat_index = np.argpartition(v.ravel(), -max_stixels)[-max_stixels:]
    ind = np.unravel_index(flat_index, v.shape)
    return np.array(ind).T

def calculate_sig_stixels_variance(cell_sta: np.ndarray, threshold: float=3.0) -> np.ndarray:
    """
    Compute the significant stixels based on variance of the time course.

    Parameters:
        cell_sta: STA of a cell size t,y,x,c with zero mean
        threshold: Threshold for determining significant stixels in standard deviations above the mean.

    Returns: numpy array denoting the x/y values for the significant stixels.
    """
    # Compare the variance/stdev in the first vs second half of the STA
    if cell_sta.shape[0] == 3:
        v_base = np.var(cell_sta[(0,2),...], axis=0)
        v = np.var(cell_sta[(0,1),...], axis=0)
    else:
        half_pts = np.ceil(cell_sta.shape[0]/2).astype(int)
        # Compute the baseline variance in the second half of the sta.
        v_base = np.var(cell_sta[half_pts:,...], axis=0)
        v_base = np.mean(v_base, axis=2)
        # Compute the variance along the time axis.
        v = np.var(cell_sta[:half_pts,...], axis=0)
    # Average across the color channels.
    v = np.mean(v, axis=2)
    v -= np.mean(v_base)
    # Normalize by the standard deviation
    v /= np.std(v_base)
    return np.argwhere(v >= threshold)

def calculate_sig_stixels(sta, num_frames=9, alpha=1/10000, color_channels=None):
    """
    Sam Cooler sig stix implementation.

    Parameters:
        sta: STA of size t,y,x,c, with zero mean!
        num_frames: frames to use in noise calculation
        alpha: significance level
        color_channels: color channels to use
    Returns: tuple of booleans denoting the sig stix and the p values.
    
    Flips the STA and re-orders dims, but changes it back (doesn't mutate).

    """
    # Rearrange the dims and flip.
    sta = np.flip(np.moveaxis(sta,0,-2),axis=-2)

    if color_channels is None:
        color_channels = (0,1,2)

    early_frames = np.mean(sta[...,0:num_frames,color_channels], 3) 
    std_this_cell = stats.median_abs_deviation(early_frames.flatten())

    # Mean over color, square, sum over time
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        input_frame = np.sum(np.power(np.mean(sta[...,-(num_frames+1):-1,
                                color_channels],3) / std_this_cell, 2), 2)

    # Chi squares survival = 1 - CDF
    pvalues = stats.chi2(df=num_frames).sf(input_frame) + np.finfo(float).eps 
    pvalues *= (sta.shape[0] * sta.shape[1]) # bonferroni correction

    # Put it back in order.
    sta = np.flip(np.moveaxis(sta,-2,0),axis=0)

    return pvalues < alpha, pvalues

def mask_sta(sta: np.ndarray, sig_stixels: np.ndarray) -> np.ndarray:
    """
    Arguments:
        sta: tensor of size t,y,x,c
        sig_stixels: array of stixel indices (0-indexed)
    Returns: STA with non-significant stixels set to 0.

    Sets all non significant stixels in an STA to zero.
    """

    masked_sta = np.zeros(sta.shape)

    for i in range(sig_stixels.shape[0]):
        y_stix = sig_stixels[i][1]
        x_stix = sig_stixels[i][0]
        masked_sta[:,x_stix,y_stix,:] = sta[:,x_stix,y_stix,:]

    return masked_sta

def get_frame_times(ttl_times,frames_per_ttl):
    """
    Linearly interpolates frame times based on the ttls

    Parameters:
        ttl_times: ttl times
        frames_per_ttl: frames per ttl
    
    Returns:
        Array of interpolated frame times.
    """

    frame_times = []

    for i in range(ttl_times.shape[0]-1):
        frame_times.append(np.linspace(ttl_times[i],
                               ttl_times[i+1],
                               frames_per_ttl,
                               endpoint=False))

    return np.asarray(frame_times).ravel()

def get_binned_spikes(vcd,cells):
    '''
    Parameters:
        vcd: Vision data table object
        cells: list of cells to compute STAs for
    Returns:
        A matrix of size frames by cells

    Bins the spike train using the (full resolution) frames of the monitor 
    as bin edges.
    '''
    frames_per_ttl = vcd.runtimemovie_params.framesPerTTL

    # Get the interpolated frame times.
    ttl_times = vcd.ttl_times / cfg.FS * 1000 # ms.
    frame_times = get_frame_times(ttl_times,frames_per_ttl)

    # Set the bin edges as as the frame times with an additional fake frame.
    monitor_interval = 1 / vcd.runtimemovie_params.monitorFrequency * 1000
    bin_edges = np.append(frame_times,frame_times[-1] +  monitor_interval)

    # Loop through cells and bin spike train.
    binned_spikes = []

    for cell in cells:
        spike_times = vcd.get_spike_times_for_cell(cell) / cfg.FS * 1000 # ms
        binned_spikes.append(np.histogram(spike_times,bins=bin_edges)[0])
    
    return np.asarray(binned_spikes).T

def normalize_sta(sta_tensor):
    """
    Normalizes STA according to some constants. Puts it in 0-1, with mean .5.
    Makes a copy, does not mutate.
    Parameters:
        sta_tensor: STA tensor of size cells, t, y, x, color
    Returns:
        normalized STA tensor of same size as input
    """
    sta_tensor_norm = sta_tensor.copy()
    sta_tensor_norm /= cfg.MONITOR_BIT_DEPTH
    sta_tensor_norm -= cfg.STA_NORM_MEAN
    sta_tensor_norm /= np.amax(np.abs(sta_tensor_norm),
                        axis=(1,2,3,4)).reshape((-1,1,1,1,1))
    sta_tensor_norm /= 2
    sta_tensor_norm += cfg.STA_NORM_MEAN

    return sta_tensor_norm

def compute_sta_numpy(stimulus,binned_spikes,stride,depth=30,norm=True):
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

    # Intialize STA tensor.
    num_cells = binned_spikes.shape[1]
    n_frames,height,width,n_phosphors = stimulus.shape
    sta_tensor = np.zeros((num_cells,depth,height,width,n_phosphors))

    # Vectorize stimulus.
    stimulus = np.reshape(stimulus,(n_frames,
                            height * width * n_phosphors)).astype(np.float32).T
    binned_spikes = binned_spikes.astype(np.float32)

    # Initialize STA, the offset range, and compute.
    offset_cnt = 0
    sta_tensor = np.zeros((num_cells,depth,height,width,n_phosphors))
    offset_range = np.arange(0,(depth * stride),stride)

    # Initialize progress bar.
    widgets = ['Computing STAs: ', Percentage(), ' ', Bar(marker='*',
               left='[',right=']'),' ', ETA()]
    pbar = ProgressBar(widgets=widgets, maxval=depth)
    pbar.start()

    for offset in offset_range:
        sta = np.matmul(stimulus[:,:n_frames-offset],
                        binned_spikes[offset:,:]
                       )

        # Index appropriately to avoid nans and normalize by number of spikes.
        nonzero_inds = np.nonzero(np.sum(binned_spikes[offset:,:],0))[0]
        sta[:,nonzero_inds] /= np.sum(binned_spikes[offset:,:],0)[nonzero_inds]
        sta_tensor[:,offset_cnt,:,:,:] = np.reshape(sta.T,
                                                  (num_cells,height,
                                                   width,n_phosphors))

        del sta

        pbar.update(offset_cnt)
        offset_cnt +=1

    pbar.finish()

    if not norm:
        return sta_tensor

    return normalize_sta(sta_tensor)

def compute_sta_blas(stimulus,binned_spikes,stride,depth=30,norm=True):
    """
    Computes STA by convolving the spiketrain with the stimulus. This 
    implementation relies on the BLAS interfaces that numpy wraps around.

    Parameters:
        stimulus: stimulus tensor of size t, y, x, color (8 bit ints).
        binned_spikes: spiketrain matrix of size t,cells
        stride: the stride of convolution (interval)
        depth: number of lags for convolution
        norm: boolean whether to normalize STA for nice plotting
    Returns:
        STA tensor, indexed by cell,t,y,x,color channel
    """

    # Intialize STA tensor.
    num_cells = binned_spikes.shape[1]
    n_frames,height,width,n_phosphors = stimulus.shape
    sta_tensor = np.zeros((num_cells,depth,height,width,n_phosphors))

    # Vectorize stimulus.
    stimulus = np.reshape(stimulus,(n_frames,
                                    height * width * n_phosphors)).astype(np.float32).T

    # Check that the FORTRAN type is reasonable.
    blas_obj = blas.get_blas_funcs('gemm',(stimulus,binned_spikes))

    if blas_obj.typecode != 'd':
        assert False, "The stimulus or spikes array were not of expected type."

    # Initialize STA, the offset range, and compute.
    offset_cnt = 0
    sta_tensor = np.zeros((num_cells,depth,height,width,n_phosphors))
    offset_range = np.arange(0,(depth * stride),stride)

    # Initialize progress bar.
    widgets = ['Computing STAs: ', Percentage(), ' ', Bar(marker='*',
               left='[',right=']'),' ', ETA()]
    pbar = ProgressBar(widgets=widgets, maxval=depth)
    pbar.start()

    for offset in offset_range:

        # Ensure the spike and stimulus vector keep their byte order.
        spikes = np.array(binned_spikes[offset:,:],order="F")

        assert spikes.flags.fortran,\
               "Could not get spikes vector in FORTRAN order for BLAS routines."
        assert stimulus[:,:n_frames-offset].flags.fortran,\
               "Could not get stimulus in FORTRAN order for BLAS routines."

        # Compute the STA and write to the tensor.
        sta = blas.dgemm(alpha=1.0, a=stimulus[:,:n_frames-offset],
                         b=spikes)

        # Index appropriately to avoid nans and normalize by number of spikes.
        nonzero_inds = np.nonzero(np.sum(spikes,0))[0]
        sta[:,nonzero_inds] /= np.sum(spikes,0)[nonzero_inds]

        del spikes

        sta_tensor[:,offset_cnt,:,:,:] = np.reshape(sta.T,
                                                  (num_cells,height,
                                                   width,n_phosphors))

        del sta

        pbar.update(offset_cnt)
        offset_cnt +=1

    pbar.finish()

    if not norm:
        return sta_tensor

    return normalize_sta(sta_tensor)

def compute_sta_torch(stimulus,binned_spikes,stride,depth=30,norm=False):
    """
    Computes STA by convolving the spiketrain with the stimulus. This
    implementation uses PyTorch. It doesn't use GPUs because of the memory
    intensiveness, but it is extremely fast relative to numpy even on CPU.

    Parameters:
        stimulus: stimulus tensor of size t, y, x, color (8 bit ints).
        binned_spikes: spiketrain matrix of size t,cells
        stride: the stride of convolution (interval)
        depth: number of lags for convolution
        norm: boolean whether to normalize STA for nice plotting
    Returns:
        STA tensor, indexed by cell,t,y,x,color channel
    """

    # Intialize STA tensor.
    num_cells = binned_spikes.shape[1]
    n_frames,height,width,n_phosphors = stimulus.shape
    sta_tensor = np.zeros((num_cells,depth,height,width,n_phosphors))

    # Vectorize stimulus.
    stimulus = np.reshape(stimulus,(n_frames,
                                    height * width * n_phosphors)).astype(np.uint8).T
    
    # Initialize STA, the offset range, and compute.
    offset_cnt = 0
    sta_tensor = np.zeros((num_cells,depth,height,width,n_phosphors))
    offset_range = np.arange(0,(depth * stride),stride)

    # Initialize progress bar.
    widgets = ['Computing STAs: ', Percentage(), ' ', Bar(marker='*',
               left='[',right=']'),' ', ETA()]
    pbar = ProgressBar(widgets=widgets, maxval=depth)
    pbar.start()

    # Set up torch stuff.
    device = torch.device('cpu')
    binned_spikes = torch.tensor(binned_spikes).to(device).to(torch.float32)
    stimulus = torch.tensor(stimulus).to(device).to(torch.float32)

    for offset in offset_range:

        with torch.no_grad():
            sta = torch.matmul(stimulus[:,:n_frames-offset],
                            binned_spikes[offset:,:])

            # Index appropriately to avoid nans and normalize by number of spikes.
            nonzero_inds = torch.nonzero(torch.sum(binned_spikes[offset:,:],
                                                   dim=0),as_tuple=False)
            sta[:,nonzero_inds] /= torch.sum(binned_spikes[offset:,:],
                                            dim=0)[nonzero_inds]

        sta = sta.detach().cpu().numpy() 
        sta_tensor[:,offset_cnt,...] = np.reshape(sta.T,
                                                  (num_cells,height,
                                                   width,n_phosphors))
        pbar.update(offset_cnt)
        offset_cnt +=1

    pbar.finish()

    if not norm:
        return sta_tensor

    return normalize_sta(sta_tensor)

def compute_stv(stimulus: np.ndarray, binned_spikes: np.ndarray, depth: int=61):
    """
    Computes the spike-triggered variance (STV) by convolving the spiketrain with the stimulus. This
    implementation uses PyTorch. It doesn't use GPUs because of the memory
    intensiveness, but it is extremely fast relative to numpy even on CPU.

    Parameters:
        stimulus: stimulus tensor of size t, y, x, color (8 bit ints).
        binned_spikes: spiketrain matrix of size t,cells
        stride: the stride of convolution (interval)
        depth: number of lags for convolution
    Returns:
        STV tensor, indexed by cell,t,y,x,color channel
    """

    # Intialize STA tensor.
    stride = 1
    num_cells = binned_spikes.shape[1]
    n_frames,height,width,n_phosphors = stimulus.shape
    sta_tensor = np.zeros((num_cells,depth,height,width,n_phosphors))

    # Vectorize stimulus.
    stimulus = np.reshape(stimulus,(n_frames,
                                    height * width * n_phosphors)).astype(np.uint8).T
    
    # Initialize STA, the offset range, and compute.
    offset_cnt = 0
    sta_tensor = np.zeros((num_cells,depth,height,width,n_phosphors))
    offset_range = np.arange(0,(depth * stride),stride)

    # Initialize progress bar.
    widgets = ['Computing STVs: ', Percentage(), ' ', Bar(marker='*',
               left='[',right=']'),' ', ETA()]
    pbar = ProgressBar(widgets=widgets, maxval=depth)
    pbar.start()

    # Set up torch stuff.
    device = torch.device('cpu')
    binned_spikes = torch.tensor(binned_spikes).to(device).to(torch.float32)
    stimulus = torch.tensor(stimulus).to(device).to(torch.float32)

    for offset in offset_range:

        with torch.no_grad():
            sta = torch.matmul(stimulus[:,:n_frames-offset],
                            binned_spikes[offset:,:])

            # Index appropriately to avoid nans and normalize by number of spikes.
            nonzero_inds = torch.nonzero(torch.sum(binned_spikes[offset:,:],
                                                   dim=0),as_tuple=False)
            sta[:,nonzero_inds] /= torch.sum(binned_spikes[offset:,:],
                                            dim=0)[nonzero_inds]

        sta = sta.detach().cpu().numpy() 
        sta_tensor[:,offset_cnt,...] = np.reshape(sta.T,
                                                  (num_cells,height,
                                                   width,n_phosphors))
        pbar.update(offset_cnt)
        offset_cnt +=1

    pbar.finish()

    return sta_tensor

def compute_timecourses(stas,norm=False,compute_mean=False,alpha=1/10000):
    """
    Computes STA timecourses by summing the primaries of the significant stixels.

    Parameters:
        STA tensor: size num_cells, t,  height,  width, color 
        norm: whether to center at 0 and divide by L2 norm (True) or just return
              the raw sum (False)
        compute_mean: whether to normalize the sum by the stixel number or not.
        alpha: significance level for significant stixels
    Returns:
        Array of size num_cells x, t,colors of each time course. If there was no
        time course for a cell, it is a matrix of zeros.

    Assumes that this time course is not flipped, and is flipped internally (i.e. same order
    as from the sta calculation in this module)
    """
    n_cells,sta_depth,_,_,n_phosphors = stas.shape
    timecourse_matrix = np.zeros((n_cells,sta_depth,n_phosphors))


    for i in range(n_cells):
        cell_sta = stas[i,...]
        sig_stixels = np.argwhere(calculate_sig_stixels(cell_sta - np.mean(cell_sta),
                                                        alpha=alpha)[0])

        if sig_stixels.shape[0] == 0:
            continue
    
        for stixel in sig_stixels:
            timecourse_matrix[i,...] += cell_sta[:,stixel[0],stixel[1],:]
        
        if compute_mean:
            timecourse_matrix[i,...] / sig_stixels.shape[0]

    if not norm:
        return timecourse_matrix
    
    # Mean subtract along each color channel and L2 normalize.
    timecourse_matrix -= np.mean(timecourse_matrix,axis=1,keepdims=True)

    # Get nonzero indices to avoid nans.
    nonzero_inds = ~np.all(timecourse_matrix == 0,axis=(1,2))
    timecourse_norm = linalg.norm(timecourse_matrix,axis=(1,2))[:,None,None]
    timecourse_matrix[nonzero_inds,...] /= timecourse_norm[nonzero_inds,...]

    return timecourse_matrix

def compute_spatial_maps(stas,celltypes_dict):
    """
    Computes spatial map by convolving the STA with the time course. For non-SBCs, the 
    grayscaled timecourse is used; for SBCs, only the Blue phosphor is.

    Parameters:
        stas: tensor of size num_cells,t,y,x,color
        celltypes_dict: dict mapping each cell to a type
    Returns:
        tensor of size num_cells,y,x,c
    """
    n_cells,_,y,x,c = stas.shape
    spatial_maps = []
    timecourses = compute_timecourses(stas,norm=False,compute_mean=True)
    cellids = sorted(list(celltypes_dict.keys()))

    for i in range(n_cells):
        cell = cellids[i]
        celltype = celltypes_dict[cell].lower()
        cell_sta = stas[i,...]
        sig_stixels = np.argwhere(calculate_sig_stixels(cell_sta -
                                                    np.mean(cell_sta))[0])

        if sig_stixels.shape[0] == 0:
            spatial_maps.append(np.zeros((y,x,c)))
            continue

        timecourse = timecourses[i,...]

        # If SBC, take blue, otherwise, grayscale. 
        if "sbc" in celltype or "blue" in celltype:
            timecourse = timecourse[...,-1]
        else:
            timecourse = np.mean(timecourse,axis=-1)

        timecourse -= np.mean(timecourse.ravel())

        # If the cell is "off", flip polarity.
        if "off" in celltype:
            timecourse *= -1
        
        spatial_map = np.dot(cell_sta.T - np.mean(cell_sta.ravel()),timecourse).T

        # Normalize to get back in a nice plotting format.
        spatial_map = (spatial_map /
                       np.max(np.abs(spatial_map)) / 2) +  cfg.STA_NORM_MEAN
        spatial_maps.append(spatial_map)
    
    return np.asarray(spatial_maps)

def compute_filters_svd(cell_sta, bin_rate=120.0, compute_rank_1=True):
    nt, ny, nx, n_phosphors = cell_sta.shape
    space_filter = np.zeros((ny, nx, n_phosphors))
    time_filter = np.zeros((nt,n_phosphors))
    r2 = np.zeros((n_phosphors,))
    lobe_pts = range(np.round(bin_rate*0.05).astype(int)-1,np.round(bin_rate*0.1).astype(int)-1)
    for i in range(n_phosphors):
        # Pass the STA through a Gaussian spatial filter.
        s = gaussian_filter(cell_sta[...,i], sigma=0.75)
        s_mat = s.reshape((nt, ny*nx))
        # Extract the spatiotemporal features using SVD.
        space, singular_values, time = np.linalg.svd(s_mat.T, full_matrices=False)
        # Get the space filter.
        space_filter[:,:,i] = space[:,0].reshape(ny,nx)
        time_filter[:,i] = time[0,:]
        # Divide by the max L2 norm.
        l2 = np.max(np.linalg.norm(time_filter[:,i], axis=0))
        time_filter[:,i] /= l2
        if compute_rank_1:
            # Get the rank-1 approximation
            r1 = np.outer(time[0,:], space[:,0]) * singular_values[0]
            # Get the r^2 value
            r2[i] = np.corrcoef(s_mat.ravel(), r1.ravel())[0,1]
    return space_filter, time_filter, r2

def compute_sta_params_fast(stas, snr_threshold=2.0):
    """
    Computes STA timecourses by summing the primaries of the significant stixels.

    Parameters:
        STA tensor: size num_cells, t,  height,  width, color 
        norm: whether to center at 0 and divide by L2 norm (True) or just return
              the raw sum (False)
        compute_mean: whether to normalize the sum by the stixel number or not.
        alpha: significance level for significant stixels
    Returns:
        Array of size num_cells x, t,colors of each time course. If there was no
        time course for a cell, it is a matrix of zeros.

    Assumes that this time course is not flipped, and is flipped internally (i.e. same order
    as from the sta calculation in this module)
    """
    n_cells,sta_depth,x,y,n_phosphors = stas.shape
    timecourse_matrix = np.zeros((n_cells,sta_depth,n_phosphors))
    spatial_maps = np.zeros((n_cells,x,y,n_phosphors))
    gauss_params = np.ones((n_cells,6))

    half_time = np.floor(sta_depth/2).astype(int) + 1

    # Calculate SNR.
    np.seterr(all='ignore')
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        for i in tqdm(range(n_cells), desc ="Computing RF parameters"):
            cell_sta = stas[i,...]
            time_course = calculate_time_course_svd(cell_sta, threshold=2.0)
            space_map = compute_spatial_map(cell_sta, time_course)

            rf_map = space_map.copy()
            for ii in range(rf_map.shape[2]):
                if np.nanmax(-rf_map[:,:,ii]) > np.nanmax(rf_map[:,:,ii]):
                    rf_map[:,:,ii] = -rf_map[:,:,ii]
            
            if np.ndim(rf_map) == 3:
                rf_map = np.nanmean(rf_map[:,:,(0,2)], axis=2)
            # Remove values less than zero.
            rf_map[rf_map < 0.0] = 0.0
            rf_map = rf_map / np.nanmax(rf_map)
            # Find the optimal threshold for the RF map.
            threshold = find_rf_threshold(rf_map)
            # Apply the threshold to the RF map.
            rf_map[rf_map < threshold] = 0.0
            rf_map[rf_map >= threshold] = 1.0
            sig_stixels = np.argwhere(rf_map > 0.0)
            if sig_stixels.shape[0] == 0:
                continue
            
            # time_course, sig_stixels = calculate_map_and_time_course(cell_sta, threshold=2.0)
            snr = np.divide(np.std(time_course[1:half_time,(0,2)],axis=0), np.std(time_course[half_time:,(0,2)],axis=0))
            snr[np.isnan(snr)] = 0.0 # Set NaN's to zero.
            snr = np.nanmax(snr)
            if (not (time_course == 0).all()) and (snr > snr_threshold) and (sig_stixels.shape[0] != 0):
                space_map = compute_spatial_map(cell_sta, time_course)
                try:
                    amplitude, theta, x0, y0, width_x, width_y = analysis_utils.gaussian_moments(space_map[...,0], sig_stixels)
                except:
                    amplitude=0.0; theta=0.0; x0=1.0; y0=1.0; width_x=1.0; width_y=1.0
                timecourse_matrix[i,...] = time_course
                spatial_maps[i,...] = space_map
                # amplitude, theta, xo, yo, sigma_x, sigma_y
                gauss_params[i,...] = [amplitude, theta, x0, y0, width_x, width_y]

    return timecourse_matrix, spatial_maps, gauss_params

def compute_rank1_rf_fit(cell_sta: np.ndarray, sig_stixels: np.ndarray) -> np.ndarray:
    """
    Computes the rank-1 approximation of the STA.

    Parameters:
        cell_sta: The cell's STA (size: t,y,x,color).
        sig_stixels: The significant stixels (size: n,2).

    Returns:
        The Pearson correlation between the raw STA and the rank-1 approximation of the STA.
    """
    # Fit the sta to the rank_1 appromation.
    nt, ny, nx, n_phosphors = cell_sta.shape
    r2 = np.zeros((n_phosphors,))

    if sig_stixels.shape[0] == 0:
        return r2
        
    s = np.zeros(cell_sta.shape)
    for i in range(n_phosphors):
        s[...,i] = gaussian_filter(cell_sta[...,i], sigma=1.0)

    time_filter = compute_timecourse_from_index(s, sig_stixels)
    s_map = compute_spatial_map(s, time_filter)

    for i in range(n_phosphors):
        s_mat = s_map[...,i].reshape((ny*nx,1))
        r1 = np.outer(time_filter[:,i].T, s_mat)
        r2[i] = np.corrcoef(s[...,i].reshape((nt,ny*nx)).ravel(), r1.ravel())[0,1]
    return r2

def compute_polygon_area(p: np.ndarray) -> float:
    """Compute the area of a polygon using the shoelace formula.
    https://en.wikipedia.org/wiki/Shoelace_formula
    Parameters:
        p: The vertices of the polygon. (np.ndarray)
    Returns:
        area: The area of the polygon. (float)
    """
    area = 0.0
    i_max = p.shape[0]-1
    for i in range(i_max):
        area += (p[i][0] * p[i+1][1]) - (p[i+1][0] * p[i][1])
    area += (p[i_max][0] * p[0][1]) - (p[0][0] * p[i_max][1])
    return np.abs(area / 2.0)

def compute_polygon_centroid(p: np.ndarray) -> Tuple[float, float]:
    """Compute the centroid of a polygon using the shoelace formula.
    https://en.wikipedia.org/wiki/Shoelace_formula
    Parameters:
        p: The vertices of the polygon. (np.ndarray)
    Returns:
        centroid_x: The x-coordinate of the centroid. (float)
        centroid_y: The y-coordinate of the centroid. (float)
    """  
    area = compute_polygon_area(p)

    if p.shape[0] < 1:
        return 1.0, 1.0
    elif (p.shape[0] < 3) | (area < 1.0):
        return np.mean(p[:,0]), np.mean(p[:,1])
    
    centroid_x = 0.0
    centroid_y = 0.0
    i_max = p.shape[0]-1
    for i in range(i_max):
        centroid_x += (p[i][0] + p[i+1][0]) * ((p[i][0] * p[i+1][1]) - (p[i+1][0] * p[i][1]))
        centroid_y += (p[i][1] + p[i+1][1]) * ((p[i][0] * p[i+1][1]) - (p[i+1][0] * p[i][1]))
    centroid_x += (p[i_max][0] + p[0][0]) * ((p[i_max][0] * p[0][1]) - (p[0][0] * p[i_max][1]))
    centroid_y += (p[i_max][1] + p[0][1]) * ((p[i_max][0] * p[0][1]) - (p[0][0] * p[i_max][1]))
    centroid_x /= (area * 6.0)
    centroid_y /= (area * 6.0)
    return np.abs(centroid_x), np.abs(centroid_y)

def compute_polygon_orientation(p: np.ndarray) -> float:
    # Compute the orientation of a polygon.
    if p.shape[0] < 2:
        return 0.0
    m, _ = np.polyfit(p[:,0], p[:,1], 1)
    return np.arctan(m)/np.pi*180

def compute_hull_parameters(hull_vertices: np.ndarray) -> np.ndarray:
    """Compute the parameters of the convex hull.
    Parameters:
        hull_vertices: The vertices of the convex hull [n,2]. (np.ndarray)
    Returns:
        x0: The x-coordinate of the centroid. (float)
        y0: The y-coordinate of the centroid. (float)
        sigma_x: The standard deviation of the RF in the x-direction. (float)
        sigma_y: The standard deviation of the RF in the y-direction. (float)
        orientation: The orientation of the hull in degrees. (float)
    """
    # Compute the parameters of the convex hull. This is meant to be a more reliable
    # way of computing the RF parameters than the Gaussian fit.
    np.seterr(all="ignore")
    warnings.filterwarnings('ignore')
    # Default values.
    x0 = 1.0
    y0 = 1.0
    sigma_x = 1.0
    sigma_y = 1.0
    orientation = 0.0
    if np.sum(hull_vertices) > 0.0:
        hull_vertices = hull_vertices[np.argwhere(hull_vertices[:,0] > 0)[:,0],:]
        x0, y0 = compute_polygon_centroid(hull_vertices)
        orientation = compute_polygon_orientation(hull_vertices)
        sigma_x = np.max([np.std(hull_vertices[:,0]), 0.5])
        sigma_y = np.max([np.std(hull_vertices[:,1]), 0.5])
        if sigma_x > 3.0*sigma_y:
            sigma_x = 3.0*sigma_y
        elif sigma_y > 3.0*sigma_x:
            sigma_y = 3.0*sigma_x
    return np.array([x0, y0, sigma_x, sigma_y, orientation])

def find_rf_threshold(rf_map: np.ndarray, thresholds: np.ndarray=np.logspace(np.log10(1.25), np.log10(7), 20)/10) -> float:
    """Find the optimal threshold for the RF map.
    Parameters:
        rf_map: The RF map to test. (np.ndarray)
        thresholds: The thresholds to test (np.ndarray)
    Returns:
        threshold: The optimal threshold. (float)
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        if np.ndim(rf_map) == 3:
            rf_map = np.nanmean(rf_map, axis=2)
        rf_map = rf_map / np.nanmax(rf_map)
    rf_map = np.nan_to_num(rf_map, nan=0.0, posinf=0.0, neginf=0.0)
    
    mean_island_area = np.zeros_like(thresholds)
    for idx, threshold in enumerate(thresholds):
        map_bool = rf_map >= threshold
        area = np.count_nonzero(map_bool)
        if area > 0:
            _, island_count = ndimage.label(map_bool)
            mean_island_area[idx] = area / island_count
    # Filter threshold options to avoid low thresh noisey RFs at the low mean RF size:
    decreasing = np.diff(mean_island_area) <= 0
    decreasing = np.concatenate([[True], decreasing])
    start_index = np.argmin(decreasing)
    mean_island_area = savgol_filter(mean_island_area, 5, 3, mode='nearest')

    best_index = np.argmax(mean_island_area[start_index:]) + start_index
    return thresholds[best_index]


def compute_spatiotemporal_maps(sta: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Computes STA timecourses by summing the primaries of the significant stixels.

    Parameters:
        STA tensor: size num_cells, t,  height,  width, color 
    Returns:
        Array of size num_cells x, t,colors of each time course. If there was no
        time course for a cell, it is a matrix of zeros.

    Assumes that this time course is not flipped, and is flipped internally (i.e. same order
    as from the sta calculation in this module)
    """
    n_cells,sta_depth,x,y,n_phosphors = sta.shape
    timecourse_matrix = np.zeros((n_cells,sta_depth,n_phosphors))
    spatial_maps = np.zeros((n_cells,x,y,n_phosphors))
    significance_maps = np.zeros((n_cells,x,y))
    hull_area = np.zeros((n_cells,))
    hull_vertices = np.zeros((n_cells,50,2))
    hull_parameters = np.zeros((n_cells,5))
    rank1_r2 = np.zeros((n_cells,3))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        for i in tqdm(range(n_cells), desc ="Computing RF parameters"):
            cell_sta = sta[i,...]

            # Take an initial pass through with a threshold of 2.0
            # time_course = calculate_time_course(cell_sta, threshold=2.0)
            if cell_sta.shape[0] > 3:
                time_course = calculate_time_course_count(cell_sta, max_stixels=3)
                # time_course = calculate_time_course_svd(cell_sta, threshold=2.0)
                space_map = compute_spatial_map(cell_sta, time_course)
            else:
                _, x, y, _ = cell_sta.shape
                time_course = np.zeros((cell_sta.shape[0], cell_sta.shape[3]))
                try:
                    sig_stixels, best_threshold, _ = find_best_map_threshold(cell_sta)
                    space_map = np.zeros((x, y, cell_sta.shape[3]))
                    if sig_stixels.shape[0] > 0:
                        significance_map = map_from_points(sig_stixels, x, y)
                        for stixel in sig_stixels:
                            time_course += cell_sta[:,stixel[0],stixel[1],:]
                        for cc in range(cell_sta.shape[3]):
                            space_map[:,:,cc] = cell_sta[1,:,:,cc].copy()
                            space_map[:,:,cc][significance_map == 0] = 0.0
                except Exception as e:
                    print('Error computing time course: ' + str(e))

            timecourse_matrix[i,...] = time_course
            spatial_maps[i,...] = space_map
            
            if cell_sta.shape[0] > 3:
                try:
                    significance_map, sig_stixels  = compute_significance_map(cell_sta, threshold=2.0)
                except Exception as e:
                    print(e)
                    sig_stixels = np.zeros((0,2))
                    significance_map = np.zeros((x,y))

            if sig_stixels.shape[0] == 0:
                continue
            
            significance_maps[i,...] = significance_map
            
            rank1_r2[i,:] = compute_rank1_rf_fit(cell_sta, sig_stixels)

            rf_map = space_map.copy()
            for ii in range(rf_map.shape[2]):
                if np.nanmax(-rf_map[:,:,ii]) > np.nanmax(rf_map[:,:,ii]):
                    rf_map[:,:,ii] = -rf_map[:,:,ii]
            
            if np.ndim(rf_map) == 3:
                rf_map = np.nanmean(rf_map[:,:,(0,2)], axis=2)
            # Remove values less than zero.
            rf_map[rf_map < 0.0] = 0.0
            rf_map = rf_map / np.nanmax(rf_map)
            # Find the optimal threshold for the RF map.
            threshold = find_rf_threshold(rf_map)
            # Apply the threshold to the RF map.
            rf_map[rf_map < threshold] = 0.0
            rf_map[rf_map >= threshold] = 1.0

            h_vert, h_area = compute_hull(rf_map) # compute_hull(significance_map)
            hull_area[i] = h_area
            hull_vertices[i,:h_vert.shape[0],:] = h_vert

            # # Recompute the time course and spatial map based on the hull estimate.
            # time_course = compute_timecourse_in_hull(cell_sta, h_vert)
            # if time_course is not None:
            #     timecourse_matrix[i,...] = time_course
            #     space_map = compute_spatial_map(cell_sta, time_course)
            #     spatial_maps[i,...] = space_map

            # Get the parameters of the convex hull.
            hull_parameters[i,:] = compute_hull_parameters(h_vert)

            _, thresh = cv2.threshold((255*significance_map).astype(np.uint8),0,255,0)
            cntrs = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
            cntrs = cntrs[0] if len(cntrs) == 2 else cntrs[1]
            if cntrs[0].shape[0] > 5:
                ellipse = cv2.fitEllipse(cntrs[0])
                hull_parameters[i,0] = ellipse[0][0]
                hull_parameters[i,1] = ellipse[0][1]
                hull_parameters[i,4] = ellipse[2]

            # Flip the spatia map sign if it's negative.
            for jj in range(timecourse_matrix[i,...].shape[1]):
                if (-np.min(timecourse_matrix[i,:,jj]) == np.max(np.abs(timecourse_matrix[i,:,jj]))):
                    spatial_maps[i,:,:,jj] *= -1

    return timecourse_matrix, spatial_maps, significance_maps, hull_parameters, hull_area, hull_vertices, rank1_r2

def compute_spatiotemporal_maps_legacy(stas, use_svd=True, compute_mean=False,alpha=1/10000, bin_rate=120.0):
    """
    Computes STA timecourses by summing the primaries of the significant stixels.

    Parameters:
        STA tensor: size num_cells, t,  height,  width, color 
        norm: whether to center at 0 and divide by L2 norm (True) or just return
              the raw sum (False)
        compute_mean: whether to normalize the sum by the stixel number or not.
        alpha: significance level for significant stixels
    Returns:
        Array of size num_cells x, t,colors of each time course. If there was no
        time course for a cell, it is a matrix of zeros.

    Assumes that this time course is not flipped, and is flipped internally (i.e. same order
    as from the sta calculation in this module)
    """
    n_cells,sta_depth,x,y,n_phosphors = stas.shape
    timecourse_matrix = np.zeros((n_cells,sta_depth,n_phosphors))
    spatial_maps = np.zeros((n_cells,x,y,n_phosphors))
    significance_maps = np.zeros((n_cells,x,y))
    hull_area = np.zeros((n_cells,))
    hull_vertices = np.zeros((n_cells,50,2))
    hull_parameters = np.zeros((n_cells,5))
    gauss_params = np.ones((n_cells,6))
    gabor_params = np.ones((n_cells,7))
    gauss_r2 = np.zeros((n_cells,1))
    gabor_r2 = np.zeros((n_cells,1))
    rank1_r2 = np.zeros((n_cells,3))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        for i in tqdm(range(n_cells), desc ="Computing RF parameters"):
            # print('Processing cell ' + str(i))
            cell_sta = stas[i,...]

            # Take an initial pass through with a threshold of 2.0
            time_course = calculate_time_course(cell_sta, threshold=2.0)
            space_map = compute_spatial_map(cell_sta, time_course)

            timecourse_matrix[i,...] = time_course
            spatial_maps[i,...] = space_map
            

            significance_map, sig_stixels  = compute_significance_map(cell_sta, threshold=2.0)

            significance_maps[i,...] = significance_map

            if sig_stixels.shape[0] == 0:
                continue
            
            rank1_r2[i,:] = compute_rank1_rf_fit(cell_sta, sig_stixels)

            h_vert, h_area = compute_hull(significance_map)
            hull_area[i] = h_area
            hull_vertices[i,:h_vert.shape[0],:] = h_vert

            # Get the parameters of the convex hull.
            hull_parameters[i,:] = compute_hull_parameters(h_vert)[0]

            # Compute the rotation angle of the points about the mean to guess at the theta
            if sig_stixels.shape[0] > 1:
                # xy_centered = sig_stixels.astype(float)
                # xy_centered[:,0] -= np.mean(xy_centered[:,0])
                # xy_centered[:,1] -= np.mean(xy_centered[:,1])
                # theta_guess = np.arctan(np.mean(xy_centered[:,1]/xy_centered[:,0]))/np.pi*180
                theta_guess = hull_parameters[i,4]
            else:
                theta_guess = 0.0

            # Fit a 2D Gaussian
            try:
                s_map = np.mean(spatial_maps[i,:,:,(0,2)], axis=0)
                s_map /= np.max(np.abs(s_map))
                p_gauss, r2 = analysis_utils.fit_gaussian_rf(s_map, theta_guess)
                gauss_params[i,:] = p_gauss
                gauss_r2[i] = r2
            except:
                pass
                # print('Gaussian fit failed for ' + str(i))

            # Fit a 2D Gabor
            try:
                s_map = np.mean(spatial_maps[i,:,:,(0,2)], axis=0)
                s_map /= np.max(np.abs(s_map))
                p_gabor, r2 = analysis_utils.fit_gabor_rf(s_map, theta_guess)
                gabor_params[i,:] = p_gabor
                gabor_r2[i] = r2
            except:
                pass
                # print('Gabor fit failed for ' + str(i))

            # Flip the spatia map sign if it's negative.
            for jj in range(timecourse_matrix[i,...].shape[1]):
                if (-np.min(timecourse_matrix[i,:,jj]) == np.max(np.abs(timecourse_matrix[i,:,jj]))):
                    spatial_maps[i,:,:,jj] *= -1

    return timecourse_matrix, spatial_maps, significance_maps, hull_area, hull_vertices, gauss_params, gabor_params, gauss_r2, gabor_r2, rank1_r2

def mask_sta(sta, sig_stixels,norm=True):
    """
    Arguments:
        sta: tensor of size t,y,x,color
        sig_stixels: array of stixel indices (0-indexed)
        norm: boolean indicating to mean subtract and L2 normalize among the sig stixels.
    Returns: STA with non-significant stixels set to 0.

    Sets all non significant stixels in an STA to zero.
    """
    masked_sta = np.zeros(sta.shape)
    stixel_vals = []

    for stixel in sig_stixels:
        masked_sta[:,stixel[0],stixel[1],:] = sta[:,stixel[0],stixel[1],:]
        stixel_vals.append(sta[:,stixel[0],stixel[1],:])
    
    stixel_vals = np.asarray(stixel_vals)

    if not norm: 
        return masked_sta
    
    # Normalize based on the cached values. 
    for stixel in sig_stixels:
        masked_sta[:,stixel[0],stixel[1],:] -= np.mean(stixel_vals.ravel())
        masked_sta[:,stixel[0],stixel[1],:] /= linalg.norm(stixel_vals.ravel())

    return masked_sta

def window_tensor(tensor,sig_stixels,vectorize=True):
    """
    Returns a windowed tensor based on the extrema of the sig stixels.
    
    Parameters:
        tensor: tensor (either stimulus or sta) of size t,y,x,color
        sig_stixels: array of sig stixels
        vectorize: boolean indicating whether to vectorize the tensor
    Returns:
        a truncated tensor based on the significant stixels.
    """

    # Get extrema of the window based on significant stixels.
    n_frames,height,width,n_phosphors = tensor.shape
    min_y = np.min(sig_stixels[:,0])
    max_y = np.minimum(np.max(sig_stixels[:,0])+1,height) # for boundaries.
    min_x = np.min(sig_stixels[:,1])
    max_x = np.minimum(np.max(sig_stixels[:,1])+1,width) # for boundaries.

    # Chop the tensor and vectorize. 
    tensor_window = tensor[:,min_y:max_y,min_x:max_x,:].copy()

    if not vectorize:
        return tensor_window

    _,window_height,window_width,_ = tensor_window.shape

    return np.reshape(tensor_window,(n_frames,
                                window_height * window_width * n_phosphors))

def compute_generator_signal_torch(stimulus, stas, stride, use_mask=False, alpha=1/10000, device_type='gpu'):
    """
    Computes the generator signal by convolving the STA with the stimulus. This
    is implemented in PyTorch (on CPU) for massive speedup relative to numpy.

    Parameters:
        stimulus: stimulus tensor of size t,y,x,color
        stas: tensor of size cells,n_frames,y,x,color
        stride: stride of the convolution
        alpha: significance level for sig stixels.
    Returns:
        generator signal matrix of size cells x frames
    """
    n_frames,_,_,_ = stimulus.shape
    n_cells,depth,height,width,n_phosphors = stas.shape

    # widgets = ['Computing generator signals: ', Percentage(), ' ', Bar(marker='*',
    #            left='[',right=']'),' ', ETA()]
    # pbar = ProgressBar(widgets=widgets, maxval=depth)
    # pbar.start()

    if device_type == 'cpu':
        device = torch.device('cpu')
    else:
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    masked_sta = np.reshape(stas.copy(),(n_cells,depth,
            height * width * n_phosphors))
    time_samples = stimulus.shape[0]
    stimulus_window = np.reshape(stimulus.copy(),(time_samples,
            height * width * n_phosphors))
    
    generator = torch.zeros(n_cells, n_frames).to(device).to(torch.float32)

    # Convert to floating point
    stimulus_window = torch.tensor(stimulus_window).to(device).to(torch.float32)
    masked_sta = torch.tensor(masked_sta).to(device).to(torch.float32)

    for offset in range(depth):
        with torch.no_grad():
            generator[:, (offset * stride):] +=\
                torch.matmul(masked_sta[:,offset, :],
                            stimulus_window[:n_frames - (offset * stride), :].T)
            # pbar.update(offset)
    generator = generator.detach().cpu().numpy()

    # pbar.finish()
    if device.type == 'cuda':
        torch.cuda.empty_cache()
    
    # Clear some memory
    del masked_sta, stimulus_window
    if device.type == 'cuda':
        torch.cuda.empty_cache()
    gc.collect()
    
    return np.asarray(generator)

def compute_generator_signal_torch_old(stimulus, stas, stride, use_mask=False, alpha=1/10000, device_type='gpu'):
    """
    Computes the generator signal by convolving the STA with the stimulus. This
    is implemented in PyTorch (on CPU) for massive speedup relative to numpy.

    Parameters:
        stimulus: stimulus tensor of size t,y,x,color
        stas: tensor of size cells,n_frames,y,x,color
        stride: stride of the convolution
        alpha: significance level for sig stixels.
    Returns:
        generator signal matrix of size cells x frames
    """
    n_frames,_,_,_ = stimulus.shape
    widgets = ['Computing generator signals: ', Percentage(), ' ', Bar(marker='*',
               left='[',right=']'),' ', ETA()]
    pbar = ProgressBar(widgets=widgets, maxval=stas.shape[0])
    pbar.start()
    generator_list = []
    depth = stas.shape[1]

    if device_type == 'cpu':
        device = torch.device('cpu')
    else:
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    # Normalize to 0,1 if required: assumes that this is in 8 bit format.
    # if np.max(stas.ravel()) > 1:
    #     stas /= cfg.MONITOR_BIT_DEPTH # to 0,1

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        for i in range(stas.shape[0]):
            cell_sta = stas[i,...]
            sig_stixels = np.argwhere(calculate_sig_stixels(cell_sta - np.mean(cell_sta), alpha=alpha)[0])

            if sig_stixels.shape[0] == 0:
                generator_list.append(np.asarray([np.nan for _ in range(n_frames)]))
                pbar.update(i)
                continue
            
            # Get the masked STA and window the stimulus and STA.
            if use_mask:
                masked_sta = mask_sta(cell_sta,sig_stixels)
                masked_sta = window_tensor(masked_sta,sig_stixels)
                stimulus_window = window_tensor(stimulus,sig_stixels)
            else:
                masked_sta = cell_sta.copy()
                n_frames,window_height,window_width,n_phosphors = masked_sta.shape
                masked_sta = np.reshape(masked_sta,(n_frames,
                                window_height * window_width * n_phosphors))
                n_frames,window_height,window_width,n_phosphors = stimulus.shape
                stimulus_window = np.reshape(stimulus.copy(),(n_frames,
                                window_height * window_width * n_phosphors))

            # L2 normalize.
            masked_sta /= np.linalg.norm(masked_sta)
            generator = torch.zeros(n_frames).to(device).to(torch.float32)

            # Convert to floating point
            stimulus_window = torch.tensor(stimulus_window).to(device).to(torch.float32)
            masked_sta = torch.tensor(masked_sta).to(device).to(torch.float32)

            for offset in range(depth):
                with torch.no_grad():
                    generator[(offset * stride):] +=\
                        torch.matmul(stimulus_window[:n_frames - (offset * stride),:],
                                    masked_sta[offset,:])

            pbar.update(i)
            generator = generator.detach().cpu().numpy()
            generator_list.append(generator)

    pbar.finish()
    if device.type == 'cuda':
        torch.cuda.empty_cache()
    
    # Clear some memory
    del generator, masked_sta, stimulus_window
    gc.collect()
    
    return np.asarray(generator_list)

def bin_nonlinearity_count(prediction, response, num_bins=25, order='descend'):
    if order == 'descend':
        inds = np.argsort(-prediction)
    else:
        inds = np.argsort(prediction)
    x_sort = prediction[inds]
    y_sort = response[inds]
    vals_per_bin = np.floor(len(x_sort) / num_bins).astype(int)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        x_bin = np.mean(np.reshape(x_sort[:num_bins*vals_per_bin], (num_bins, vals_per_bin)).T, axis=0)
        y_bin = np.mean(np.reshape(y_sort[:num_bins*vals_per_bin], (num_bins, vals_per_bin)).T, axis=0)
    return x_bin, y_bin

def bin_nonlinearity_hist(prediction, response, num_bins=25):
    edges = np.histogram_bin_edges(prediction, bins=num_bins)
    x_bin = np.zeros((num_bins,))
    y_bin = np.zeros((num_bins,))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        for b in range(num_bins):
            inds = np.argwhere(np.logical_and(prediction >= edges[b],
                                            prediction < edges[b+1])).flatten()
            x_bin[b] = np.nanmean(prediction[inds])
            y_bin[b] = np.nanmean(response[inds])
    return x_bin, y_bin

def bin_nonlinearities(generator, response, method='count', num_bins=25, normalize=True):
    """
    Computes nonlinearity for LNP model by relating firing rate and generator signal.

    Parameters:
        generator: matrix of size n_cells, n_frames
        binned_spikes: matrix of size n_cells, n_frames
    Returns:
        tuple mean generator and mean spike counts matrices, each of size 
        num_cells,n_frames
    """
    # Normalize the generator signals by the nonlinearity.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        if normalize:
            generator /= np.std(generator, axis=1, keepdims=True)
        np.zeros((generator.shape[0], num_bins))
        x_bin = np.zeros((generator.shape[0], num_bins))
        y_bin = np.zeros((generator.shape[0], num_bins))
        for i in range(generator.shape[0]):
            if np.all(np.isfinite(generator[i,:])):
                if method == 'hist':
                    x_bin[i,:], y_bin[i,:] = bin_nonlinearity_hist(generator[i,:], response[i,:], num_bins=num_bins)
                else:
                    x_bin[i,:], y_bin[i,:] = bin_nonlinearity_count(generator[i,:], response[i,:], num_bins=num_bins)
    return x_bin, y_bin

def compute_nonlinearity(generator_signals, binned_spikes,
                         monitor_fs, stdz_g=True, step=12):
    """
    Computes nonlinearity for LNP model by relating firing rate and generator signal.

    Parameters:
        generator_signals: matrix of size n_cells, n_frames
        binned_spikes: matrix of size n_cells, n_frames
        monitor_fs: refresh freq. (full res.) of the monitor in Hz
        stdz_g: boolean indicating whether to do standard normalization on 
                generator signal (recommended!)
        step: step sized used when estimating nonlinearity
    Returns:
        tuple mean generator and mean spike counts matrices, each of size 
        num_cells,n_frames
    """
    # Standard normalization on generator signal if indicated.
    if stdz_g:
        std = np.std(generator_signals, axis=1, keepdims=True)
        generator_signals /= std

    mean_spike_counts = []
    mean_generator = []
    percentiles = np.linspace(0,100,step)

    for i in range(generator_signals.shape[0]):
        spike_counts_tmp = []
        mean_generator_tmp = []

        for j in range(len(percentiles)-1):
            lower = np.percentile(generator_signals[i,:], percentiles[j])
            inds = np.argwhere(np.logical_and(generator_signals[i,:] >= lower,
                                    generator_signals[i,:] < upper)).flatten()
            
            # Compute the mean in these bins.
            f = np.mean(binned_spikes[i,inds]) * monitor_fs # to get into Hz
            g = np.mean(generator_signals[i,inds])
            spike_counts_tmp.append(f)
            mean_generator_tmp.append(g)
        
        mean_spike_counts.append(np.asarray(spike_counts_tmp))
        mean_generator.append(np.asarray(mean_generator_tmp))
    
    return np.asarray(mean_generator),np.asarray(mean_spike_counts)
