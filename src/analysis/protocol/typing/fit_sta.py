import numpy as np
from scipy.stats import chi2, norm
from scipy.ndimage import convolve
from lnp import compute_spatial_map, map_from_points, compute_hull, compute_hull_parameters
import numpy as np
from scipy.optimize import curve_fit, minimize

# Helper function for robust standard deviation
def robust_std(data: np.ndarray, method: str='mad'):
    """
    Compute robust measures of standard deviation.

    Parameters:
        data (array-like): Input data.
        method (str): Method for robust standard deviation:
            - 'mad' : Median Absolute Deviation (MAD) scaled to match std
            - 'iqr' : Interquartile Range scaled to match std
            - 'std' : Basic standard deviation (not robust)
        scale_factor (float): Scaling factor for MAD (default = 1.4826 for Gaussian).
        iqr_factor (float): Scaling factor for IQR to estimate std (default = 0.7413 for Gaussian).
    
    Returns:
        float: Robust standard deviation.
    """
    data = np.asarray(data)
    if method == 'mad': # Median absolute deviation (MAD)
        median = np.median(data, axis=0)
        deviation = np.abs(data - median)
        return 1.4826 * np.median(deviation, axis=0)
    elif method == 'iqr':
        # Interquartile Range
        q75, q25 = np.percentile(data, [75, 25])
        return 0.7413 * (q75 - q25)
    elif method == 'std':
        # Standard deviation (non-robust)
        return np.std(data, ddof=1)
    else:
        raise ValueError("Invalid method. Choose 'mad', 'iqr', or 'std'.")

# Helper function to normalize vectors to unit sphere
def normalize_to_unit_sphere(vectors):
    norms = np.linalg.norm(vectors, axis=1, keepdims=True)
    return vectors / norms

# find_local_maxima(rf_test,'thresh',abs_thresh,'radius',params.radius,'return','matrix');
def find_local_maxima(rf_test, threshold: float, radius: float=1, return_type: str='indices'):
    search_locations = np.argwhere(rf_test > threshold)
    # Go through each location and check the surrounding pixels for local maxima.
    for loc in search_locations:
        pass
    raise NotImplementedError("Local maxima selection is not implemented.")

# Main function
def significant_stixels(sta: np.ndarray, selection_type: str='threshold', threshold: float=5, radius: float=1, strength_type: str='inner', spatial_filter: np.ndarray=None):
    """ Compute significant spatiotemporal pixels in a spatiotemporal array.
    Parameters:
        sta: Spatiotemporal array of shape (t, y, x, c).
        selection_type: Method for selecting significant stixels. Options are 'threshold' and 'maxima'.
        threshold: Threshold for selecting significant stixels.
        radius: Radius for selecting significant stixels.
        strength_type: Method for computing the strength of the stixel. Options are 'vector length' and 'inner'.
        spatial_filter: Spatial filter to apply to the STA before computing significant stixels.
    Returns:
        sig_stixels: Array of significant stixels.
        params: Dictionary of parameters.
        rf_strength: Array of RF strengths. 
    """
    # Handle empty STA
    if sta is None or sta.size == 0:
        return None, {}, None

    ## Check this code!!
    # Apply spatial filter
    if spatial_filter is not None:
        filt = spatial_filter / np.sqrt(np.sum(spatial_filter ** 2))
        cell_sta = convolve(sta, filt, mode='reflect')
    else:
        cell_sta = sta.copy()

    n = cell_sta.shape[0]
    rf = np.std(cell_sta, axis=0, ddof=1)
    rf = norm.ppf(chi2.cdf((n - 1) * rf.ravel() ** 2, df=(n - 1)))
    rf[rf == np.inf] = 5.5
    rf = rf.reshape(cell_sta.shape[1], cell_sta.shape[2], cell_sta.shape[3])

    # Re-center each color's distribution at 0
    for ii in range(cell_sta.shape[3]):
        rf[:,:,ii] -= np.median(rf[:,:,ii])
        
    rf_renorm = np.zeros(rf.shape)

    # Normalize by the noise standard deviation
    noise_sigmas = np.zeros(3)
    for ii in range(3):
        noise_sigmas[ii] = robust_std(rf[:,:,ii].ravel())
        rf_renorm[:,:,ii] = rf[:,:,ii] / noise_sigmas[ii]

    # Compute RF strength
    if strength_type == 'vector length':
        rf_strength = np.sqrt(np.sum(rf_renorm ** 2, axis=2))
        rf_test = rf_strength
        f_thresh = lambda th: np.sqrt(chi2.ppf(norm.cdf(th), df=rf.shape[2]))
    elif strength_type == 'inner':
        v = np.array([1, 1, 1])
        rf_dot_prod = rf_renorm.reshape(-1, rf.shape[2]) @ v
        rf_dot_prod /= robust_std(rf_dot_prod)
        rf_strength = np.abs(rf_dot_prod.reshape(rf.shape[0], rf.shape[1]))
        rf_test = rf_strength
        f_thresh = lambda th: -norm.ppf(norm.cdf(-th) / 2)
    elif strength_type == 'inner or':
        v = np.array([1, 1, 1])
        rf_dot_prod = rf_renorm.reshape(-1, rf.shape[2]) @ v
        rf_dot_prod /= robust_std(rf_dot_prod)
        rf_strength = np.abs(rf_dot_prod.reshape(rf.shape[0], rf.shape[1]))
        rf_test = rf_strength
        f_thresh = lambda th: -norm.ppf(norm.cdf(-th) / 2)
    else:
        raise ValueError(f"Unknown strength type: {strength_type}")    

    # Threshold for significant stixels.
    abs_thresh = f_thresh(threshold) if selection_type in ['threshold', 'maxima'] else None

    if selection_type == 'threshold':
        sig_stixels = np.argwhere(rf_test > abs_thresh)
    elif selection_type == 'maxima':
        # sig_stixels = find_local_maxima(rf_test,'thresh',abs_thresh,'radius',params.radius,'return','matrix');
        raise NotImplementedError("Local maxima selection is not implemented.")

    return sig_stixels, noise_sigmas, rf_strength

# Get the time course from the significant stixels.
def time_course_from_sta(cell_sta: np.ndarray, sig_stixels: np.ndarray, color_channel: str=None) -> np.ndarray:
    """
    Get the time course from the significant stixels.

    Parameters:
        cell_sta: The spatiotemporal receptive field
        sig_stixels: Coordinates of the significant stixels

    Returns: The time course from the significant stixels.
    """
    time_course = np.zeros((cell_sta.shape[0], cell_sta.shape[3]))
    if (sig_stixels is None) or (sig_stixels.shape[0] == 0):
        return time_course
    for stixel in sig_stixels:
        time_course += cell_sta[:,stixel[0],stixel[1],:]
    time_course /= sig_stixels.shape[0]
    if color_channel is not None:
        color_idx = 0 if color_channel == 'red' else 1 if color_channel == 'green' else 2
        time_course = time_course[:, color_idx]
    return time_course

def linear_filter_function(t, *params):
    """ Temporal filter function. 
    Arguments:
        t (np.ndarray): Time vector.
        params (list): List of parameters.
    Returns:
        f (np.ndarray): Filter output.
    """
    f = params[0] * (((t/params[1]) ** params[5]) / 
        (1+( (t/params[1]) ** params[5] ))) * np.exp(-((t/params[2]))) * np.cos(((2*np.pi*t) / params[3])+(2*np.pi*params[4]/360))
    return f

def temporal_filter_function(params, t):
    """ Temporal filter function. 
    Arguments:
        t (np.ndarray): Time vector.
        params (list): List of parameters.
    Returns:
        f (np.ndarray): Filter output.
    """
    f = params[0] * (((t/params[1]) ** params[5]) / 
        (1+((t/params[1])**params[5]))) * np.exp(-((t/params[2]))) * np.cos(((2*np.pi*t) / params[3])+(2*np.pi*params[4]/360))
    return f

def temporal_cost_function(params, t, data):
    """Cost function to minimize: sum of squared residuals."""
    model = temporal_filter_function(params, t)
    return np.sum((data - model)**2)

def fit_temporal_filters(time_course, bin_rate=120.0):
    t = np.arange(0, time_course.shape[0]) / bin_rate
    time_fits = np.zeros((time_course.shape))
    # Define the boundaries.
    bounds = [(-1000, 1000), (0.01, 0.3), (0.005, 0.5), (0.01, 1.0), (-360, 360), (4.0, 10.0)]
    # Define the fitting options.
    options = {'maxiter': 1e05, 'disp': False, 'maxfev': 1e05}
    # non_bds = ([-1000, 0.01, 0.005, 0.01, -360, 4.0], 
    #        [1000, 0.3, 0.5, 1.0, 360.0, 10.0])
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
        # initial_guess = [sgn*15, filter_peak/bin_rate/2, 0.06, 0.6, 90.0, 5.0]
        initial_guess = p0 = [sgn*15, 0.061433021205849,  0.016320868351090,  0.560874309169834, 52.432435563319949, 4.112884082863220]
        result = minimize(temporal_cost_function, initial_guess, args=(t, time_course[:,j].T), bounds=bounds, options=options, method='Nelder-Mead')
        gauss_params = result.x
        # gauss_params, _ = curve_fit(f=temporal_filter_function, xdata=t, ydata=time_course[:,j].T, p0=p0, bounds=non_bds, max_nfev=1e05)
        time_fits[:,j] = temporal_filter_function(gauss_params, t).T
        t_params[j,:]  = gauss_params
    return time_fits, t_params

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
        # p0 = [sgn*7, filter_peak/bin_rate, 0.06, 0.6, 90.0, 5.0]
        p0 = [sgn*15, 0.061433021205849,  0.016320868351090,  0.560874309169834, 52.432435563319949, 4.112884082863220]
        # if lobe_pts > 0.0:
        #     p0 = [62.562745558862289,  0.061433021205849,  0.016320868351090,  0.560874309169834, 52.432435563319949, 4.112884082863220]
        # else:
        #     p0 = [-62.562745558862289,  0.061433021205849,  0.016320868351090,  0.560874309169834, 52.432435563319949, 4.112884082863220]
        gauss_params, _ = curve_fit(f=linear_filter_function, xdata=t, ydata=time_course[:,j].T, p0=p0, bounds=non_bds, max_nfev=1e05)
        time_fits[:,j] = linear_filter_function(t, *gauss_params).T
        t_params[j,:]  = gauss_params
    return time_fits, t_params

# Define the 2D Gaussian function
def gaussian_2d(params, x, y):
    """2D Gaussian function."""
    if len(params) == 6:
        amplitude, xo, yo, sigma_x, sigma_y, theta = params
        surround_weight = 0.0
        surround_scale = 0.0
    elif len(params) == 8:
        amplitude, xo, yo, sigma_x, sigma_y, theta, surround_weight, surround_scale = params
    # amplitude, xo, yo, sigma_x, sigma_y, theta, = params
    xo, yo = float(xo), float(yo)  # Ensure xo and yo are floats
    a = (np.cos(theta)**2) / (2 * sigma_x**2) + (np.sin(theta)**2) / (2 * sigma_y**2)
    b = -(np.sin(2 * theta)) / (4 * sigma_x**2) + (np.sin(2 * theta)) / (4 * sigma_y**2)
    c = (np.sin(theta)**2) / (2 * sigma_x**2) + (np.cos(theta)**2) / (2 * sigma_y**2)
    g = amplitude * np.exp(-(a * ((x - xo)**2) +
                                      2 * b * (x - xo) * (y - yo) +
                                      c * ((y - yo)**2)))
    if len(params) == 6:
        return g
    if len(params) == 8:
        a = (np.cos(theta)**2) / (2 * (surround_scale*sigma_x)**2) + (np.sin(theta)**2) / (2 * (surround_scale*sigma_y)**2)
        b = -(np.sin(2 * theta)) / (4 * (surround_scale*sigma_x)**2) + (np.sin(2 * theta)) / (4 * (surround_scale*sigma_y)**2)
        c = (np.sin(theta)**2) / (2 * (surround_scale*sigma_x)**2) + (np.cos(theta)**2) / (2 * (surround_scale*sigma_y)**2)
        g2 = amplitude * np.exp(-(a * ((x - xo)**2) +
                                        2 * b * (x - xo) * (y - yo) +
                                        c * ((y - yo)**2)))
        return g + surround_weight * g2

# Define the spatial cost function to minimize
def spatial_cost_function(params, x, y, data):
    """Cost function to minimize: sum of squared residuals."""
    model = gaussian_2d(params, x, y)
    return np.sum((data - model)**2)

# Fit a 2D Gaussian to data
def fit_2d_gaussian(x, y, data, initial_guess, bounds=None):
    """Fit a 2D Gaussian to the given data using the Nelder-Mead method.
    Parameters:
        x (array-like): X-coordinates of the data.
        y (array-like): Y-coordinates of the data.
        data (array-like): Data to fit.
        initial_guess (array-like): Initial guess for the parameters.
        bounds (array-like): Bounds for the parameters.
    Returns:
        result: The optimization result.
    """
    # Make sure the initial guess is within bounds.
    if bounds is not None:
        for ii, (low, high) in enumerate(bounds):
            if initial_guess[ii] < low:
                initial_guess[ii] = low
            if initial_guess[ii] > high:
                initial_guess[ii] = high
    result = minimize(spatial_cost_function, initial_guess, args=(x, y, data), bounds=bounds, method='Nelder-Mead')
    return result

def points_in_ellipse(x, y, center, axes, theta):
    """
    Find the points that fall within an ellipse.
    
    Parameters:
        x (array-like): X-coordinates of the points.
        y (array-like): Y-coordinates of the points.
        center (tuple): (x_center, y_center) coordinates of the ellipse center.
        axes (tuple): (semi_major_axis, semi_minor_axis) of the ellipse.
        theta (float): Rotation angle of the ellipse in radians.
    
    Returns:
        array: Indices of the points that fall within the ellipse.
    """
    # Translate points to the ellipse's center
    x_centered = x - center[0]
    y_centered = y - center[1]
    
    # Rotate points to align with the ellipse's axes
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    x_rotated = x_centered * cos_theta + y_centered * sin_theta
    y_rotated = -x_centered * sin_theta + y_centered * cos_theta
    
    # Check if points are within the ellipse equation
    semi_major, semi_minor = axes
    ellipse_eq = (x_rotated / semi_major) ** 2 + (y_rotated / semi_minor) ** 2
    inside_indices = np.where(ellipse_eq <= 1)[0]  # Indices where the condition is met
    points_inside = list(zip(x[inside_indices], y[inside_indices]))
    # points_inside = np.array(points_inside)
    return points_inside

def estimate_spatial_params(sig_stixels: np.ndarray) -> np.ndarray:
    """ Estimate the 2D Gaussian parameters from the significant stixels.
    Parameters:
        sig_stixels: Array of significant stixels.
    Returns:
        Array of estimated parameters: [x_mean, y_mean, sd_x, sd_y, rotation_angle].
    """
    if sig_stixels.size == 0:
        return np.zeros(5)
    elif sig_stixels.shape[0] == 1:
        return np.array([sig_stixels[0,1], sig_stixels[0,0], 1, 1, 0])
    try:
        centered_stixels = sig_stixels.copy()
        centered_stixels = centered_stixels[:,::-1]
        # Center the significant stixels.
        xy_mean = np.mean(centered_stixels, axis=0)
        centered_stixels = centered_stixels - xy_mean
        U, S, Vt = np.linalg.svd(centered_stixels[:,::-1], full_matrices=True, compute_uv=True, hermitian=False)
        eigenvalues = S**2 / (centered_stixels.shape[0]-1)
        V = Vt.T
        # rotation_angle = V[1,0] * np.pi
        rotation_angle = np.arctan2(V[1, 0], V[0, 0])
        sd_x = np.sqrt(eigenvalues[0])
        sd_y = np.sqrt(eigenvalues[1])
    except Exception as e:
        return np.array([sig_stixels[0,1], sig_stixels[0,0], 1, 1, 0])
    return np.array([xy_mean[0], xy_mean[1], sd_x*1.2, sd_y*1.2, rotation_angle])

def fit_spatial_rf(spatial_rf: np.ndarray, sig_stixels: np.ndarray, fit_surround: bool=False) -> np.ndarray:
    """ Fit a 2D Gaussian to the spatial RF.
    Parameters:
        spatial_rf: The spatial RF.
        sig_stixels: The significant stixels.
    Returns:
        Array of estimated parameters: [amplitude, x_mean, y_mean, sd_x, sd_y, rotation_angle].
    """
    x = np.arange(spatial_rf.shape[1])
    y = np.arange(spatial_rf.shape[0])
    x, y = np.meshgrid(x, y)
    if np.max(spatial_rf) > np.abs(np.min(spatial_rf)):
        spatial_amp = np.max(spatial_rf)
    else:
        spatial_amp = np.min(spatial_rf)
    # Estimate the spatial parameters.
    start_params = estimate_spatial_params(sig_stixels)
    if fit_surround:
        initial_guess = [spatial_amp, start_params[0], start_params[1], start_params[2], start_params[3], start_params[4], -0.15, 2.0]
        bounds = [(-1, 1), (0, spatial_rf.shape[1]), (0, spatial_rf.shape[0]), (0.2, spatial_rf.shape[1]/3), (0.2, spatial_rf.shape[0]/3), (-np.pi, np.pi), (-0.5, 0.0), (1.2, 2.5)]
    else:
        initial_guess = [spatial_amp, start_params[0], start_params[1], start_params[2], start_params[3], start_params[4]]
        bounds = [(-1, 1), (0, spatial_rf.shape[1]), (0, spatial_rf.shape[0]), (0.2, spatial_rf.shape[1]/3), (0.2, spatial_rf.shape[0]/3), (-np.pi, np.pi)]
    result = fit_2d_gaussian(x, y, spatial_rf, initial_guess, bounds)
    spatial_params = result.x
    return spatial_params

def compute_rf_parameters(cell_sta: np.ndarray, fit_type: str='spatial') -> np.ndarray:
    """ Compute the RF parameters.
    Parameters:
        cell_sta: The spatiotemporal receptive field.
        fit_type: The type of fit to perform. Options are 'spatial' and 'full'.
    Returns:
        Array of estimated parameters: [amplitude, x_mean, y_mean, sd_x, sd_y, rotation_angle].
    """
    sig_stixels, _, _ = significant_stixels(cell_sta)
    time_course = time_course_from_sta(cell_sta, sig_stixels)
    space_map = compute_spatial_map(cell_sta, time_course)
    spatial_rf = np.mean(space_map, axis=2)
    spatial_rf /= np.max(np.abs(spatial_rf))
    # Fit the spatial RF.
    gauss_params = fit_spatial_rf(spatial_rf, sig_stixels=sig_stixels)
    # Find the points inside the ellipse
    x = np.arange(spatial_rf.shape[1])
    y = np.arange(spatial_rf.shape[0])
    x, y = np.meshgrid(x, y)
    points_inside = points_in_ellipse(x.ravel(), y.ravel(), (gauss_params[1], gauss_params[2]), (gauss_params[3], gauss_params[4]), gauss_params[5])
    points_inside = np.array(points_inside)
    points_inside = points_inside[:,(1,0)]
    # Find the intersection of the points inside the ellipse and the significant stixels
    intersecting_points = set(zip(*points_inside.T)) & set(zip(*sig_stixels.T))#set(map(tuple, sig_stixels)) & set(map(tuple, points_inside))
    intersecting_points = np.array(list(intersecting_points))
    # If the intersection is empty, use the significant stixels.
    if intersecting_points.shape[0] > 0:
        # Recompute the temporal filter based on the fit.
        time_course = time_course_from_sta(cell_sta, intersecting_points)
        # Re-fit the spatial RF.
        space_map = compute_spatial_map(cell_sta, time_course)
        spatial_rf = np.mean(space_map, axis=2)
        spatial_rf /= np.max(np.abs(spatial_rf))
        gauss_params = fit_spatial_rf(spatial_rf, sig_stixels=sig_stixels)
    if fit_type == 'spatial':
        pass
    else: # Fit the full spatiotemporal RF
        # Fit the time course.
        time_fits, t_params = fit_linear_filters(time_course, bin_rate=120.0)
        pass
    return gauss_params





