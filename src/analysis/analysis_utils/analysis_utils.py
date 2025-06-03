"""
Module of functions for analysis of spiking properties: EIs, ACFs, etc.
"""
import sys, os
# sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import numpy.ma as ma
from numpy import pi
# from gaussfitter.mpfit import mpfit
import scipy
from scipy.optimize import curve_fit
from scipy.interpolate import splrep, sproot, splev

def compute_pairwise_correlation(m_1: np.ndarray, m_2: np.ndarray, mask_invalid: bool=False) -> np.ndarray:
    """ Compute the correlation between the EI maps of the two chunks at all possible shifts by one electrode.
    Parameters:
        m_1: The first matrix.
        m_2: The second matrix.
    Returns:
        corr: The correlation between the EI maps of the two chunks at all possible shifts by one electrode.
    """
    # Make sure the matrices are the same size.
    if m_1.shape[1] != m_2.shape[1]:
        raise ValueError('Matrices must have the same number of columns.')
    if mask_invalid:
        m_1 = ma.masked_invalid(m_1)
        m_2 = ma.masked_invalid(m_2)
    # Subtract the mean.
    with np.errstate(invalid='ignore'):
        m_1 = m_1 - np.nanmean(m_1,axis=1)[:,np.newaxis]
        m_2 = m_2 - np.nanmean(m_2,axis=1)[:,np.newaxis]
    if mask_invalid:
        numerator = ma.dot(m_1, m_2.T)
    else:
        numerator = np.matmul(m_1, m_2.T)
    # Compute the L2 norm along the columns.
    norm_1 = np.linalg.norm(np.nan_to_num(m_1), axis=1)[:,np.newaxis]
    norm_2 = np.linalg.norm(np.nan_to_num(m_2), axis=1)[:,np.newaxis]
    denominator = np.matmul(norm_1, norm_2.T)
    with np.errstate(divide='ignore', invalid='ignore'):
        corr = numerator / denominator
        corr[denominator == 0] = 0
    return corr

# Cosine of argument in degrees
def cosd(x: float) -> float:
    """ Cosine of argument in degrees. 
    Parameters:
        x: The argument in degrees.
    Returns:
        The cosine of the argument.
    """
    return np.cos( np.deg2rad(x) )

# Sine of argument in degrees
def sind(x: float) -> float:
    """ Sine of argument in degrees.
    Parameters:
        x: The argument in degrees.
    Returns:
        The sine of the argument.
    """
    return np.sin( np.deg2rad(x) )

# Full-width at half-max
def fwhm(x: np.ndarray, y: np.ndarray):
    """ Full-width at half-max.
    Parameters:
        x: The x values.
        y: The y values.
    Returns:
        The full-width at half-max.
    """
    try:
        s = splrep(x, y - 0.5)
        roots = sproot(s)
        out = np.abs(roots[1] - roots[0])
    except:
        out = 1.0
    return out

def gaussian_moments(data: np.ndarray, sig_stixels: np.ndarray = None):
    """ Get the Gaussian moments. 
    Parameters:
        data: The data (np.ndarray).
        sig_stixels: The significant stixels (np.ndarray).
    Returns:
        amplitude: The amplitude.
        theta: The angle.
        x0: The x center.
        y0: The y center.
        width_x: The width in the x direction.
        width_y: The width in the y direction.
    """
    # Guess at the center of the receptive field.
    max_idx = np.unravel_index(np.argmax(np.abs(data), axis=None), data.shape)
    # Get the amplitude at the maximum.
    amplitude = data[max_idx]
    col = data[int(max_idx[0]), :] / amplitude
    row = data[:, int(max_idx[1])] / amplitude
    # Get the full-width at half-max along the row and column dimensions.
    width_x = fwhm(np.arange(0,len(col))-max_idx[1], col)
    width_y = fwhm(np.arange(0,len(row))-max_idx[0], row)

    # Guess at rotation angle.
    if sig_stixels is not None:
        if sig_stixels.shape[0] > 1:
            xy_centered = sig_stixels.astype(float)
            xy_centered[:,0] -= np.mean(xy_centered[:,0])
            xy_centered[:,1] -= np.mean(xy_centered[:,1])
            # Fit a line to get the slope.
            m, _ = np.polyfit(xy_centered[:,0],xy_centered[:,1],1)
            theta = np.arctan(m)/pi*180 # Convert to degrees
        else:
            theta = 0.0
    else:
        theta = 0.0
    return amplitude, theta, max_idx[1], max_idx[0], width_x, width_y

# Define a 2D Gaussian function
def gaussian2d(xy, amplitude, theta, xo, yo, sigma_x, sigma_y):
    (y,x) = xy # Weird python convention to reverse x,y np.indices
    xo = float(xo)
    yo = float(yo)    
    a = (cosd(theta)**2)/(2*sigma_x**2) + (sind(theta)**2)/(2*sigma_y**2)
    b = -(sind(2*theta))/(4*sigma_x**2) + (sind(2*theta))/(4*sigma_y**2)
    c = (sind(theta)**2)/(2*sigma_x**2) + (cosd(theta)**2)/(2*sigma_y**2)
    g = amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2)))
    return g

# Get an initial estimate of the Gabor parameters for the fit.
def gabor_moments(data, sig_stixels = None):
    # Get the Gaussian moments.
    amplitude, theta, x0, y0, width_x, width_y = gaussian_moments(data, sig_stixels)
    sigma = (width_x + width_y) / 2.0
    wavelength = 1.5
    phase = 0.0
    return amplitude, theta, x0, y0, sigma, wavelength, phase

# Gabor function.
def gabor(xy, amplitude, theta, x0, y0, sigma, wavelength, phase): #gabor(xy, *params):
    # amplitude, theta, x0, y0, sigma, gamma, wavelength, phase = params
    (y,x) = xy # Weird python convention to reverse x,y np.indices
    x0 = float(x0)
    y0 = float(y0)
    g = amplitude * np.exp(-((x-x0)**2 + (y-y0)**2)/(2*sigma**2))
    c = np.cos((cosd(theta) * (x-x0) + sind(theta) * (y-y0)) / wavelength + phase)
    # c = np.cos((cosd(theta) * 2*pi*(x-x0)/np.min(x.shape) + sind(theta) * 2*pi*(y-y0)/np.min(x.shape)) * 7 + phase)
    s = g * c # Gabor is the product of the Gaussian and the cosine
    if s.any() == None:
        s = np.zeros(s.shape)
    return s

def fit_gabor(M, *args):
    x, y = M
    arr = np.zeros(x.shape)
    arr = gabor(M, *args)
    return arr.ravel()

def fit_gaussian(M, *args):
    x, y = M
    arr = np.zeros(x.shape)
    arr = gaussian2d(M, *args)
    return arr.ravel()

# Compute the goodness of fit
def compute_gof(xdata, ydata, func, *popt):
    residuals = ydata - func(xdata, *popt)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum( (ydata - np.mean(ydata))**2 )
    r_squared = 1 - (ss_res / ss_tot)
    return r_squared, ss_res # Return r^2 and mse

# Make sure the start parameters fall between the upper and lower bounds.
def check_start_params(params, non_bds):
    for j in range(len(params)):
        if params[j] < non_bds[0][j]: # Check lower bound
            params[j] = non_bds[0][j]
        elif params[j] > non_bds[1][j]:
            params[j] = non_bds[1][j]
    return params

# amplitude, theta, xo, yo, sigma_x, sigma_y
def fit_gaussian_rf(s_map, theta_guess = 0.0):
    non_bds = ([-1.0, -180.0, 1.0, 1.0, 0.25, 0.25], 
        [1.0, 180.0, np.max(s_map.shape).astype('float'), np.max(s_map.shape).astype('float'), np.max(s_map.shape).astype('float')/3.0, np.max(s_map.shape).astype('float')/3.0])
    # Normalize the RF.
    s_map /= np.max(np.abs(s_map))
    params = gaussian_moments(s_map) # Initial parameter estimates
    params = list(params)
    params[1] = theta_guess
    # Check the start parameters.
    params = check_start_params(params, non_bds)
    # popt, _ = curve_fit(f=linear_filter_function, xdata=t, ydata=time_course[:,0].T, p0=p0, bounds=non_bds, max_nfev=1e05)
    p_gauss, _ = curve_fit(f=fit_gaussian, xdata=np.indices(s_map.shape), ydata=s_map.ravel(), p0=params, bounds=non_bds, max_nfev=1e05)
    r_2,mse = compute_gof(np.indices(s_map.shape),s_map,gaussian2d,*p_gauss)
    return p_gauss, r_2

# amplitude, theta, x0, y0, sigma, wavelength, phase
def fit_gabor_rf(s_map, theta_guess = 0.0):
    non_bds = ([-1.0, -180.0, 1.0, 1.0, 0.25, 1.5, -np.pi], 
        [1.0, 180.0, np.max(s_map.shape).astype('float'), np.max(s_map.shape).astype('float'), np.max(s_map.shape).astype('float')/3.0, np.max(s_map.shape).astype('float')/5.0,np.pi])
    # Normalize the RF.
    s_map /= np.max(np.abs(s_map))
    # Gabor fit
    params = gabor_moments(s_map) # Initial parameter estimates
    params = list(params)
    params[1] = theta_guess
    # Check the start parameters.
    params = check_start_params(params, non_bds)
    # p_gabor, _ = curve_fit(fit_gabor, np.indices(s_map.shape), s_map.ravel(), params)
    p_gabor,_ = curve_fit(f=fit_gabor, xdata=np.indices(s_map.shape), ydata=s_map.ravel(), p0=params, bounds=non_bds, max_nfev=1e05)
    # Compute the goodness-of-fit
    r_2,mse = compute_gof(np.indices(s_map.shape),s_map,gabor,*p_gabor)
    return p_gabor, r_2

# y,x = np.indices(s_map.shape)
# xdata = np.vstack((X1.ravel(), Y1.ravel())).astype(float)
# popt, pcov = curve_fit(fit_gabor, xdata, s_map.ravel(), params)
# popt, pcov = curve_fit(fit_gabor, np.indices(s_map.shape), s_map.ravel(), params)


# def fit_gaussian_rf(s_map):
#     x = np.arange(0,s_map.shape[0])+1
#     y = np.arange(0,s_map.shape[1])+1
#     X, Y = np.meshgrid(x, y)
#     # Ravel the meshgrids of X, Y points to a pair of 1-D arrays.
#     xdata = np.vstack((X.ravel(), Y.ravel()))

#     # Normalize.
#     s_map /= np.max(np.abs(s_map))

#     # Get the initial guess.
#     max_idx = np.unravel_index(np.argmax(np.abs(s_map), axis=None), s_map.shape)
#     p0 = [s_map[max_idx], max_idx[1], max_idx[0], 5, 5, 0]
#     p0 = [s_map[max_idx],0, max_idx[1], max_idx[0], 5, 5]

#     bounds = ((-1.0, 1.0, 1.0, 0.5, 0.5, 0),(1,78.0,58.0,30,30,3.14159))

#     # Fit
#     popt, pcov = curve_fit(_gaussian2d, xdata, s_map.ravel(), p0, bounds=bounds)
#     return None

class Analysis():
    # def __init__(self, 
    #         experimentName: str):
    #     self.experimentName = experimentName
    #     self.M = Metadata(experimentName)
    #     self.D = Dataset(experimentName)
    #     self.S = Stimulus()
    
    # Cosine of argument in degrees
    def cosd(self, x):
        return np.cos( np.deg2rad(x) )
    
    # Sine of argument in degrees
    def sind(self, x):
        return np.sin( np.deg2rad(x) )

    # Define a 2D Gaussian function
    def gaussian2d(self, xy, amplitude, xo, yo, sigma_x, sigma_y, theta):
        (x, y) = xy
        xo = float(xo)
        yo = float(yo)    
        a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
        b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
        c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
        g = amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                                + c*((y-yo)**2)))
        return g.ravel()
    
    def _gaussian2d(self, xy, *args):
        return self.gaussian2d(xy,*args)

    def fit_gaussian_rf(self, s_map):
        x = np.arange(0,s_map.shape[0])+1
        y = np.arange(0,s_map.shape[1])+1
        X, Y = np.meshgrid(x, y)
        # Ravel the meshgrids of X, Y points to a pair of 1-D arrays.
        xdata = np.vstack((X.ravel(), Y.ravel()))
        # Get the initial guess.
        max_idx = np.unravel_index(np.argmax(s_map, axis=None), s_map.shape)
        p0 = [1, max_idx[1], max_idx[0], 5, 5, 0]

        # Fit
        popt, pcov = curve_fit(self._gaussian2d, xdata, s_map.ravel(), p0)
        return None

    def twoD_Gaussian(xy, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
        (x, y) = xy
        xo = float(xo)
        yo = float(yo)    
        a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
        b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
        c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
        g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                                + c*((y-yo)**2)))
        return g.ravel()

    def gaussian(self, height, center_x, center_y, width_x, width_y, rotation):
        """Returns a gaussian function with the given parameters"""
        width_x = float(width_x)
        width_y = float(width_y)

        rotation = np.deg2rad(rotation)
        center_x = center_x * np.cos(rotation) - center_y * np.sin(rotation)
        center_y = center_x * np.sin(rotation) + center_y * np.cos(rotation)

        def rotgauss(x,y):
            xp = x * np.cos(rotation) - y * np.sin(rotation)
            yp = x * np.sin(rotation) + y * np.cos(rotation)
            g = height*np.exp(
                -(((center_x-xp)/width_x)**2+
                  ((center_y-yp)/width_y)**2)/2.)
            return g
        return rotgauss

    def moments(self, data):
        """Returns (height, x, y, width_x, width_y)
        the gaussian parameters of a 2D distribution by calculating its
        moments """
        total = data.sum()
        X, Y = np.indices(data.shape)
        x = (X*data).sum()/total
        y = (Y*data).sum()/total
        col = data[:, int(y)]
        width_x = np.sqrt(abs((np.arange(col.size)-y)**2*col).sum()/col.sum())
        row = data[int(x), :]
        width_y = np.sqrt(abs((np.arange(row.size)-x)**2*row).sum()/row.sum())
        height = data.max()
        return height, x, y, width_x, width_y, 0.0

    
    def fitgaussian(self, data):
        """Returns (height, x, y, width_x, width_y)
        the gaussian parameters of a 2D distribution found by a fit"""
        params = self.moments(data)
        errorfunction = lambda p: np.ravel(self.gaussian(*p)(*np.indices(data.shape)) - data)
        p, success = scipy.optimize.leastsq(errorfunction, params)
        return p, success

