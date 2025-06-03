"""
Module for fitting functions to the data, both estim and vstim.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import scipy as sp
import config as cfg
from scipy.optimize import curve_fit,minimize
import math
from scipy import stats
from scipy.stats import median_absolute_deviation as mad
import lnp
import warnings
from scipy.optimize import OptimizeWarning
import pdb

def get_abc_params_gaussian2d(theta,sigma_x,sigma_y):
    """Computes a, b, and c parameters for 2D Gaussian.
    See https://en.wikipedia.org/wiki/Gaussian_function
    @param theta: Angle of the data to fit (radians).
    @param sigma_x: Standard deviation in the x dimension.
    @param sigma_y: Standard deviation in the y dimension.
    @return a,b,c: Tuple with three param values.
    """
    a = ((math.cos(theta)**2 / (2 * (sigma_x**2))) +
        (math.sin(theta)**2 / (2 * (sigma_y**2))))

    b = (-(math.sin(2 * theta) / (4 * (sigma_x**2))) +
        (math.sin(2 * theta) / (4 * (sigma_y**2))))

    c = ((math.sin(theta)**2 / (2 * (sigma_x**2))) +
        (math.cos(theta)**2 / (2 * (sigma_y**2))))

    return a,b,c

def gaussian2d(xy,A,mu_x,mu_y,sigma_x,sigma_y,theta):
    """
    2D Gaussian function
    TODO: document
    """
    x,y = xy
    a,b,c = get_abc_params_gaussian2d(theta,sigma_x,sigma_y)
    
    return A * np.exp(-((a * ((x - mu_x)**2)) +
                     (2 * b * (x - mu_x) * (y - mu_y)) +
                     (c * ((y - mu_y)**2))))

def init_params_gaussian2d(spatial_map,sig_stixels):
    """
    Gets initial guess of 2D Gaussian parameters for fitting STA.
    Sets the initial guess of theta to be 0.
    Parameters:
        spatial_map: spatial map of cell of size y,x,color
        sig_stixels: array of sig stixels from full STA
    Returns:
        list of parameters: A, mu_x,mu_y,sigma_x,sigma_y,theta
    """
    
    # Mask and grayscale the mean subtracted map
    masked_map = np.mean(
                    lnp.mask_sta(
                        spatial_map[None,...] - np.mean(spatial_map),
                        sig_stixels,norm=False
                    ),
                    axis=-1,keepdims=True
                 )

    # Compute stats on the map to get initial guesses.
    mu_x = np.median(sig_stixels[:,1])
    mu_y = np.median(sig_stixels[:,0])
    sigma_x = mad(sig_stixels[:,0])
    sigma_y = mad(sig_stixels[:,1])

    # Safeguard against small initial stds.
    if sigma_x < 1:
        sigma_x = 2
    
    if sigma_y < 1:
        sigma_y = 2

    A = np.max(np.abs(spatial_map))

    # Flip depending on the polarity.
    if np.abs(np.max(masked_map)) < np.abs(np.min(masked_map)):
        A *= -1

    return [A,mu_x,mu_y,sigma_x,sigma_y,0]

def fit_gaussian2d(spatial_map,sig_stixels,grayscale=True,
                   channel=None,blur=False):
    """
    Fits a 2D Gaussian to the STA spatial map. By default, the map 
    is grayscaled by averaging over the channels. If channel is
    indicated, then the Gaussian is fitted over that single channel.
    Note: multiple channel fitting is not supported in this function.
    
    Parameters: 
        spatial_map: spatial map of STA, of size y,x,color
        sig_stixels: array of sig stixels
        grayscale: boolean indicating greyscale
        channel: channel over which to fit; grayscale must be set to False
        blur: whether to filter the STA (not implemented yet)
    Returns:
        the optimal parameters from fitting
    """

    # Depending on the channels/grayscale arguments, index accordingly.
    if (not grayscale and channel is None) or\
       (grayscale and channel is not None):
        assert False,"incompatible args for grayscale and channels."

    if grayscale:
        map_to_fit = np.mean(spatial_map,axis=-1,keepdims=True)
    else:
        map_to_fit = spatial_map[...,channel][...,None]

    # Get initial guess of params for the map. 
    initial_params = init_params_gaussian2d(map_to_fit,sig_stixels)
    x = np.arange(0,spatial_map.shape[1],1)
    y = np.arange(0,spatial_map.shape[0],1)
    xx,yy = np.meshgrid(x,y)

    # Mean subtract the map and fit.
    map_to_fit -= np.mean(map_to_fit)

    try:
        popt,_ = curve_fit(
                    gaussian2d,(xx.ravel(),yy.ravel()),
                    map_to_fit.ravel(),
                    p0=initial_params,maxfev=10000
                )
    except:
        print('could not fit Gaussian')
        return None
    #pdb.set_trace()
    
    return popt

def diff_casc_lfps(t: int or np.ndarray , p1: float, p2: float,
           tau1: float, tau2: float, n1: int or float,
           n2: int or float) -> int or np.ndarray:
    '''
    Difference of cascades of low-pass filters to fit STA time course with,
    as in Chichilnisky and Kalmar, 2002 (J Neuro)
    Parameters:
        t: Frame of the time course (sample number)
        p1: Proportional to one of the peaks
        p2: Proportional to one of the peaks
        tau1: Time constant to one of the peaks
        tau2: Time constant to one of the peaks
        n1: Width of one of the peaks
        n2: Width of one of the peaks
    Returns: Function evaluated at some t and the above parameters.
    '''
    tf1 = p1 * ((t / tau1)**n1) * np.exp(-n1 * ((t / tau1) - 1))
    tf2 = p2 * ((t / tau2)**n2) * np.exp(-n2 * ((t / tau2) - 1))

    return tf1 - tf2

def init_params_diff_casc_lfps(timecourse,n1_init=15,n2_init=6,
                      p1_scalar = 1.05,p2_scalar = .95,
                      tau1_scalar=1.1,tau2_scalar = 1.25,
                      min_trough_peak_ratio = .15):
    """
    Initialize the parameters for fitting a time course. Originally from 
    heuristics in a MATLAB implementation written by Greg Field.

    Parameters:
        timecourse: time reversed STA time course
        ... a list of default parameters.
    Returns:
        array of parameters.
    """

    # Initialize parameters, based on time course.
    n1 = n1_init
    n2 = n2_init
    max_tc = np.max(timecourse)
    max_frame = np.argmax(timecourse)
    min_tc = np.min(timecourse)
    min_frame = np.argmin(timecourse)

    # Set things up, depending on the polarity of the cell.
    if min_frame < max_frame:
        peak_frame = min_frame + 1
        trough_frame = max_frame + 1
        peak_scale = min_tc
        trough_scale = max_tc
    else:
        peak_frame = max_frame + 1
        trough_frame = min_frame + 1
        peak_scale = max_tc
        trough_scale = min_tc

    p1 = peak_scale * p1_scalar
    p2 = trough_scale * p2_scalar
    tau1 = peak_frame * tau1_scalar

    if np.abs(trough_scale / peak_scale) > min_trough_peak_ratio:
        tau2 = trough_frame
    else:
        tau2 = trough_frame * tau2_scalar

    return [p1, p2, tau1, tau2, n1, n2]

def check_init_params_diff_casc_lfps(initial_params):
    """
    Checks the initial parameters for fitting a differnece of cascades 
    of low pass filters to the time course. MODIFIES the array!

    Puts bounds on the parameters to avoid crazy fits. To compute these values,
    I looked at the distributions of popt's found across roughly 10000 cells
    and took the median and MAD STD to come up with some ad hoc boundaries.
    I allow tau parameters to just span the whole positive number set because
    those tend to not be an issue. This is very ad hoc, but it's reasonable.

    Parameters:
        initial_params: initial parameters
    Returns:
        updated list of params from bound checking.
    """
    lower_bounds = cfg.LOWER_BOUNDS_LFP_FIT
    upper_bounds = cfg.UPPER_BOUNDS_LFP_FIT

    for i in range(len(initial_params)):
        if initial_params[i] < lower_bounds[i]:
            initial_params[i] = lower_bounds[i]

        if initial_params[i] > upper_bounds[i]:
            initial_params[i] = upper_bounds[i]
    
    return initial_params

def fit_sta_timecourse(timecourse, max_iter=100000,
                       check_init_params=True,bound_opt=True):
    '''
    Fit an STA time course with a difference of cascades of low pass filters.
    Only fits a single channel, so pass in the color channel of interest outside
    this function.

    Parameters:
        timecourse: Timecourse vector (array), time reversed (so that the 
        time of spike is the 0th index, NOT the -1.
        max_iter: max numnber of iterations in the optimization
        check_init_params: boolean indicating whether to check params bounds
        bound_opt: boolean indicating whether to bound the optimization 
    Returns:
        The fit to the function, evaluated at the data points from ttf
        Optimal parameters found.

    Initialization of the parameters is critical for a good fit. Here, the 
    parameters are initialized in a similar way to a MATLAB implementation of this 
    function, written by GDF. This Python
    version was based on code originally written by bhaishahster.
    '''

    # Get initial parameters and check the bounds if indicated.
    initial_params = init_params_diff_casc_lfps(timecourse) 

    if check_init_params:
        initial_params = check_init_params_diff_casc_lfps(initial_params)

    if bound_opt:
        bounds = (cfg.LOWER_BOUNDS_LFP_FIT,cfg.UPPER_BOUNDS_LFP_FIT)
    else:
        bounds = ()

    # Get the fit and evaluate at the popt
    try:
        popt, _ = curve_fit(
                    diff_casc_lfps, np.arange(1,timecourse.shape[0]+1),
                    timecourse, p0=initial_params, maxfev=max_iter,
                    bounds=bounds
                  )
    except:
        print("couldn't fit time course")
        return None,None

    timecourse_fit = diff_casc_lfps(np.arange(1,timecourse.shape[0]+1), *popt)

    return popt,timecourse_fit

def lfp_anon(p1, p2, tau1, tau2, n1, n2):
    '''
    Difference of cascades of low-pass filters to fit STA time course with,
    as in Chichilnisky and Kalmar, 2002 (J Neuro). This is the anonymous 
    callable object for use in computing zero cross, etc.
    Parameters:
        t: Frame of the time course (sample number)
        p1: Proportional to one of the peaks
        p2: Proportional to one of the peaks
        tau1: Time constant to one of the peaks
        tau2: Time constant to one of the peaks
        n1: Width of one of the peaks
        n2: Width of one of the peaks
    Returns: Function evaluated at some t and the above parameters.
    '''
    return lambda t: p1 * ((t / tau1)**n1) * np.exp(- n1 *((t / tau1) - 1)) - \
        p2 * ((t / tau2)**n2) * np.exp(- n2 * ((t / tau2) - 1))


def gss(f,a,b,type_search,tol=math.sqrt(np.finfo(np.float).eps)):
    '''
    Implementation of the Golden Section Search routine to numerically solve for
    the extrema of a smooth function. Assumes within a time window, there exists
    only a single extrema. Algorithm will fail otherwise. See
    https://en.wikipedia.org/wiki/Golden-section_search and Numerical Recipes in C
    book for more details; this implemetation is based largely off these two.
    @param f Callable function object.
    @param a Lower bound of search
    @param b Upper bound of search
    @param type_search Type of search to be used. This implementation can compute
        either the minimum (type_search = 'min') or the maximum
        (type_search = 'max') of a function.
    @param tol Tolerance parameter of the difference in the two bounds. Default is
        square root of the machine's floating point precision (reccomendation
        from Numerical Recipes in C)
    @return The estimated value at which the extrema occurs, and the value of
            extrema.
    '''

    # Compute Golden ratio.
    gr = (math.sqrt(5) + 1) / 2

    # Find midpoints.
    c = b - (b - a) / gr
    d = a + (b - a) / gr

    # Iteratively search for extrema.
    while abs(c - d) > tol:

        if (f(c) < f(d) and type_search == 'min') or\
           (f(c) > f(d) and type_search == 'max'):
            b = d
        else:
            a = c

        # Recompute midpoints.
        c = b - (b - a) / gr
        d = a + (b - a) / gr

    # Once all checks out, return the mean and the function evaluated there.
    return (b + a) / 2,f((b + a) / 2)

def bissection_search_root(f,a,b,max_iter=1000):
    '''
    Implementation of the Bissection search root finding method to numerically solve
    for when a time course function equals zero ('time of zero crossing'). See
    https://en.wikipedia.org/wiki/Bisection_method and Numerical Recipes in C for
    further details. Implemented iteratively, until MAX_ITER is reached. Tolerance
    parameter used is for the difference in the bound. Default is the
    floating point precision of the machine multiplied by (a + b) / 2.
    @param f Callable function object
    @param a Lower bound on search
    @param b Upper bound on search
    @param max_iter maximum number of iterations
    @return Time at which function evaluates to zero, and the found value (close to 0)
    @raises OptimizeWarning if max_iter is exceeded, returns None,None.
    '''

    # Tolerance parameter.
    tol = np.finfo(np.float).eps * (a + b) / 2

    iter_cnt = 1 # Counter.

    while True:
        c = (a + b) / 2

        # Test convergence, otherwise update bounds.
        if f(c) == 0 or abs((b - a) / 2) < tol:
            #pdb.set_trace()
            return c,f(c)
        iter_cnt +=1

        if np.sign(f(a)) == np.sign(f(c)): a = c
        else: b = c

        if iter_cnt == max_iter:
            msg_out = 'Maximum number of iterations exceeded. Root not found.'
            warnings.warn(OptimizeWarning(msg_out))
            return None,None

def nl(x,alpha,beta,gamma):
    """ 
    Functional form of the nonlinearity in the LNP model, of the form
        f(x) = alpha * L(beta * x + gamma)
    Parameters:
        x: generator signal, either single scalar or array
        alpha, beta, gamma: free parameters
    Returns:
        output of the function
    """
    return alpha * stats.logistic.cdf((beta * x) + gamma)

def fit_nonlinearity(mean_generator,mean_spike_counts,abs_max_g,n_eval=100,
                     p0=[1, 1, 1],maxfev=10000):
    """
    Function to fit a nonlinearity with a parameterized function.
    Parameters:
        mean_generator: array of mean generator values
        mean_spike_counts: array of mean spike_counts
        abs_max_g: maximum value of generator signal to evaluate
        n_eval: number of data points to evaluate nonlinearity at
        p0: initial guess of popt
        maxfev: max num. iterations
    Returns:
        3 tuple of evaluated generator values, firing rate and the popt
    """

    # Create bounds on optimization: alpha is max FR, which is > 0 and < 200 (Hz).
    lower_bounds = [0,-np.inf,-np.inf]
    upper_bounds = [200,np.inf,np.inf]

    try:
        popt,_ = curve_fit(nl,mean_generator,mean_spike_counts,p0=p0,
                            maxfev=maxfev,bounds=(lower_bounds,upper_bounds))
        mean_generator_fit = np.linspace(-abs_max_g,abs_max_g,n_eval)
        mean_spike_counts_fit = nl(mean_generator_fit,*popt)
    except:
        print("could not fit nonlinearity")
        return None,None,None

    return mean_generator_fit,mean_spike_counts_fit,popt

def clean_probs(probs,amplitudes,thresh_lb=0.5,
                thresh_ub=3.8,flip_thr=0.25,
                zero_thr=0.15):
    """
    Utility to clean up probabilities from g-sort. This deals with noisy 
    amplitudes and helps with the 0-1 disambiguation problem. Adapted from
    Praful Vasireddy code.

    The goal here is to find the amplitude that is closest to threshold 
    (~.5 probability), flip 0s above this to 1, and 1s below this to 0. This
    removes noisy data points and cleans up sigmoid fitting. Uses a bunch of 
    heuristics.

    Parameters:
        probs: array of probabilities
        amplitudes: array of current amplitudes
        thresh_lb: lower bound of amplitudes
        thresh_ub: upper bound of amplitudes
        flip_thr: min probability required for considering cleaning
        zero_thr: threshold probability, below which it's basically 0.
    Returns:
        cleaned probabilities
    """

    thresh_inds = np.where(
                 (amplitudes <= thresh_ub) &
                 (amplitudes >= thresh_lb)
                 )[0]

    if np.amax(probs[thresh_inds]) >= flip_thr:
        zero_inds = np.where(probs <= zero_thr)[0]
        '''
        Find the amplitude with probability closest to 0.5. Flip all 0s above 
        this amplitude to 1s, force all 0s below this amp to exactly 0
        '''
        thr_idx = thresh_inds[np.argmin(np.absolute(probs[thresh_inds] - 0.5))]
        zero_inds_1 = zero_inds[zero_inds > thr_idx]
        zero_inds_0 = zero_inds[zero_inds < thr_idx]

        probs_cleaned = probs.copy()
        probs_cleaned[zero_inds_1] = 1
        probs_cleaned[zero_inds_0] = 0

        return probs_cleaned

    return probs

def sigmoid(x,k,x0):
    """
    Functional form of sigmoid to fit activation curves with.
    Parameters:
        x: input data
        k: free param
        x0: free param
    Returns:
        sigmoid evaluated at point(s) x and free params
    """
    return 1 / (1 + np.exp(-k * (x - x0)))

def fit_sigmoid_mle(x,y,init_params=[1,1,1],method='Nelder-Mead'):
    """
    Fits sigmoid to data using MLE. Defines an inner function for 
    negative log likelihood computation.
    Parameters:
        x: input data
        y: output data
        init_params: initial params for log likelihood function
        method: minimization method; Nelder-Mead works very well
    Returns:
        2-tuple of popt and the sigmoid fit from minimization.
    """

    def ll_sigmoid(params):
        """
        Log likelihood function for sigmoid fitting.
        Parameters:
            params: a tuple of the following:
                k: free param
                x0: free param
                sd: std deviation for the logistic pdf.
        Returns:
            negative log likelihood
        """
        k,x0,sd = params
        y_pred = sigmoid(x,k,x0)

        return -np.sum(stats.norm.logpdf(y, loc=y_pred, scale=sd))

    results = minimize(ll_sigmoid, init_params, method=method)
    popt = results.x

    # Only care about the first two popt for sigmoid fitting.
    fit = sigmoid(x,*popt[0:2])

    return popt[0:2],fit