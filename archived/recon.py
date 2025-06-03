"""
Module of functions for doing reconstruction analysis, including
greedy stimulation algorithm simulation.
"""
import os
os.environ["OMP_NUM_THREADS"] = "30"
import scipy as sp
import numpy as np
import cv2
import src.lnp as lnp
import src.config as cfg
import src.fitting as fit
import src.io_util as io
import src.fitting as fit
from numpy import linalg
import pdb
import torch
import cvxpy as cp

def get_random_stimulus(num_stimuli,stixel_size,
                        n_pixels_x,n_pixels_y,
                        seed=None,binary=True):
    """
    Generates random white noise stimuli tensors.
    Parameters:
        num_stimuli: number of stimuli
        stixel_size: pixels per stixel on the stimulus
        n_pixels_x: number of x pixels
        n_pixels_y: number of y pixels
        seed: random number gen seed
        binary: binary noise (only supported for now)
    Returns:
        Stimulus tensor of size num_stimuli, n_pixels_y, n_pixels_x
    """

    if not binary: raise ValueError('Non-binary stimuli not yet implemented.')

    # Grab the number of stixels from the stixel size.
    num_stixels_y = int(n_pixels_y / stixel_size)
    num_stixels_x = int(n_pixels_x / stixel_size)

    # Randomly generate frames.
    if seed is not None: sp.random.seed(seed)
    random_stimuli_tmp = sp.random.randn((num_stixels_y * num_stixels_x),
                                      num_stimuli) > 0
    random_stimuli = np.empty(np.shape(random_stimuli_tmp))

    for i in range(np.shape(random_stimuli_tmp)[0]):
        for j in range(np.shape(random_stimuli_tmp)[1]):
            random_stimuli[i,j] = int(random_stimuli_tmp[i,j])
    random_stimuli -= 0.5

    return random_stimuli * (0.48 / 0.5)

def get_random_stimuli_jitter(n_stimuli,stixel_size,x_dim,y_dim,
                              n_pixels_x=640,n_pixels_y=320,seed=None,
                              factor=2):
    """
    Adds jitter to the random stimuli; useful for training linear 
    reconstruction filters. Calls the get_random_stimulus function under the 
    hood and applies a jitter.
    Parameters:
        n_stimuli: number of stimuli
        stixel_size: pixels per stixel
        x_dim: size of x dim
        y_dim: size of y dim
        n_pixels_x: number of x pixels
        n_pixels_y: number of y pixels
        seed: rng seed
        factor: amount of extra stimulus to get for jitter
    Return:
        Stimulus tensor of size num_stimuli, n_pixels_y, n_pixels_x
    """
    stimulus = get_random_stimulus(n_stimuli,stixel_size,
                                   n_pixels_x=int(n_pixels_x * factor),
                                   n_pixels_y=int(n_pixels_y * factor),
                                   binary=True,seed=seed)
    stimulus_reshape = []

    # Compute the field size and resample.
    field_x = int((n_pixels_x * factor) / stixel_size)
    field_y = int((n_pixels_y * factor) / stixel_size)

    for i in range(n_stimuli):
        tmp = cv2.resize(np.reshape(stimulus[:,i],(field_y,field_x)),
                                  dsize=(int(x_dim * factor),
                                         int(y_dim * factor)),
                                  interpolation=cv2.INTER_NEAREST)
        stimulus_reshape.append(tmp)

    # Now, for each stimulus, get a random window to simulate "jitter".
    stimulus_jitter = []
    np.random.seed(seed)

    for i in range(n_stimuli):
        random_x = np.random.randint(0,(int(x_dim * factor)))
        random_y = np.random.randint(0,(int(y_dim * factor)))

        if random_x + x_dim < x_dim * factor:
            x_window = np.arange(random_x,random_x + x_dim)
        else:
            x_window = np.arange(random_x - x_dim,random_x)

        if random_y + y_dim < y_dim * factor:
            y_window = np.arange(random_y,random_y + y_dim)
        else:
            y_window = np.arange(random_y - y_dim,random_y)

        # Chop the stimulus according to the new window
        tmp = stimulus_reshape[i][y_window,:][:,x_window]
        stimulus_jitter.append(tmp)

    return np.moveaxis(np.c_[[tmp for tmp in stimulus_jitter]],0,-1)

def get_gaussian_filters(gaussian_popt,original_stixel):
    """
    Evaluates 2D Gaussian on a stixel 4 grid. Assumes that the fits were 
    obtained on the spatial maps normalized in 0-1 range. If indicated,
    clips anything below (abs) 1/255 to zero.

    Parameters:
        gaussian_popt: list of Gaussian params for each cell
        original_stixel: original stixel size. Used to to rescale params.
    Returns:
        list of evaluated Gaussian maps
    """

    # Get the scale based on the original stixel size.
    scale = original_stixel / cfg.RESAMPLE_STIX

    # For each cell, evaluate the Gaussian fit at stixel size 4.
    x = np.arange(0,cfg.RESAMPLE_FIELD_X)
    y = np.arange(0,cfg.RESAMPLE_FIELD_Y)
    xx,yy = np.meshgrid(x,y)
    gaussian_filters = []

    for cc in range(len(gaussian_popt)):

        if gaussian_popt[cc] is None:
            gaussian_filters.append([])
            continue

        # Scale the mean and std. 
        popt = np.asarray(gaussian_popt[cc].copy())
        popt[1:5] *= scale
        gaussian_filter = fit.gaussian2d((xx,yy),*popt)
        gaussian_filter[np.abs(gaussian_filter) < 1/255] = 0
        gaussian_filters.append(gaussian_filter)
    
    return gaussian_filters

def learn_recon_filters(stixel,encoder):
    """
    Learns reconstruction filters using least squares regression

    TODO: document
    """

    n_cells,n_pixels_y,n_pixels_x = encoder.shape
    gaussian_filters_matrix = np.reshape(
                         encoder,
                        (n_cells,n_pixels_x * n_pixels_y)
                        )
    stimulus = get_random_stimuli_jitter(cfg.N_TRAIN_STIMULI,stixel,
                                n_pixels_x,n_pixels_y,seed=cfg.TRAIN_SEED)
    stimulus_vector = np.reshape(stimulus,(n_pixels_x * n_pixels_y,
                                cfg.N_TRAIN_STIMULI))
    R = np.transpose(np.maximum(gaussian_filters_matrix@stimulus_vector,0))
    S = stimulus_vector.T
    n_pixels = int(n_pixels_x * n_pixels_y)
    n_cells = gaussian_filters_matrix.shape[0]
    W = np.zeros((n_cells,n_pixels))

    for i in range(n_pixels):

        # Grab cells with non-zero filters.
        cell_inds_sta = np.argwhere(np.abs(
                            gaussian_filters_matrix[:,i]) > 0
                        ).flatten()

        if cell_inds_sta.shape[0] == 0:
            continue

        R_t = R[:,cell_inds_sta]

        S = stimulus_vector[i,:].T
        W[cell_inds_sta,i] = linalg.inv(R_t.T@R_t)@R_t.T@S

    return W

def learn_recon_filters_ns(encoder):
    """
    TODO: doc
    """
    n_cells,n_pixels = encoder.shape
    gaussian_filters_matrix = np.reshape(
                            encoder,
                            (n_cells,n_pixels)
                            )
    
    # Get stixel size 2 NS tensors and downsample 2x.
    stimulus = io.get_raw_movie(os.path.join(cfg.IMAGENET_PARENT,
                                cfg.IMAGENET_TRAIN),
                                cfg.N_TRAIN_STIMULI
                        )
    stimulus = np.asarray([cv2.resize(stimulus[i,...],
                          dsize=(cfg.RESAMPLE_FIELD_X,cfg.RESAMPLE_FIELD_Y)
                          )
                          for i in range(stimulus.shape[0])
                    ])

    # Normalize the stimulus, center it and take only one color.
    stimulus = stimulus / 255
    stimulus -= np.mean(stimulus)
    stimulus = stimulus[...,0] # doesn't matter since grayscaled
    n_stimuli,n_pixels_y,n_pixels_x = stimulus.shape
    stimulus_vector = np.reshape(stimulus,(n_stimuli,
                                n_pixels_x * n_pixels_y)).T
    R = np.transpose(np.maximum(gaussian_filters_matrix@stimulus_vector,0))
    n_pixels = int(n_pixels_x * n_pixels_y)
    n_cells = gaussian_filters_matrix.shape[0]
    #W = linalg.inv(R.T@R)@R.T@stimulus_vector
    W = np.zeros((n_cells,n_pixels))

    for i in range(n_pixels):

        # Grab cells with non-zero filters.
        cell_inds_sta = np.argwhere(np.abs(
                            gaussian_filters_matrix[:,i]) > 0
                        ).flatten()

        if cell_inds_sta.shape[0] == 0:
            continue

        R_t = R[:,cell_inds_sta]

        S = stimulus_vector[i,:]
        W[cell_inds_sta,i] = linalg.inv(R_t.T@R_t)@R_t.T@S

    return W

def get_encoding_filters(stixel,gaussian_filters_tensor):
    """
    Gets the encoding model filters.

    TODO: document
    """
    n_cells,n_pixels_y,n_pixels_x = gaussian_filters_tensor.shape
    gaussian_filters_matrix = np.reshape(
                    gaussian_filters_tensor,
                    (n_cells,n_pixels_x * n_pixels_y)
                )

    # Get the random stimuli and scale the filters to get desired response.
    stimulus = get_random_stimuli_jitter(
                        cfg.N_TRAIN_STIMULI,stixel,
                        n_pixels_x,n_pixels_y,
                        seed=cfg.PRE_TRAIN_SEED
                    )
    stimulus_vector = np.reshape(
                        stimulus,
                        (n_pixels_x * n_pixels_y,cfg.N_TRAIN_STIMULI)
                    )
    R = np.maximum(gaussian_filters_matrix@stimulus_vector,0)
    mean_firing_rate = np.mean(R.ravel())
    scale = cfg.SPIKES_PER_FLASH / mean_firing_rate

    return gaussian_filters_matrix * scale

def get_encoding_filters_ns(gaussian_filters_tensor):
    """
    Gets encoding model filters 

    TODO: document
    """

    n_cells,n_pixels_y,n_pixels_x = gaussian_filters_tensor.shape
    gaussian_filters_matrix = np.reshape(
                    gaussian_filters_tensor,
                    (n_cells,n_pixels_x * n_pixels_y)
                )

    # Get stixel size 2 NS tensors and downsample 2x.
    stimulus = io.get_raw_movie(os.path.join(cfg.IMAGENET_PARENT,
                                cfg.IMAGENET_PRETRAIN),
                                cfg.N_TRAIN_STIMULI
                        )
    stimulus = np.asarray([cv2.resize(stimulus[i,...],
                          dsize=(cfg.RESAMPLE_FIELD_X,cfg.RESAMPLE_FIELD_Y)
                          )
                          for i in range(stimulus.shape[0])
                    ])

    # Normalize the stimulus, center it and take only one color.
    stimulus = stimulus / 255
    stimulus -= np.mean(stimulus)
    stimulus = stimulus[...,0] # doesn't matter since grayscaled
    stimulus_vector = np.reshape(
                            stimulus,
                            (cfg.N_TRAIN_STIMULI,n_pixels_x * n_pixels_y)
                     ).T
    R = np.maximum(gaussian_filters_matrix@stimulus_vector,0)
    mean_firing_rate = np.mean(R.ravel())
    scale = cfg.SPIKES_PER_FLASH / mean_firing_rate

    return gaussian_filters_matrix * scale

def center_filters(filter_tensor,n_pixels_y,n_pixels_x):
    """
    Centers either the encoding or decoding filters. Takes the sum over cells,
    computes the median x and y pixel and sets it to the ~center of the 
    visual field

    Parameters:
        filter_tensor: tensor of size n_cells,n_pixels_y,n_pixels_x
        n_pixels_y: number of y pixels
        n_pixels_x: number of x pixels
    Returns:
        a centered filter tensor.
    """

    # Take the sum over cells and find nonzero pixels.
    nonzero_stixels = np.argwhere(np.sum(filter_tensor,axis=0) != 0)
    med_y = np.median(nonzero_stixels[:,0])
    med_x = np.median(nonzero_stixels[:,1])

    # Compute the difference from midpoint and roll it.
    diff_y = int(n_pixels_y/2 - med_y)
    diff_x = int(n_pixels_x/2 - med_x)

    return np.roll(filter_tensor,[diff_y,diff_x],axis=(1,2))

def get_filter_dict(cellids,gaussian_popt,original_stixel):
    """
    Function to learn reconstruction filters from LNP model. Simulates
    responses to white noise stimuli using the LNP model and applies least
    squares regression to learn reconstruction filters.

    TODO: DOCUMENT
    """

    # Get the Gaussian maps and get the cells that have nontrivial filters.
    gaussian_filters = get_gaussian_filters(gaussian_popt,original_stixel)
    cells_reconstruction = [cellids[i] for i in range(len(gaussian_filters))
                            if gaussian_filters[i] != []]
    gaussian_filters_tensor = []

    for gaussian_filter in gaussian_filters:
        
        if gaussian_filter == []:
            continue
        gaussian_filters_tensor.append(gaussian_filter)
    gaussian_filters_tensor = np.asarray(gaussian_filters_tensor)
    gaussian_filters_tensor = normalize_filters(gaussian_filters_tensor)

    # Learn encoding and decoding filters for each stixel size. 
    filter_dict = dict()
    filter_dict['encoding_filters'] = dict()
    filter_dict['decoding_filters'] = dict()
    filter_dict['cells_reconstruction'] = cells_reconstruction

    for stixel in cfg.STIXEL_SIZES:

        print(f'learning filters for stixel {stixel}')

        # Scale the Gaussian filters to get the encoder and learn decoder.
        encoder = get_encoding_filters(
                                stixel,
                                gaussian_filters_tensor
                            )
        n_cells,n_pixels_y,n_pixels_x = gaussian_filters_tensor.shape
        decoder = learn_recon_filters(
                                stixel,
                                encoder.reshape((n_cells,
                                                 n_pixels_y,
                                                 n_pixels_x)),
                            )
        filter_dict['encoding_filters']["%s"%str(stixel)] = encoder
        filter_dict['decoding_filters']["%s"%str(stixel)] = decoder

    return filter_dict

def normalize_filters(gaussian_filters):
    """
    Mean subtracts and L2-normalizes the Gaussian filters among nonzero 
    elements.
    TODO: DOCUMENT
    """
    n_cells,n_pixels_y,n_pixels_x = gaussian_filters.shape
    gaussian_filters = np.reshape(gaussian_filters,(n_cells,
                                  n_pixels_x * n_pixels_y))
    normalized_filters = []

    for cc in range(gaussian_filters.shape[0]):
        sig_stix = np.argwhere(gaussian_filters[cc,:] != 0).flatten()
        normalized_filter = gaussian_filters[cc,:].copy()

        # TODO: figure out if right thing to do.
        #normalized_filter[sig_stix] -= np.mean(normalized_filter[sig_stix])
        normalized_filter[sig_stix] /= linalg.norm(normalized_filter[sig_stix])
        normalized_filters.append(normalized_filter)
    
    normalized_filters = np.asarray(normalized_filters)

    return normalized_filters.reshape(n_cells,n_pixels_y,n_pixels_x)

def get_responses_lnp(gaussian_filters_tensor,stimulus,nl_popt):
    """
    Gets the responses from the LNP model.
    TODO: document
    """

    # Vectorize stimuli and filters 
    n_cells,n_pixels_y,n_pixels_x = gaussian_filters_tensor.shape
    gaussian_filters = np.reshape(
                            gaussian_filters_tensor,
                            (n_cells,n_pixels_y * n_pixels_x)
                        )
    n_stimuli = stimulus.shape[0]
    stimulus_vector = np.reshape(stimulus,(n_stimuli,n_pixels_x * n_pixels_y))

    # Compute generator signals and look up the FR from nl.
    generators = gaussian_filters@stimulus_vector.T
    responses = []

    for cc in range(len(nl_popt)):
        responses.append(fit.nl(generators[cc],*nl_popt[cc]))
    
    return np.asarray(responses) * (cfg.N_MS_POST_FLASH / 1000)

def get_filter_dict_ns(cellids,gaussian_popt,nl_popt,original_stixel_size):
    """
    Learns reconstruction filters for natural scenes using the LNP model.

    TODO: DOCUMENT
    """

    # Get the Gaussian filters and normalize them.
    gaussian_filters = get_gaussian_filters(gaussian_popt,original_stixel_size)
    '''
    cells_reconstruction = [cellids[i] for i in range(len(gaussian_filters))
                            if gaussian_filters[i] != [] and nl_popt[i]
                            is not None]
    '''
    cells_reconstruction = [cellids[i] for i in range(len(gaussian_filters))
                            if gaussian_filters[i] != []]
    nl_popt = [i for i in nl_popt if i is not None]
    gaussian_filters_tensor = []

    for gaussian_filter in gaussian_filters:
        
        if gaussian_filter == []:
            continue

        gaussian_filters_tensor.append(gaussian_filter)

    gaussian_filters_tensor = np.asarray(gaussian_filters_tensor)
    gaussian_filters_tensor = normalize_filters(gaussian_filters_tensor)

    # Center the Gaussian filters to the middle of visual field.
    _,n_pixels_y,n_pixels_x = gaussian_filters_tensor.shape
    gaussian_filters_tensor = center_filters(
                                    gaussian_filters_tensor,
                                    n_pixels_y,n_pixels_x
                              )
    encoding_filters = get_encoding_filters_ns(gaussian_filters_tensor)

    filter_dict = dict()
    filter_dict['encoding_filters'] = encoding_filters
    filter_dict['decoding_filters'] = np.ndarray
    filter_dict['cells_reconstruction'] = cells_reconstruction
    filter_dict['nl_popt'] = nl_popt

    decoder = learn_recon_filters_ns(encoding_filters)
    filter_dict['decoding_filters'] = decoder

    return filter_dict

def get_sig_stixels_encoder(encoder,excluded_cell_inds=None):
    """
    Gets sig stixels by taking union of all nonzero encoder pixels.

    TODO: doc
    """
    sig_stixels_all = set()

    for i in range(encoder.shape[0]):
        
        if excluded_cell_inds is not None and i in excluded_cell_inds:
            continue

        sig_stixels = np.nonzero(encoder[i,:].ravel())[0]

        for s in sig_stixels:
            sig_stixels_all.add(s)

    sig_stixels_all = np.asarray(sorted(list(sig_stixels_all)))

    return sig_stixels_all

def get_dictionary_variance(dictionary_matrix,decoder):
    """
    Computes variance of the dictionary

    TODO: document.
    """
    n_cells,n_pixels = decoder.shape
    decoder_norm_sq = linalg.norm(decoder,axis=1)**2
    var = np.matmul(dictionary_matrix * (1 - dictionary_matrix),
                    decoder_norm_sq
                )

    return var

def get_optimal_recon(decoder,stimulus,discretize=True):
    """
    Gets the optimal reconstruction by solving for the optimal responses under
    squared loss. Discretizes the spike vector if indicated.

    TODO: document
    """
    optimal_decoded_stimuli = []
    n_cells = decoder.shape[0]
    
    for i in range(stimulus.shape[1]):
        x = cp.Variable(n_cells)
        objective = cp.Minimize(
                        cp.sum_squares(
                            decoder.T @ x - stimulus[...,i].ravel()
                        )
                    )
        constraints = [x >= 0] # non-negative firing rate
        prob = cp.Problem(objective,constraints)
        prob.solve();

        if discretize:
            x.value = np.floor(x.value)

        optimal_decoded_stimuli.append(decoder.T@x.value)
    
    return np.stack(optimal_decoded_stimuli,axis=1)

def n_mse(stimuli,decoded_stimuli,sig_stixels):
    """
    Computes normalized squared error between the stimulus and decoded stimulus
    within the set of nonzero pixels.

    TODO: document
    """
    return (linalg.norm(
                stimuli[sig_stixels,:]
                - decoded_stimuli[sig_stixels,:],
                axis=0)**2 / linalg.norm(stimuli[sig_stixels,:],axis=0)**2
            )

def get_frac_incorrect_pix(stimuli,decoded_stimuli,sig_stixels):
    """
    Computes the fraction of incorrect pixels between the reconstruction 
    and the original stimulus within a specified set of pixels.

    TODO: document
    """
    return (np.argwhere(np.sign(stimuli[sig_stixels]) !=
                np.sign(decoded_stimuli[sig_stixels])).shape[0] / 
                sig_stixels.shape[0])

def get_sign_diff_map(stimuli,decoded_stimuli,sig_stixels):
    """
    Get the sign difference map.

    TODO: document
    """
    sign_map = np.zeros(stimuli.shape)
    diff_inds = np.argwhere(np.sign(stimuli[sig_stixels]) !=
                np.sign(decoded_stimuli[sig_stixels])).flatten()
    sign_map[sig_stixels[diff_inds]] = -1
    same_inds = np.argwhere(np.sign(stimuli[sig_stixels]) ==
                np.sign(decoded_stimuli[sig_stixels])).flatten()
    sign_map[sig_stixels[same_inds]] = 1

    return sign_map

def greedy_stim(dictionary_matrix,decoder,var_dict,
                sig_stixels_all,test_stimulus_vector): 
    """
    Function to run the greedy stimulation algorithm. Uses torch on GPU 
    for high performance relative to numpy/CPU.

    TODO: document
    """
    n_cells,n_pixels = decoder.shape
    decoded_stimuli = np.zeros((n_pixels,1))
    responses_all = np.zeros(n_cells)
    decoded_stimuli_partial = decoder.T@dictionary_matrix.T
    error_log = []
    element_log = []

    # Initialize the clock and refractory log
    t = 0
    cnt = 0
    refractory_log = np.zeros(n_cells)
    
    # Put on the GPU.
    decoder = torch.DoubleTensor(decoder).to(cfg.GPU)
    decoded_stimuli = torch.DoubleTensor(decoded_stimuli).to(cfg.GPU)
    decoded_stimuli_partial = torch.DoubleTensor(decoded_stimuli_partial
                                            ).to(cfg.GPU)
    test_stimulus_vector = torch.DoubleTensor(
                            test_stimulus_vector.reshape(n_pixels,
                                                        1)).to(cfg.GPU)
    var_dict = torch.DoubleTensor(var_dict).to(cfg.GPU)
    dictionary_matrix = torch.DoubleTensor(dictionary_matrix).to(cfg.GPU)
    responses_all = torch.DoubleTensor(responses_all).to(cfg.GPU)
    refractory_log = torch.DoubleTensor(refractory_log).to(cfg.GPU)

    while True:

        # Get perception, take the mean squared error and add on the variance.
        cumul_p = torch.add(decoded_stimuli,decoded_stimuli_partial)
        error = torch.add((
                torch.norm(
                torch.sub(cumul_p[sig_stixels_all],
                        test_stimulus_vector[sig_stixels_all]),dim=0)**2 / 
                    torch.norm(test_stimulus_vector[sig_stixels_all])**2),
                                var_dict
                        )
        elements = torch.argsort(error)

        # Get one with minimal error. 
        j = 0

        while torch.any(
            refractory_log[torch.where(
            dictionary_matrix[elements[j],:] 
            > cfg.MIN_SPIKE_PROB)[0]] > t):
            j +=1

            # If this brings us uphill, hold constant.
            if error[elements[j]] > error_log[-1]:
                j = torch.where(elements 
                                == elements.shape[0]-1)[0].cpu().numpy()[0]
                break

            # If we reach the end, break out. 
            if j == elements.shape[0]-1:
                break
                
        chosen_element = elements[j]
        error_log.append(error[chosen_element].item())
        element_log.append(chosen_element.item())

        # Add to the decoded stimuli
        decoded_stimuli += torch.reshape(
                           torch.matmul(
                           torch.t(decoder),
                           torch.t(dictionary_matrix[chosen_element,:])),
                           (n_pixels,1)
                           )
        responses_all += dictionary_matrix[chosen_element,:]
        
        if cnt % 50 == 0:
            print("stimulation: %s, error: %.3f"%(cnt,error[chosen_element]))

        # Update refractory log based on stimulation.
        refractory_log[torch.where(
                    dictionary_matrix[chosen_element,:] 
                    > cfg.MIN_SPIKE_PROB
                    )[0]] = cfg.REFRACTORY_PERIOD + t
        t += cfg.STIMULATION_DT
        cnt +=1

        # Exit the loop when the no-op element is chosen for a critical period. 
        if (len(element_log) > cfg.REFRACTORY_WINDOW and 
            np.all(np.asarray(element_log[-cfg.REFRACTORY_WINDOW:]) 
            == element_log[-1]) and
            torch.all(dictionary_matrix[element_log[-1],:] == 0.0)):
            print('converged after %s total stimulations.'%str(cnt))
            break
    
    # Put everything on CPU to avoid memory issues and write to dictionary.
    greedy_dict = dict()
    greedy_dict['decoded_stimuli'] = decoded_stimuli.cpu().numpy()
    greedy_dict['error_log'] = error_log
    greedy_dict['element_log'] = element_log
    greedy_dict['responses_all'] = responses_all.cpu().numpy()
    
    return greedy_dict