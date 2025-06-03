import numpy as np
import os
import sys
import whitenoise.random_noise as rn
import config as cfg
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import rawmovie as rm

"""
Miscellaneous io utilities.
"""

def get_raw_movie(raw_movie_path,n_frames):
    """
    Utility for getting natural scenes raw movies from disk 
    Parameters:
        raw_movie_path: path to stimulus
        n_frames: number of unique frames (int)
    Returns: 
        stimulus tensor of size n,y,x,c 
    """
    assert os.path.isfile(raw_movie_path), "Stimulus provided not found."

   # Initialize the stimulus object.
    rm_obj = rm.RawMovieReader(raw_movie_path)
    stimulus_tensor,_ = rm_obj.get_frame_sequence(0,n_frames)

    return stimulus_tensor

def get_movie_xml_str(vcd,seed,contrast,independent):
    """
    Constructs stimulus XML str of the form
     RGB/BW-stixel-interval-contrast-seed.xml
    Parameters:
        vcd: vision object
        seed: seed of stimulus
        contrast: contrast of stimulus
    Returns:
        stimulus movie str
    """
    stixel = int(vcd.runtimemovie_params.pixelsPerStixelX)
    interval = int(vcd.runtimemovie_params.interval)
    
    if independent:
        movie_xml_str = f'RGB-{stixel}-{interval}-{contrast}-{seed}.xml'
    else:
        movie_xml_str = f'BW-{stixel}-{interval}-{contrast}-{seed}.xml'
    
    return movie_xml_str

def get_celltypes_dict(class_path,lower=True):
    """
    Gets the celltype dictionary mapping IDs to celltypes from
    a text file.

    Parameters:
        class_path: full path to text file of cell types
        lower: boolean indicating whether to lowercase the strings

    Returns:
        dictionary mapping IDs to celltype.
    """

    f = open(class_path)
    celltypes_dict = dict()

    for j in f:
        tmp = ""

        for jj,substr in enumerate(j.split()[1:]):
            tmp +=substr

            if jj < len(j.split()[1:])-1:
                tmp += " "

        if lower:
            tmp = tmp.lower()

        celltypes_dict[int(j.split()[0])] = tmp

    f.close()

    return celltypes_dict

def get_stimulus(movie_xml_str,n_frames,resample=True,
                normalize=False,center=False,grayscale=False):
    """
    Gets a white noise visual stimulus from an xml str

    Parameters:
        movie_xml_str: white noise movie xml
        n_frames: number of frames
        resample: boolean indicating whether to upsample based on interval.
        normalize: boolean to normalize in 0,1
        center: boolean to mean subtract.
        grayscale: boolean indicating whether to grayscale.
    Returns:
        tensor of size frames x height x width x color channels
    """

    # Check that the stimulus exists, and initialize object.
    if cfg.XML_EXT not in movie_xml_str:
        movie_xml_str += cfg.XML_EXT

    movie_xml_path = os.path.join(cfg.MOVIE_PATH,movie_xml_str)

    if not os.path.isfile(movie_xml_path):
        raise OSError("%s not found."%movie_xml_str)

    # Load the stimulus and resample if necessary.
    rn_obj = rn.RandomNoiseFrameGenerator.construct_from_xml(movie_xml_path)

    # TODO: FIXME: use ericwu updated code; kludge to deal with overflow
    #stimulus = rn_obj.generate_block_of_frames(n_frames)
    stimulus = []

    for i in range(n_frames):
        stimulus.append(rn_obj.generate_next_frame())
    
    stimulus = np.asarray(stimulus)

    if normalize:
        stimulus /= cfg.MONITOR_BIT_DEPTH

    if center:
        stimulus -= np.mean(stimulus.ravel())

    if grayscale:
        stimulus = np.mean(stimulus,axis=-1,keepdims=True)

    interval = int(movie_xml_str.split('-')[2]) # Assumes RGB-pix-int structure

    if interval == 1 or not resample:
        return stimulus

    return np.repeat(stimulus,interval,axis=0)