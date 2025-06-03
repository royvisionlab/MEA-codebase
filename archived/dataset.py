import enum
import sys
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import visionloader as vl
import os
import config as cfg
import io_util as io
import importlib as il
import lnp 
import fitting as fit
from scipy.io import loadmat
import multiprocessing as mp
import pdb
from progressbar import *
from scipy.stats import gmean
import spiking as spkg
import recon
import cv2
il.reload(cfg);
il.reload(fit);

"""
Dataset object: contains critical analysis info for each data set.
Wrapper around visionloader with other bells and whistles.
"""

class Dataset(object):
    
    def __init__(self,dataset,wn_datarun=None,estim_datarun=None,
                 movie_xml_str=None):
        """
        Constructor method
        Parameters:
            dataset: dataset name
            wn_datarun: white noise datarun
            estim_datarun: estim data run
            movie_xml_str: white noise movie XML
        """
        self.dataset = dataset
        self.wn_datarun = wn_datarun
        self.estim_datarun = estim_datarun
        self.movie_xml_str = movie_xml_str
        self.wn_datapath = None
        self.estim_datapath = None
        self.vcd = None
        self.celltypes_dict = None
        self.stas = None
        self.binned_spikes = None
        self.cellids = None
        self.n_frames = None
        self.stimulus = None
        self.generator_signals = None
        self.nonlinearity = None
        self.spatial_maps = None
        self.timecourses = None
        self.sta_spatial_fit = None
        self.sig_stixels = None
        self.cellids = None
        self.timecourse_fit = None
        self.nonlinearity_fit = None
        self.time_zero = None
        self.gsort_results = None
        self.dictionary = None
        self.sta_time_ms = None
        self.rf_sizes = None
        self.acfs = None
        self.acvs = None
        self.bundles = None
        self.selectivity_dict = None
        self.recon_filter_dict = None
        self.recon_filter_dict_ns = None
        self.pruned_dictionary = None
        self.greedy_results = None

        if wn_datarun is not None:
            self.wn_datapath= os.path.join(cfg.ANALYSIS_PARENT,
                                            dataset,wn_datarun)
        
        if estim_datarun is not None:
            self.estim_datapath = os.path.join(cfg.ANALYSIS_PARENT,
                                            dataset,estim_datarun)
    
    def set_vcd(self,include_params=True,
                     include_ei=True,
                     include_neurons=True,
                     include_runtimemovie_params=True):
        """
        Sets the vision data table object
        """
        if self.wn_datapath is not None:
            self.vcd = vl.load_vision_data(self.wn_datapath,
                        os.path.basename(self.wn_datapath),
                        include_params=include_params,
                        include_ei=include_ei,
                        include_neurons=include_neurons,
                        include_runtimemovie_params=include_runtimemovie_params)
        
        return None

    def get_vcd(self):
        """
        Returns the vcd
        """
        if self.vcd is None:
            self.set_vcd()
        
        return self.vcd
        
    def set_celltypes_dict(self):
        """
        Sets the cell type dictionary from the text file on disk.
        """
        class_path = os.path.join(self.wn_datapath,
                        cfg.CLASS_FNAME%os.path.basename(self.wn_datarun))
        self.celltypes_dict = io.get_celltypes_dict(class_path)

        return None

    def get_celltypes_dict(self):
        """
        Gets the dict.
        """
        if self.celltypes_dict is None:
            self.set_celltypes_dict()
        
        return self.celltypes_dict

    def get_cellids(self):
        """
        Wrapper around visionloader
        """
        if self.cellids is not None:
            return self.cellids

        if self.vcd is None:
            self.set_vcd()
        
        self.cellids = sorted(self.vcd.get_cell_ids()) # sort by convention!

        return self.cellids

    def set_cellids(self):
        """
        Sets to the object.
        """
        if self.cellids is None:
            _ = self.get_cellids()
        
        return None

    def get_binned_spikes(self):
        """
        Wrapper around lnp.py. Gets binned spikes and also makes it play 
        nice with interval.
        """

        if self.binned_spikes is not None:
            return self.binned_spikes

        if self.vcd is None:
            self.set_vcd()

        if self.cellids is None: 
            self.cellids = sorted(self.vcd.get_cell_ids())

        # Get the binned spikes and truncate to be an even number of frames.
        if self.binned_spikes is None:
            self.binned_spikes = lnp.get_binned_spikes(self.vcd,self.cellids)
        
        if self.n_frames is None:
            self.n_frames = self.binned_spikes.shape[0]

        # Make sure number of frames is multiple of interval. 
        interval = int(self.vcd.runtimemovie_params.interval) 

        while self.n_frames % interval != 0:
            self.n_frames -= 1
        
        # Once it's even, truncate the binned spikes.
        self.binned_spikes[0:self.n_frames,:]

        return None

    def set_binned_spikes(self):
        """
        Sets in the object.
        """
        if self.binned_spikes is None:
            _ = self.get_binned_spikes()

        return None

    def compute_stas(self,write=True):
        """
        Wrapper around STA calculation in lnp.py, sets the tensor to the object.
        """

        if self.vcd is None:
            self.set_vcd()

        if self.cellids is None: 
            self.cellids = sorted(self.vcd.get_cell_ids())

        # Get the binned spikes and truncate to be an even number of frames.
        if self.binned_spikes is None:
            self.set_binned_spikes()
        
        interval = int(self.vcd.runtimemovie_params.interval) 

        if self.stimulus is None:
            self.stimulus = io.get_stimulus(self.movie_xml_str,
                                            int(self.n_frames / interval))

        self.stas = lnp.compute_sta_torch(self.stimulus,
                                          self.binned_spikes,
                                          interval)

        if write:
            dirout = os.path.join(cfg.LNP_PARENT,
                                    self.dataset,self.wn_datarun)

            if not os.path.isdir(dirout):
                os.makedirs(dirout)
            
            fnameout = os.path.join(dirout,
                                cfg.STA_FNAME%os.path.basename(self.wn_datarun))
            np.save(fnameout,self.stas)

        return None
    
    def compute_timecourses(self,write=True):
        """
        Wrapper around lnp.py. Normalizes timecourses by default.
        """
        if self.stas is None:
            self.set_stas()

        self.timecourses = lnp.compute_timecourses(self.stas,norm=True)

        if write:
            dirout = os.path.join(cfg.LNP_PARENT,
                                    self.dataset,self.wn_datarun)

            if not os.path.isdir(dirout):
                os.makedirs(dirout)
            
            fnameout = os.path.join(dirout,
                                cfg.TC_FNAME%os.path.basename(self.wn_datarun))
            np.save(fnameout,self.timecourses)
        
        return None

    def compute_generator_signals(self,write=True):
        """
        lnp.py wrapper
        """

        if self.vcd is None:
            self.set_vcd()

        if self.cellids is None: 
            self.cellids = sorted(self.vcd.get_cell_ids())

        # Get the binned spikes and truncate to be an even number of frames.
        if self.binned_spikes is None:
            self.get_binned_spikes()
        
        interval = int(self.vcd.runtimemovie_params.interval) 

        if self.stimulus is None:
            self.stimulus = io.get_stimulus(self.movie_xml_str,
                                            int(self.n_frames / interval))
        
        self.generator_signals = lnp.compute_generator_signal_torch(
                                     self.stimulus,self.stas,interval
                                 )   

        if write:
            dirout = os.path.join(cfg.LNP_PARENT,
                                    self.dataset,self.wn_datarun)

            if not os.path.isdir(dirout):
                os.makedirs(dirout)
            
            fnameout = os.path.join(dirout,
                                cfg.GS_FNAME%os.path.basename(self.wn_datarun))
            np.save(fnameout,self.generator_signals)

            fnameout = os.path.join(dirout,
                                cfg.BS_FNAME%os.path.basename(self.wn_datarun))
            np.save(fnameout,self.binned_spikes)

        return None
    
    def compute_nonlinearity(self,stdz_g=True,write=True):
        """
        Wrapper around lnp.py. 
        """

        if self.vcd is None:
            self.set_vcd()

        if self.binned_spikes is None:
            self.set_binned_spikes()
        
        if self.generator_signals is None:
            self.set_generator_signals()

        monitor_fs = int(self.vcd.runtimemovie_params.monitorFrequency)
        mean_g,mean_fr = lnp.compute_nonlinearity(self.generator_signals,
                                                     self.binned_spikes.T,
                                                     monitor_fs,stdz_g=stdz_g)
        
        nl_dict = dict()
        nl_dict['mean_generator_signals'] = mean_g
        nl_dict['mean_spike_counts'] = mean_fr
        self.nonlinearity = nl_dict

        if write:
            dirout = os.path.join(cfg.LNP_PARENT,
                        self.dataset,self.wn_datarun)

            if not os.path.isdir(dirout):
                os.makedirs(dirout)
            
            fnameout = os.path.join(dirout,
                                cfg.NL_FNAME%os.path.basename(self.wn_datarun))
            np.save(fnameout,self.nonlinearity)
    
    def get_nonlinearity(self):
        """
        Gets from disk.
        """
        if self.nonlinearity is not None:
            return self.nonlinearity
        
        nl_path = os.path.join(cfg.LNP_PARENT,self.dataset,
                            self.wn_datarun,
                            cfg.NL_FNAME%os.path.basename(self.wn_datarun))

        try:
            self.nonlinearity = np.load(nl_path,allow_pickle=True).item()
        except:
            print(f"couldn't load nonlinearity for {self.dataset},\
                                                     {self.wn_datarun}")

            return None
        
        return self.nonlinearity

    def set_nonlinearity(self):
        """
        Sets to the object.
        """
        if self.nonlinearity is None:
            _ = self.get_nonlinearity()
        
        return None

    def get_generator_signals(self):
        """
        Gets generator signals from disk.
        """
        if self.generator_signals is not None:
            return self.generator_signals
        
        gs_path = os.path.join(cfg.LNP_PARENT,self.dataset,
                            self.wn_datarun,
                            cfg.GS_FNAME%os.path.basename(self.wn_datarun))

        try:
            self.generator_signals = np.load(gs_path)
        except:
            print(f"couldn't load gen. signals for {self.dataset}, {self.wn_datarun}")

            return None
        
        return self.generator_signals

    def set_generator_signals(self):
        """
        Sets in the object
        """
        if self.generator_signals is None:
            _ = self.get_generator_signals()
        
        return None

    def get_stas(self):
        """
        Loads existing STAs from disk. Does NOT compute them.
        """
        if self.stas is not None:
            return self.stas

        stas_path = os.path.join(cfg.LNP_PARENT,self.dataset,
                                 self.wn_datarun,
                                 cfg.STA_FNAME%os.path.basename(self.wn_datarun))
        
        try:
            self.stas = np.load(stas_path)
        except:
            print(f"couldn't load STAs for {self.dataset}, {self.wn_datarun}")

            return None

        return self.stas

    def set_stas(self):
        """
        Sets stas in the object, returns None.
        """
        _ = self.get_stas()

        return None

    def get_timecourses(self):
        """
        Loads timecourses.
        """
        if self.timecourses is not None:
            return self.timecourses
        
        tcs_path = os.path.join(cfg.LNP_PARENT,self.dataset,
                                 self.wn_datarun,
                                 cfg.TC_FNAME%os.path.basename(self.wn_datarun))
        
        try:
            self.timecourses = np.load(tcs_path)
        except:
            print(f"couldn't load timecourses for {self.dataset}, {self.wn_datarun}")

            return None

        return self.timecourses
    
    def set_timecourses(self):
        """
        Sets timecoursees in the object, returns None.
        """
        if self.timecourses is None:
            _ = self.get_timecourses()

        return None

    def compute_spatial_maps(self,write=True):
        """
        Wrapper around lnp.py. Passes a normalized STA to get a spatial map.
        """
        if self.celltypes_dict is None:
            self.set_celltypes_dict()
        
        if self.stas is None:
            self.set_stas()
        
        self.spatial_maps = lnp.compute_spatial_maps(lnp.normalize_sta(self.stas),
                                                     self.celltypes_dict)

        if write:
            dirout = os.path.join(cfg.LNP_PARENT,
                        self.dataset,self.wn_datarun)

            if not os.path.isdir(dirout):
                os.makedirs(dirout)
            
            fnameout = os.path.join(dirout,
                                cfg.SM_FNAME%os.path.basename(self.wn_datarun))
            np.save(fnameout,self.spatial_maps)
        
        return None

    def get_spatial_maps(self):
        """
        Gets maps from disk
        """
        if self.spatial_maps is not None:
            return self.spatial_maps

        sms_path = os.path.join(cfg.LNP_PARENT,self.dataset,
                                 self.wn_datarun,
                                 cfg.SM_FNAME%os.path.basename(self.wn_datarun))
        
        try:
            self.spatial_maps = np.load(sms_path)
        except:
            print(f"couldn't load spatial maps for {self.dataset}, {self.wn_datarun}")

            return None
        
        return self.spatial_maps
    
    def set_spatial_maps(self):
        """
        Sets to the object.
        """
        if self.spatial_maps is None:
            _ = self.get_spatial_maps()
        
        return None

    def get_sig_stixels(self):
        """
        Wrapper around lnp.

        Returns a list of sig stixels, where each element is for a cell
        """
        if self.sig_stixels is not None:
            return self.sig_stixels
        
        if self.stas is None:
            self.set_stas()
        
        sig_stixels_list = []

        for i in range(self.stas.shape[0]):
            sig_stixels = lnp.calculate_sig_stixels(
                          self.stas[i,...] - np.mean(self.stas[i,...])
                          )
            sig_stixels_list.append(np.argwhere(sig_stixels[0]))
        
        self.sig_stixels = sig_stixels_list

        return self.sig_stixels

    def set_sig_stixels(self):
        """
        Sets to the object
        """
        if self.sig_stixels is None:
            _ = self.get_sig_stixels()
        
        return None

    def fit_sta_spatial(self,write=True):
        """
        Fits spatial component of STA with a 2D Gaussian.
        For each cell, gets the spatial component of STA and fits 2D Gaussian.
        """

        if self.vcd is None:
            self.set_vcd()
        
        if self.celltypes_dict is None:
            self.set_celltypes_dict()

        if self.spatial_maps is None:
            self.set_spatial_maps()

        if self.sig_stixels is None:
            self.set_sig_stixels()

        # Loop through each of the spatial maps and get Gaussian fits.
        popt_list = []

        for i in range(self.spatial_maps.shape[0]):
            sig_stixels = self.sig_stixels[i]

            if sig_stixels.shape[0] == 0:
                popt_list.append(None)
                continue

            # Fit the 2D Gaussian and cache the parameters.
            popt = fit.fit_gaussian2d(
                     self.spatial_maps[i,...],
                     sig_stixels
                   )
            popt_list.append(popt)
        
        self.sta_spatial_fit = popt_list

        if write:
            dirout = os.path.join(cfg.LNP_PARENT,
                                    self.dataset,self.wn_datarun)

            if not os.path.isdir(dirout):
                os.makedirs(dirout)
            
            fnameout = os.path.join(
                               dirout,
                               cfg.STA_FIT_FNAME%os.path.basename(self.wn_datarun)
                               )
            
            tmp = dict()
            tmp['sta_spatial_fit'] = self.sta_spatial_fit
            np.save(fnameout,tmp,allow_pickle=True)
    
        return None

    def get_sta_spatial_fit(self):
        """
        Gets from disk.
        """
        if self.sta_spatial_fit is not None:
            return self.sta_spatial_fit

        spatial_fit_path = os.path.join(
                            cfg.LNP_PARENT,self.dataset,
                            self.wn_datarun,
                            cfg.STA_FIT_FNAME%os.path.basename(self.wn_datarun)
                            )
        
        try:
            self.sta_spatial_fit = np.load(
                                        spatial_fit_path,
                                        allow_pickle=True
                                   )
            self.sta_spatial_fit = self.sta_spatial_fit.item()['sta_spatial_fit']

        except:
            print(f"couldn't load STA fits for {self.dataset},\
                                       {self.wn_datarun}")
            return None
        
        return self.sta_spatial_fit

    def set_sta_spatial_fit(self):
        """
        Sets to the object
        """
        if self.sta_spatial_fit is None:
            _ = self.get_sta_spatial_fit()
        
        return None

    def compute_rf_size(self):
        """
        Takes the geometric mean of the two std parameters in the Gaussian 
        fit. Multiplies by the constant in the config, the stixel size and 
        the microns per pixel. Just one reasonable way of doing this.
        """
        if self.rf_sizes is not None:
            return None

        if self.sta_spatial_fit is None:
            self.set_sta_spatial_fit()
        
        if self.vcd is None:
            self.set_vcd()
        
        stixel_size = self.vcd.runtimemovie_params.pixelsPerStixelX

        # Loop through the Gaussian fits and compute gmean
        rf_sizes = []

        for i in range(len(self.sta_spatial_fit)):
            sta_spatial_fit = self.sta_spatial_fit[i]

            if sta_spatial_fit is None:
                rf_sizes.append(None)
                continue
            
            # Get the STD params.
            sigma_x = sta_spatial_fit[3]
            sigma_y = sta_spatial_fit[4]
            rf_size = gmean(
                        np.array([sigma_x,sigma_y])
                      ) * stixel_size * cfg.MICRONS_PER_PIXEL * cfg.N_SIGMAS_RF
            rf_sizes.append(rf_size)
        
        self.rf_sizes = rf_sizes

        return None
        
    def fit_sta_timecourses(self,write=True):
        """
        Fits time courses: wrapper around fitting.py. For every non-blue 
        cell type, fits to green; otherwise, fits to blue.
        """

        if self.timecourses is None:
            self.set_timecourses()

        if self.celltypes_dict is None:
            self.set_celltypes_dict()

        if self.cellids is None:
            self.set_cellids()
        
        # Loop through each time course and fit.
        popt_list = []
        fit_list = []

        for i in range(self.timecourses.shape[0]):

            # Choose channel based on celltype.
            celltype = self.celltypes_dict[self.cellids[i]]

            if "blue" in celltype or "sbc" in celltype:
                channel = 2 # b
            else:
                channel = 1 # g

            popt,tc_fit = fit.fit_sta_timecourse(self.timecourses[i,:,channel])
            popt_list.append(popt)
            fit_list.append(tc_fit)
        
        tmp = dict()
        tmp['popt'] = popt_list
        tmp['fits'] = fit_list

        self.timecourse_fit = tmp

        if write:
            dirout = os.path.join(cfg.LNP_PARENT,
                                    self.dataset,self.wn_datarun)

            if not os.path.isdir(dirout):
                os.makedirs(dirout)
            
            fnameout = os.path.join(
                               dirout,
                               cfg.TC_FIT_FNAME%os.path.basename(self.wn_datarun)
                               )
            np.save(fnameout,self.timecourse_fit,allow_pickle=True)
    
        return None

    def get_sta_timecourses_fits(self):
        """
        Gets from disk.
        """
        if self.timecourse_fit is not None:
            return self.timecourse_fit
        
        timecourse_fit_path = os.path.join(
                    cfg.LNP_PARENT,self.dataset,
                    self.wn_datarun,
                    cfg.TC_FIT_FNAME%os.path.basename(self.wn_datarun)
                    )
        
        try:
            self.timecourse_fit = np.load(
                                        timecourse_fit_path,
                                        allow_pickle=True
                                   )
            self.timecourse_fit = self.timecourse_fit.item()

        except:
            print(f"couldn't load timecourse fits for {self.dataset},\
                                       {self.wn_datarun}")
            return None
        
        return self.timecourse_fit

    def set_sta_timecourses_fits(self):
        """
        Sets to the object
        """
        if self.timecourse_fit is None:
            _ = self.get_sta_timecourses_fits()
        
        return None

    def fit_nonlinearity(self,write=True):
        """
        Fitting wrapper
        """
        if self.nonlinearity is None:
            self.set_nonlinearity()
        
        # Loop through each cell and fit the nonlinearity.
        gen_fit_list = []
        fr_fit_list = []
        popt_list = []
        mean_g = self.nonlinearity['mean_generator_signals']
        mean_fr = self.nonlinearity['mean_spike_counts']

        for i in range(mean_g.shape[0]):

            if np.all(np.isnan(mean_g[i])):
                gen_fit_list.append(None)
                fr_fit_list.append(None)
                popt_list.append(None)
                continue
            
            # Get the largest generator signal value.
            max_g = np.max(np.abs(mean_g[i])) * cfg.MAX_G_MULT
            gen_fit,fr_fit,popt = fit.fit_nonlinearity(
                                        mean_g[i],mean_fr[i],max_g
                                   )
            gen_fit_list.append(gen_fit)
            fr_fit_list.append(fr_fit)
            popt_list.append(popt)
        
        tmp = dict()
        tmp['g'] = gen_fit_list
        tmp['fr'] =  fr_fit_list 
        tmp['popt'] = popt_list

        self.nonlinearity_fit = tmp
        
        if write:
            dirout = os.path.join(cfg.LNP_PARENT,
                        self.dataset,self.wn_datarun)

            if not os.path.isdir(dirout):
                os.makedirs(dirout)
            
            fnameout = os.path.join(
                               dirout,
                               cfg.NL_FIT_FNAME%os.path.basename(
                                                    self.wn_datarun
                                                )
                               )
            np.save(fnameout,self.nonlinearity_fit,allow_pickle=True)
        
        return None

    def get_nonlinearity_fit(self):
        """
        Gets from disk
        """
        if self.nonlinearity_fit is not None:
            return self.nonlinearity_fit

        nonlinearity_fit_path = os.path.join(
                        cfg.LNP_PARENT,self.dataset,
                        self.wn_datarun,
                        cfg.NL_FIT_FNAME%os.path.basename(self.wn_datarun)
                        )
        
        try:
            self.nonlinearity_fit = np.load(
                                        nonlinearity_fit_path,
                                        allow_pickle=True
                                   )
            self.nonlinearity_fit = self.nonlinearity_fit.item()

        except:
            print(f"couldn't load nonlinearity fits for {self.dataset},\
                                       {self.wn_datarun}")
            return None
        
        return self.nonlinearity_fit
    
    def set_nonlinearity_fit(self):
        """
        Sets to object
        """
        if self.nonlinearity_fit is None:
            _ = self.get_nonlinearity_fit()
        
        return None
    
    def get_sta_time_ms(self):
        """
        Gets the time of the time course based on interval, depth and 
        refresh period of monitor. Returns this in the native format of the 
        STA (i.e. 0 time is first).
        """
        if self.sta_time_ms is not None:
            return self.sta_time_ms
        
        if self.timecourses is None:
            self.set_timecourses()
        
        if self.vcd is None:
            self.set_vcd()
        
        monitor_fs = self.vcd.runtimemovie_params.monitorFrequency
        interval = self.vcd.runtimemovie_params.interval
        depth = self.timecourses.shape[1]
        self.sta_time_ms = (np.arange(0,depth) * (1 / monitor_fs) * 
                            interval * 1000)
        
        return self.sta_time_ms
    
    def set_sta_time_ms(self):
        """
        Sets to the object.
        """
        if self.sta_time_ms is None:
            _ = self.get_sta_time_ms()
        
        return None

    def compute_time_zero(self,write=True):
        """
        Computes time of zero crossing for the time course fits. Updates the file
        by writing to a new field in the dict.
        """
        if self.timecourse_fit is None:
            self.set_sta_timecourses_fits()
        
        if self.celltypes_dict is None:
            self.set_celltypes_dict()

        if self.cellids is None:
            self.set_cellids()

        if self.vcd is None:
            self.set_vcd()
        
        # Get some monitor params needed.
        monitor_fs = self.vcd.runtimemovie_params.monitorFrequency
        interval = self.vcd.runtimemovie_params.interval

        # Iterate through cells and compute time z 
        time_zero_list = []

        for cc,cell in enumerate(self.cellids):
            popt = self.timecourse_fit['popt'][cc]

            if popt is None:
                time_zero_list.append(None)
                continue
            
            celltype = self.celltypes_dict[cell]
            t_fit = self.timecourse_fit['fits'][cc]
        
            # Resolve extrema based on cell types 
            if "on" in celltype:
                lower = np.argmax(t_fit)
                upper = np.argmin(t_fit)
            else:
                lower = np.argmin(t_fit)
                upper = np.argmax(t_fit)

            f = fit.lfp_anon(*popt)

            try:
                z,_ = fit.bissection_search_root(f,lower+1,upper+1)
                z = t_fit.shape[0] - z

                # Convert to ms
                z_ms = (((1 / (monitor_fs / interval)) *
                        1000 * (t_fit.shape[0] - 1 - z)))
                time_zero_list.append(z_ms)

            except:
                time_zero_list.append(None)
                continue
        
        tmp = dict()
        tmp['time_zero'] = time_zero_list
        self.time_zero = tmp.copy()

        del tmp

        # Add this to the timecourse fit part of object.
        if write:
            tmp = self.timecourse_fit
            tmp['z_ms'] = time_zero_list
            self.timecourse_fit = tmp

            dirout = os.path.join(cfg.LNP_PARENT,
            self.dataset,self.wn_datarun)

            if not os.path.isdir(dirout):
                os.makedirs(dirout)
            
            fnameout = os.path.join(
                               dirout,
                               cfg.TC_FIT_FNAME%os.path.basename(
                                                    self.wn_datarun
                                                )
                               )
            np.save(fnameout,self.timecourse_fit,allow_pickle=True)
        
        return None

    def compute_acfs(self,write=True):
        """
        Wrapper around spiking.py to compute acfs for all cells.
        """
        if self.binned_spikes is None:
            self.set_binned_spikes()
        
        if self.cellids is None:
            self.set_cellids()
        
        if self.vcd is None:
            self.set_vcd()
        
        # Initialize acf dict and loop through each cell, computing the ACF.
        acf_dict = dict()

        with mp.Pool(cfg.MP_N_THREADS) as pool:
            """
            acfs = pool.map(spkg.compute_acf,
                               [self.binned_spikes[:,i]
                                for i in range(self.binned_spikes.shape[1])]
                    )
            """
            acfs = pool.starmap(spkg.compute_acf,
                                [(self.vcd,cell) for cell in self.cellids]
                    )

        for cc,cell in enumerate(self.cellids):
            acf_dict[cell] = acfs[cc]
        
        self.acfs = acf_dict
        
        if write:
            dirout = os.path.join(cfg.ACF_PARENT,
            self.dataset,self.wn_datarun)

            if not os.path.isdir(dirout):
                os.makedirs(dirout)

            fnameout = os.path.join(
                    dirout,
                    cfg.ACF_FNAME%os.path.basename(self.wn_datarun)
                    )
            np.save(fnameout,self.acfs,allow_pickle=True)
        
        return None

    def get_acfs(self):
        """
        Gets from disk.
        """
        if self.acfs is not None:
            return self.acfs

        acf_path = os.path.join(
                cfg.ACF_PARENT,self.dataset,
                self.wn_datarun,
                cfg.ACF_FNAME%os.path.basename(self.wn_datarun)
                )
        
        try:
            self.acfs= np.load(acf_path,allow_pickle=True)
            self.acfs = self.acfs.item()

        except:
            print(f"couldn't load acfs {self.dataset},\
                                       {self.wn_datarun}")
            return None
        
        return self.acfs

    def set_acfs(self):
        """
        Sets to the object
        """
        if self.acfs is None:
            _ = self.get_acfs()
        
        return None

    def compute_acvs(self,write=True):
        """
        Wrapper around spiking.py. Computes axon conduction velocity.
        """
        if self.cellids is None:
            self.set_cellids()

        if self.vcd is None:
            self.set_vcd()
        
        # Initialize a list for the velocities.
        acv_list = []

        for cc,cell in enumerate(self.cellids):
            cell_ei = self.vcd.get_ei_for_cell(cell).ei
            acv = spkg.get_axonal_conduction_velocity(
                                        cell_ei,
                                        threshold_factor=cfg.ACV_THRESHOLD
                    )
            acv_list.append(acv)
        
        # Set it to the object and write it out.
        self.acvs = dict()
        self.acvs['acvs'] = acv_list
        
        if write:
            dirout = os.path.join(cfg.ACV_PARENT,
            self.dataset,self.wn_datarun)

            if not os.path.isdir(dirout):
                os.makedirs(dirout)

            fnameout = os.path.join(
                    dirout,
                    cfg.ACV_FNAME%os.path.basename(self.wn_datarun)
                    )
            np.save(fnameout,self.acvs,allow_pickle=True)
        
        return None

    def get_acvs(self):
        """
        Gets from disk.
        """
        if self.acvs is not None:
            return self.acvs

        acv_path = os.path.join(
                cfg.ACV_PARENT,self.dataset,
                self.wn_datarun,
                cfg.ACV_FNAME%os.path.basename(self.wn_datarun)
                )

        try:
            self.acvs = np.load(acv_path,allow_pickle=True)
            self.acvs = self.acvs.item()
        except:
            print(f"couldn't load acvs {self.dataset},\
                                       {self.wn_datarun}")
            return None
        
        return self.acvs

    def set_acvs(self):
        """
        Sets to the object.
        """
        if self.acvs is None:
            _ = self.get_acvs()
        
        return None

    def get_gsort_results(self):
        """
        Just loads from disk.
        """
        if self.gsort_results is not None:
            return self.gsort_results

        fnamein = os.path.join(
                cfg.PARENT_GSORT,
                self.dataset,
                self.estim_datarun,
                self.wn_datarun,
                cfg.GSORT_FNAME
              )
        try: 
            self.gsort_results = loadmat(fnamein)
        except:
            print("could not load gsort for {self.dataset},\
                                         {self.estim_datarun}")
            return None
        
        return self.gsort_results

    def set_gsort_results(self):
        """
        Sets to the object.
        """
        if self.gsort_results is None:
            _ = self.get_gsort_results()
        
        return None

    def make_dictionary(self,write=True):
        """
        Assembles dictionary for the ON and OFF parasol cells from the gsort
        output. Basically just indexes accordingly. In a dict obj, writes a 
        tensor for analysis and a matrix for greedy algorithm.
        """
        if self.gsort_results is None:
            self.set_gsort_results()

        if self.cellids is None:
            self.set_cellids()
        
        if self.celltypes_dict is None:
            self.set_celltypes_dict()
        
        # Get cell inds of interest and index.
        cell_inds = []

        for cc,cell in enumerate(self.cellids):

            if "parasol" not in self.celltypes_dict[cell]:
                continue
                
            cell_inds.append(cc)
        
        cell_inds = np.asarray(cell_inds)

        # Get gsort cosine probabilities and index.
        filtered_probs = self.gsort_results['filtered_probs']
        dictionary_tensor = filtered_probs[cell_inds,...]
        n_cells,n_elecs,n_amps = dictionary_tensor.shape
        dictionary_matrix = np.reshape(
                                dictionary_tensor,
                                (n_cells,n_elecs * n_amps)
                            )
        dictionary = dict()
        dictionary['dictionary_matrix_raw'] = dictionary_matrix 
        dictionary['dictionary_tensor_raw'] = dictionary_tensor
        self.dictionary = dictionary

        if write:
            fnameout = os.path.join(
                    cfg.PARENT_GSORT,
                    self.dataset,
                    self.estim_datarun,
                    self.wn_datarun,
                    cfg.DICTIONARY_FNAME
                )
            np.save(fnameout,dictionary,allow_pickle=True)

        return None

    def get_dictionary(self):
        """
        Gets dictionary from disk.
        """

        if self.dictionary is not None:
            return self.dictionary

        fnamein = os.path.join(
                cfg.PARENT_GSORT,
                self.dataset,
                self.estim_datarun,
                self.wn_datarun,
                cfg.DICTIONARY_FNAME
              )

        try: 
            self.dictionary = np.load(fnamein,allow_pickle=True).item()
        except:
            print("could not load dictionary for {self.dataset},\
                                         {self.estim_datarun}")
            return None
        
        return self.dictionary
    
    def set_dictionary(self):
        """
        Sets to the object
        """
        if self.dictionary is None:
            _ = self.get_dictionary()
        
        return None

    def fit_sigmoids(self,write=True):
        """
        Fits sigmoids to the raw dictionary tensor. Writes a dictionary with 
        the fits, and the extracted thresholds. OVERWRITES existing dictionary.
        """
        if self.dictionary is None:
            self.set_dictionary()
        
        # Loop through each cell and get the fits and extract thresholds.
        dictionary_tensor_fits = []
        thresholds_all = []
        slopes_all = []
        dictionary_tensor = self.dictionary['dictionary_tensor_raw']
        x = cfg.AMPLITUDES

        # Initialize progress bar.
        widgets = ['fitting sigmoids: ', Percentage(), ' ', Bar(marker='*',
                left='[',right=']'),' ', ETA()]
        pbar = ProgressBar(widgets=widgets, maxval=dictionary_tensor.shape[0])
        pbar.start()

        for i in range(dictionary_tensor.shape[0]):
            fits = []
            thresholds = []
            slopes = []

            with mp.Pool(cfg.MP_N_THREADS) as pool:
                out = pool.starmap(
                        fit.fit_sigmoid_mle,
                        [(x, fit.clean_probs(y,x))
                         for y in dictionary_tensor[i,...]]
                      )
            
            for popt,sigmoid_fit in out:

                # Set invalid thresholds to nan.
                slope,threshold = popt

                if threshold < 0 or threshold > cfg.AMPLITUDES[-1]:
                    threshold = np.nan
                    slope = np.nan

                # Set anything below epsilon to 0
                #eps_inds = np.argwhere(sigmoid_fit < cfg.EPSILON).flatten()
                #sigmoid_fit[eps_inds] = 0
                fits.append(sigmoid_fit)
                thresholds.append(threshold)
                slopes.append(slope)
            
            fits = np.asarray(fits)
            thresholds = np.asarray(thresholds)
            slopes = np.asarray(slopes)
            dictionary_tensor_fits.append(fits)
            thresholds_all.append(thresholds)
            pbar.update(i)

        dictionary_tensor_fits = np.asarray(dictionary_tensor_fits)
        thresholds_all = np.asarray(thresholds_all)
        slopes_all = np.asarray(slopes_all)
        pbar.finish() 

        # Update the object and write out.
        tmp = self.dictionary
        tmp['dictionary_tensor'] = dictionary_tensor_fits
        n_cells,n_elecs,n_amps = dictionary_tensor.shape
        dictionary_matrix_fits = np.reshape(
                                    dictionary_tensor_fits,
                                    (n_cells,n_elecs * n_amps)
                                 )
        tmp['dictionary_matrix'] = dictionary_matrix_fits
        tmp['thresholds'] = thresholds_all
        tmp['slopes'] = slopes_all
        self.dictionary = tmp

        if write:
            fnameout = os.path.join(
                        cfg.PARENT_GSORT,
                        self.dataset,
                        self.estim_datarun,
                        self.wn_datarun,
                        cfg.DICTIONARY_FNAME
                    )
            np.save(fnameout,self.dictionary,allow_pickle=True)
        
        return None

    def get_bundles(self):
        """
        Loads bundle thresholds from disk.
        """
        if self.bundles is not None:
            return self.bundles
        
        fnamein = os.path.join(
                cfg.PARENT_BUNDLE,
                cfg.BUNDLE_FNAME%self.dataset
              )

        try: 
            tmp = np.load(fnamein).astype(int)
        except:
            print(f"could not load bundles for {self.dataset},\
                                         {self.estim_datarun}")
            return None

        # Convert this to a dictionary mapping electrode to amp index.
        bundles = dict()

        for ee,elec_ind in enumerate(tmp[:,0]):

            if tmp[ee,1] == -1:
                continue
            bundles[elec_ind+1] = tmp[ee,1]

        self.bundles = bundles

        return self.bundles

    def set_bundles(self):
        """
        Sets to the object.
        """
        if self.bundles is None:
            _ = self.get_bundles()
        
        return None

    def compute_selectivity(self,exclude_postbundle=True):
        """
        Computes selectivity index for each parasol cell based on pre-bundle
        dictionary elements. Sets to the object

        if include_postbundle is true, probabilities after bundle threshold
        are included.
        """
        self.set_dictionary()
        self.set_bundles()
        self.set_cellids()
        self.set_celltypes_dict()
        parasols = []

        for cell in self.cellids:

            if "parasol" not in self.celltypes_dict[cell]:
                continue

            parasols.append(cell)
        
        # Initialize selectivity object.
        dictionary = self.dictionary['dictionary_tensor'].copy()

        assert len(parasols) == dictionary.shape[0]
        selectivity_dict = dict()

        # Trim the dictionary based on bundle.
        if exclude_postbundle:
            dictionary_prebundle = np.zeros(dictionary.shape)

            for ee in range(dictionary.shape[1]):
                elec = ee + 1 # index to electrode

                if elec not in self.bundles or self.bundles[elec] == -1:
                    dictionary_prebundle[:,ee,:] = dictionary[:,ee,:]
                else:
                    bundle_thr = self.bundles[elec] 
                    dictionary_prebundle[:,ee,0:bundle_thr] =\
                                                dictionary[:,ee,0:bundle_thr]
            
            dictionary = dictionary_prebundle
        else:
            pass

        # Squash anything with low probability (spontaneous) to 0.
        dictionary[np.where(np.amax(dictionary,axis=-1) < 
                            cfg.MIN_SPIKE_PROB_FIT)] = 0
            
        # For each cell, compute selectivity metrics. 
        for cc in range(dictionary.shape[0]):
            cell_probs = dictionary[cc,...][None,...]
            other_inds = np.setdiff1d(np.arange(0,len(parasols)),
                                      np.array([cc]))
            other_probs = dictionary[other_inds,...]
            '''
            Take p_target * (1 - p_nontarget): max over amps, min over cells,
            max over electrodes.
            '''
            selectivity = cell_probs * (1 - other_probs)
            max_selectivity = np.max(
                                np.min(
                                    np.max(
                                        selectivity,
                                        axis=-1
                                    ),
                                    axis=0)
                              )   
            selectivity_dict[parasols[cc]] = max_selectivity
            '''
            cell_probs = np.max(dictionary[cc,...],
                                 axis=-1)[None,:,None]
            max_inds = np.argmax(dictionary[cc,...],axis=-1)
            other_probs = np.take_along_axis(
                                    dictionary,
                                    max_inds[None,:,None],axis=-1
                                  )
            selectivity = (cell_probs * (1 - other_probs)).squeeze()

            # Take the minimum over cells and maximum over electrodes.
            other_inds = np.setdiff1d(np.arange(0,len(parasols)),
                                      np.array([cc]))
            max_selectivity = np.max(np.min(selectivity[other_inds,...],
                                     axis=0))
            selectivity_dict[parasols[cc]] = max_selectivity
            '''
    
        self.selectivity_dict = selectivity_dict

        return None

    def compute_recon_filters(self,write=True):
        """
        Wrapper to get reconstruction filters for a range of stixel sizes.
        """

        # Loop over parasols and get Gaussian fits.
        self.set_cellids()
        self.set_celltypes_dict() 
        self.set_sta_spatial_fit()
        self.set_vcd()
        parasols = []
        gaussian_popt = []
        stixel_size = int(self.vcd.runtimemovie_params.pixelsPerStixelX)

        for cc,cell in enumerate(self.cellids):
            celltype = self.celltypes_dict[cell]

            if "parasol" not in celltype:
                continue
                
            parasols.append(cell)
            gaussian_popt.append(self.sta_spatial_fit[cc])
        
        self.recon_filter_dict = recon.get_filter_dict(
                                        parasols,
                                        gaussian_popt,
                                        stixel_size
                                    )
        
        if write:
            fnameout = os.path.join(
                        cfg.PARENT_RECON,
                        cfg.RECON_FNAME%self.dataset
                    )
            np.save(fnameout,self.recon_filter_dict,allow_pickle=True)
        
        return None

    def get_recon_filters(self):
        """
        Loads from disk.
        """
        if self.recon_filter_dict is not None:
            return self.recon_filter_dict

        fnamein = os.path.join(
                cfg.PARENT_RECON,
                cfg.RECON_FNAME%self.dataset
              )

        try: 
            self.recon_filter_dict = np.load(fnamein,
                                             allow_pickle=True).item()
        except:
            print(f"could not filters for {self.dataset},\
                                         {self.estim_datarun}")
            return None

        return self.recon_filter_dict

    def set_recon_filters(self):
        """
        Sets to object.
        """
        if self.recon_filter_dict is None:
            _ = self.get_recon_filters()

        return None

    def compute_recon_filters_ns(self,write=True):
        """
        Wrapper to get reconstruction filters from simulated responses to
        natural scenes stimulation.
        """

        # Loop over parasols and get Gaussian fits.
        self.set_cellids()
        self.set_celltypes_dict() 
        self.set_sta_spatial_fit()
        self.set_nonlinearity_fit()
        self.set_vcd()
        stixel_size = int(self.vcd.runtimemovie_params.pixelsPerStixelX)
        parasols = []
        gaussian_popt = []
        nl_popt = []

        for cc,cell in enumerate(self.cellids):
            celltype = self.celltypes_dict[cell]

            if "parasol" not in celltype:
                continue
                
            parasols.append(cell)
            gaussian_popt.append(self.sta_spatial_fit[cc])
            nl_popt.append(self.nonlinearity_fit['popt'][cc])
        
        self.recon_filter_dict_ns = recon.get_filter_dict_ns(
                                                parasols,
                                                gaussian_popt,
                                                nl_popt,
                                                stixel_size
                                            )
        if write:
            fnameout = os.path.join(
                        cfg.PARENT_RECON,
                        cfg.RECON_NS_FNAME%self.dataset
                    )
            np.save(fnameout,self.recon_filter_dict_ns,allow_pickle=True)
        
        return None

    def get_recon_filters_ns(self):
        """
        Loads from disk.
        """
        if self.recon_filter_dict_ns is not None:
            return self.recon_filter_dict_ns

        fnamein = os.path.join(
                cfg.PARENT_RECON,
                cfg.RECON_NS_FNAME%self.dataset
              )

        try: 
            self.recon_filter_dict_ns = np.load(fnamein,
                                                allow_pickle=True).item()
        except:
            print(f"could not load filters for {self.dataset},\
                                         {self.estim_datarun}")
            return None

        return self.recon_filter_dict_ns

    def set_recon_filters_ns(self):
        """
        Sets to object.
        """
        if self.recon_filter_dict_ns is None:
            _ = self.get_recon_filters_ns()

        return None

    def make_cnn_data_split(self,crop=True,write=True):
        """
        Makes training/test split for the CNN denoiser on top of linear
        reconstruction. Simulates responses to a held out 10,000 set of images,
        gets the linear reconstruction, and writes both the reconstructed and 
        original to disk. Crops if indicated in a region that gets rid of zeros.
        """
        self.set_recon_filters_ns()
        stimulus_train = io.get_raw_movie(
                        os.path.join(cfg.IMAGENET_PARENT,
                        cfg.IMAGENET_TRAIN_CNN),
                        cfg.N_TRAIN_STIMULI
                    )
        '''
        stimulus_train = np.asarray([cv2.resize(stimulus_train[i,...],
                        dsize=(cfg.RESAMPLE_FIELD_X,cfg.RESAMPLE_FIELD_Y)
                        )
                        for i in range(stimulus_train.shape[0])
                ])
        '''
        stimulus_train = stimulus_train[...,0]
        stimulus_train = stimulus_train / 255
        stimulus_train -= np.mean(stimulus_train.ravel())

        stimulus_test = io.get_raw_movie(
                        os.path.join(cfg.IMAGENET_PARENT,
                        cfg.IMAGENET_TEST_CNN),
                        cfg.N_TEST_STIMULI_NS
                    )
        '''
        stimulus_test = np.asarray([cv2.resize(stimulus_test[i,...],
                        dsize=(cfg.RESAMPLE_FIELD_X,cfg.RESAMPLE_FIELD_Y)
                        )
                        for i in range(stimulus_test.shape[0])
                ])
        '''
        stimulus_test = stimulus_test[...,0]
        stimulus_test =  stimulus_test / 255
        stimulus_test -= np.mean(stimulus_test.ravel())

        n_stimuli,n_pixels_y,n_pixels_x = stimulus_train.shape
        encoder = self.recon_filter_dict_ns['encoding_filters']
        decoder = self.recon_filter_dict_ns['decoding_filters']

        # Center the filters to the middle of the visual field.
        n_cells,n_pixels = decoder.shape
        decoder = np.reshape(decoder,(n_cells,n_pixels_y,n_pixels_x))
        decoder = recon.center_filters(decoder,n_pixels_y,n_pixels_x)
        decoder = np.reshape(decoder,(n_cells,n_pixels))

        # Same for encoder.
        encoder = np.reshape(encoder,(n_cells,n_pixels_y,n_pixels_x))
        encoder = recon.center_filters(encoder,n_pixels_y,n_pixels_x)
        encoder = np.reshape(encoder,(n_cells,n_pixels))

        # Simulate responses and reconstruct.
        stimulus_train = np.reshape(stimulus_train,(n_stimuli,n_pixels))
        responses = np.maximum(stimulus_train@encoder.T,0)
        decoded_stimuli_train = responses@decoder

        n_stimuli_test,n_pixels_y,n_pixels_x = stimulus_test.shape
        stimulus_test = np.reshape(stimulus_test,(n_stimuli_test,n_pixels))
        responses = np.maximum(stimulus_test@encoder.T,0)
        decoded_stimuli_test = responses@decoder

        # Crop if indicated and write out.
        decoded_stimuli_train = np.reshape(
                            decoded_stimuli_train,
                            (n_stimuli,n_pixels_y,n_pixels_x)
                         )
        stimulus_train = np.reshape(
                         stimulus_train,
                        (n_stimuli,n_pixels_y,n_pixels_x)
                        )

        decoded_stimuli_test = np.reshape(
                            decoded_stimuli_test,
                            (n_stimuli_test,n_pixels_y,n_pixels_x)
                         )
        stimulus_test = np.reshape(
                         stimulus_test,
                        (n_stimuli_test,n_pixels_y,n_pixels_x)
                 )

        if crop:
            decoded_stimuli_train = decoded_stimuli_train[
                                ...,cfg.IMAGENET_CROP_X1:cfg.IMAGENET_CROP_X2
                                ]
            stimulus_train = stimulus_train[
                        ...,cfg.IMAGENET_CROP_X1:cfg.IMAGENET_CROP_X2
                      ]

            decoded_stimuli_test = decoded_stimuli_test[
                                ...,cfg.IMAGENET_CROP_X1:cfg.IMAGENET_CROP_X2
                                ]
            stimulus_test = stimulus_test[
                        ...,cfg.IMAGENET_CROP_X1:cfg.IMAGENET_CROP_X2
                      ]

        if write:

            if not os.path.isdir(cfg.PARENT_CNN):
                os.makedirs(cfg.PARENT_CNN)

            # Sloppy but overwrite original stimuli just in case.
            fnameout = os.path.join(
                        cfg.PARENT_CNN,
                        cfg.CNN_TRAIN_Y_FNAME    
            )
            np.save(fnameout,stimulus_train)

            fnameout = os.path.join(
                        cfg.PARENT_CNN,
                        cfg.CNN_TEST_Y_FNAME    
            )
            np.save(fnameout,stimulus_test)

            # Write data set specific reconstructions.
            fnameout = os.path.join(
                        cfg.PARENT_CNN,
                        cfg.CNN_TRAIN_X_FNAME%self.dataset
                       )
            np.save(fnameout,decoded_stimuli_train)

            fnameout = os.path.join(
                        cfg.PARENT_CNN,
                        cfg.CNN_TEST_X_FNAME%self.dataset
                       )
            np.save(fnameout,decoded_stimuli_test)
        
        return None
            



    def get_greedy_results(self,recon_type):
        """
        Loads from disk.
        """
        assert recon_type in ['wn','ns','ch'],"unknown recon type"

        if recon_type in ["wn"]:
            greedy_fname = cfg.GREEDY_FNAME
        elif recon_type in ["ns"]:
            greedy_fname = cfg.GREEDY_FNAME_NS
        else:
            greedy_fname = cfg.GREEDY_FNAME_CH

        if self.greedy_results is not None:
            return self.greedy_results

        fnamein = os.path.join(cfg.PARENT_GREEDY,
                         greedy_fname%self.dataset
                    )

        try:
            self.greedy_results = np.load(fnamein,allow_pickle=True).item()
        except:
            print(f"could not load greedy results for {self.dataset},\
                                         {self.estim_datarun}")
            return None

        return self.greedy_results

    def set_greedy_results(self,recon_type):
        """
        Sets to the object. Overwrites existing field because of the multiple 
        options.
        """
        assert recon_type in ['wn','ns','ch'],"unknown recon type"

        _ = self.get_greedy_results(recon_type)

        return None

    def set_pruned_dictionary_matrix(self,filter_type='wn',
                                     exclude_postbundle=True):
        """
        Prunes dictionary based on bundle thresholds and the cells available
        for reconstruction. Also sets any part of the sigmoid less than 
        EPSILON to 0 (for numerical reasons).

        Sets to the object, doens't return
        """
        self.set_cellids()
        self.set_celltypes_dict()
        self.set_bundles()
        self.set_dictionary()
        parasols = []

        for cell in self.cellids:

            if "parasol" not in self.celltypes_dict[cell]:
                continue

            parasols.append(cell)

        # Load the appropriate dictionary
        assert filter_type in ["ns","wn"],"wrong dictionary type"

        if filter_type in ["ns"]:
            self.set_recon_filters_ns()
            recon_filters = self.recon_filter_dict_ns
        else:
            self.set_recon_filters()
            recon_filters = self.recon_filter_dict

        # Prune by the available cells. 
        cells_reconstruction = recon_filters['cells_reconstruction']
        reconstruction_inds = []

        for cc,cell in enumerate(parasols):

            if cell in cells_reconstruction:
                reconstruction_inds.append(cc)
        
        reconstruction_inds = np.asarray(reconstruction_inds)
        tmp = self.dictionary['dictionary_tensor']

        dictionary = tmp[reconstruction_inds,...]
        del tmp

        # Prune by bundle if indicated.
        if exclude_postbundle:
            dictionary_prebundle = np.zeros(dictionary.shape)

            for ee in range(dictionary_prebundle.shape[1]):
                elec = ee + 1

                if elec not in self.bundles or self.bundles[elec] == -1:
                    dictionary_prebundle[:,ee,:] = dictionary[:,ee,:]
                else:
                    bundle_thr = self.bundles[elec] 
                    dictionary_prebundle[:,ee,0:bundle_thr] =\
                                                dictionary[:,ee,0:bundle_thr]
            
            dictionary = dictionary_prebundle
        else:
            pass


        # Squash anything with low probability (spontaneous) to 0.
        dictionary[np.where(np.amax(dictionary,axis=-1) < 
                            cfg.MIN_SPIKE_PROB_FIT)] = 0

        # Reshape into a matrix.
        n_cells,n_elecs,n_amps = dictionary.shape

        dictionary_matrix = np.reshape(
                                dictionary,
                                (n_cells,n_elecs * n_amps)
                        ).T

        # Set anything below EPSILON to 0.
        dictionary_matrix[np.where(dictionary_matrix < cfg.EPSILON)] = 0
        ''''
        n_elems,n_cells = dictionary_matrix.shape
        tmp = dictionary_matrix.ravel()
        tmp[np.argwhere(tmp < cfg.EPSILON).flatten()] = 0
        dictionary_matrix = tmp.reshape((n_elems,n_cells))


        del tmp
        '''

        # Remove all zero elements and then add a single one back.
        all_inds = np.arange(0,dictionary_matrix.shape[0])
        zero_inds = np.argwhere(np.sum(
                                dictionary_matrix,axis=1
                                ) == 0).flatten()
        nonzero_inds = np.setdiff1d(all_inds,zero_inds)
    
        if zero_inds.shape[0] != 0:
            dictionary_matrix = dictionary_matrix[nonzero_inds,:]
        
        dictionary_matrix = np.c_[dictionary_matrix.T,
                                  np.zeros(dictionary_matrix.shape[1])
                                ].T 
        self.pruned_dictionary = dictionary_matrix

        return None