# -*- coding: utf-8 -*-
"""
Created on Sat Oct 11 21:47:59 2025

@author: JL
"""

from mvh_analysis.config import *
from mvh_analysis.io import *
from mvh_analysis.df_prep import *
from mvh_analysis.data_prep import *
from mvh_analysis.plot import plot_mvh, plot_spectra
from mvh_analysis.print_mvh import print_analysis_complete, print_fit_r2_with_err, print_initializing_mvh, print_wl_selected
from pathlib import Path 
import pandas as pd

class Titration:
    
    def __init__(self, config_path):
        
        cfg = config_titration(config_path)
        
        self.filepath = cfg[0]  
        self.wl = cfg[1]
        self.conc_dye_M = cfg[3]
        self.initial_volume_L = cfg[4]
        # initial guesses: K in 1e4–1e6 M^-1 typical; n around 2–3 for intercalators
        self.p0 = cfg[5]
        self.bounds = cfg[6]
        self.function = getattr(mvh, cfg[7])
        self.raw_df = load_raw_df(self.filepath)
        self.path_template = cfg[8]
        self.template = load_raw_df(self.path_template)
        print("Configuration file loaded sucessfully")
        
    def get_wavelength(self, wavelength=None):
        
        if wavelength is None:
            
            wavelength = self.wl
        else:
            self.wl = wavelength
        self.df_wl = select_wl(self.raw_df, wavelength)
    
    def get_df_with_template(self):
        
        self.df_with_template = self.template.join(self.df_wl).copy()
        
    def get_df_full(self):
        
        self.df_full = prepare_full_df(self.df_with_template ,self.initial_volume_L, self.conc_dye_M)
        
    def get_df_clean(self):
        
        self.df_full_clean = prepare_df_without_nan(self.df_full)
        
    def set_initial_df(self):
        
        self.initial_df = self.df_full_clean.copy()
    
    def set_df(self, df=None):
        
        if df == None:
            self.df = self.df_full_clean.copy()
        else:
            self.df = df.copy()
        
    def fit(self):
        
        self.K, self.n, self.K_err, self.n_err = fit_mvh(self.function, self.df["r"], self.df["r_Cf"], self.p0, self.bounds)
    
    def get_r2(self):
        
        self.r2 = calculate_R_squared(self.function, self.df["r"], self.df["r_Cf"], (self.K, self.n))
    
    def plot(self):
        
        plot_mvh(self.df, self.K, self.n, self.K_err, self.n_err, self.r2)
        
    def plot_spectra(self, xlim=None, ylim=None):
        
        plot_spectra(self.raw_df, xlim, ylim)
        
    def initialize_mvh_analysis(self, wavelength=None):
        
        print_initializing_mvh()
        self.get_wavelength(wavelength)
        print_wl_selected(self.wl)
        self.get_df_with_template()
        self.get_df_full()
        self.get_df_clean()
        self.set_initial_df()
        self.set_df()
        self.fit()
        self.get_r2()
        print_fit_r2_with_err(self.K, self.n, self.K_err, self.n_err, self.r2)
        self.plot()
        print_analysis_complete()
        
    def mvh_analysis(self):
        
        print_initializing_mvh()
        print_wl_selected(self.wl)
        self.fit()
        self.get_r2()
        print_fit_r2_with_err(self.K, self.n, self.K_err, self.n_err, self.r2)
        self.plot()
        print_analysis_complete()
        
    def initialize_mvh_analysis_after_wl(self):
        self.df_wl = self.raw_df
        print_initializing_mvh()
        print_wl_selected(self.wl)
        self.get_df_with_template()
        self.get_df_full()
        self.get_df_clean()
        self.set_df()
        self.fit()
        self.get_r2()
        print_fit_r2_with_err(self.K, self.n, self.K_err, self.n_err, self.r2)
        self.plot()
        print_analysis_complete()