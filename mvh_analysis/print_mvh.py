# -*- coding: utf-8 -*-
"""
Created on Mon Oct 13 18:00:24 2025

@author: JL
"""

def print_initializing_mvh():
    print("Initializing McGhee - von Hippel Analysis")

def print_wl_selected(wl):
    
    print(f"Wavelength selected: {wl} nm")
    
def print_fit_r2_with_err(K_fit, n_fit, K_err, n_err, r2):
    
    print(f"K = {K_fit:.2e}±{K_err:.2e} M^-1")
    print(f"n = {n_fit:.2f}±{n_err:.2f}")
    print(f"R² = {r2:.4f}")
    
def print_analysis_complete():
    
    print("McGhee - von Hippel Analysis complete")
    
def print_spectra():
    
    print("Spectra plotted")