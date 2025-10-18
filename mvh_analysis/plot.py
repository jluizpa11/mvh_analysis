# -*- coding: utf-8 -*-
"""
Created on Sun Oct 12 21:08:09 2025

@author: JL
"""

import matplotlib.pyplot as plt
import numpy as np
from mvh_analysis.mghee_von_hippel import mvh_rhs 

def plot_mvh(df, K_fit, n_fit, K_err, n_err, r_squared):
    
    r_grid = np.linspace(df["r"].min(), df["r"].max(), 200)
    y_fit = mvh_rhs(r_grid, K_fit, n_fit)
    plt.figure(figsize=(8,6))
    plt.plot(df["r"], df["r_Cf"],'o', label="Data")
    plt.plot(r_grid, y_fit,'-', label=f"Fit\nK = {K_fit:.2e}±{K_err:.2e} M^-1\nn = {n_fit:.2f}±{n_err:.2f}\nR2 = {r_squared:.4f}")
    plt.xlabel("r")
    plt.ylabel("r/Cf")
    plt.title("McGhee-von Hippel fit")
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.show()
    
def plot_spectra(df, xlim=None, ylim=None):
    
    plt.figure(figsize=(8,6))
    for column in df.columns:
        if column != "Wavelength":
            plt.plot(df["Wavelength"], df[column],'-')
    plt.xlabel("Wavelength")
    plt.ylabel("Absorbance")
    plt.title("Plot Spectra")
    if xlim is not None:
        plt.xlim(xlim[0], xlim[1])
        
    if ylim is not None:
        plt.ylim(ylim[0], ylim[1])
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.show()