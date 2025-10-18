# -*- coding: utf-8 -*-
"""
Created on Fri Oct 10 18:46:57 2025

@author: JL
"""

from mvh_analysis.io import *
from mvh_analysis.mghee_von_hippel import mvh_analysis
from mvh_analysis.config import *
from mvh_analysis.data_prep import *
from mvh_analysis.config import config_titration
from pathlib import Path
import numpy as np


def prepare_full_df_after_wl(config_path):
    

    cfg = config_titration(config_path)
    # %% Define parameters
    
    file_path = cfg[0]
    df = load_raw_df(file_path)
    initial_volume_L = cfg[4]
    conc_dye_M = cfg[3]
    
    # %% Prepare initial df and run mvh_analysis
    
    df = compute_mols_added_DNA(df)
    df = compute_mols_total_DNA(df)
    df = compute_total_volume(df, initial_volume_L)
    df = compute_conc_DNA(df)
    df = correct_abs(df, initial_volume_L)
    # %%
    dilution_factor = (df["Volume"] + initial_volume_L) / initial_volume_L
    
    dye_corrected = dilution_factor * conc_dye_M
    
    df = compute_fraction_bound(df)
    df = compute_bound_ligand(df, dye_corrected)
    df = compute_free_ligand(df, dye_corrected)
    df = compute_r(df)
    df = compute_r_Cf(df)
    
    return df

def prepare_df_without_nan(df):
    
    df = df.copy()
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["r", "Free_ligand", "r_Cf"])
    
    return df

def select_wl(data, wavelength):
    
    """
    Select a single wavelength trace from a wavelength–titration dataset.

     This function extracts the absorbance values at a specified wavelength
     from a 2D dataset (wavelength vs. titration volume) and reshapes it
     into a two-column format suitable for further analysis.

     Parameters
     ----------
     data : pandas.DataFrame
        Raw dataset with wavelengths in the first column and absorbance
        values across titration points in subsequent columns.
    wl : float
        Target wavelength (in nm) to extract.

    Returns
    -------
    pandas.DataFrame
        A two-column DataFrame with:
        - Volume (converted to liters)
        - Absorbance (float)
    """
    # Select wavelength

    data_wl = data[data.iloc[:,0] == wavelength].copy()
    data_wl = data_wl.set_index(data_wl.columns[0]).T.reset_index()
    data_wl.rename(columns={"index": "Volume"}, inplace=True)
    data_wl.columns = ["Volume", "Absorbance"]
    data_wl["Volume"] = data_wl["Volume"].astype(float)
    data_wl["Absorbance"] = data_wl["Absorbance"].astype(float)
    data_wl["Volume"] = data_wl["Volume"] * 1e-6
    
    return data_wl

def prepare_full_df(df, initial_volume_L, conc_dye_M):
    df = df.copy()
    df = compute_mols_added_DNA(df)
    df = compute_mols_total_DNA(df)
    df = compute_total_volume(df, initial_volume_L)
    df = compute_conc_DNA(df)
    df = correct_abs(df, initial_volume_L)
    # %%
    dilution_factor = (df["Volume"] + initial_volume_L) / initial_volume_L
    
    dye_corrected = dilution_factor * conc_dye_M
    
    df = compute_fraction_bound(df)
    df = compute_bound_ligand(df, dye_corrected)
    df = compute_free_ligand(df, dye_corrected)
    df = compute_r(df)
    df = compute_r_Cf(df)
    
    return df