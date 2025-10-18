# -*- coding: utf-8 -*-
"""
Created on Sun Oct  5 20:59:10 2025

@author: JL
"""
from mvh_analysis.mghee_von_hippel import *
from mvh_analysis.config import config_titration, load_config
import mvh_analysis.mghee_von_hippel as mvh
import numpy as np

def main_df(config_path):
    
    cfg = config_titration(config_path)
    ############## Variables  ######################

    filepath = cfg[0]  
    wl = cfg[1]
    conc_DNA_M = cfg[2]
    conc_dye_M = cfg[3]
    initial_volume_L = cfg[4]
    
    # initial guesses: K in 1e4–1e6 M^-1 typical; n around 2–3 for intercalators
    p0 = cfg[5]
    bounds = cfg[6]
    function = getattr(mvh, cfg[7])

    ###########################################
    
    df = prepare_data(filepath, wl, conc_DNA_M, initial_volume_L, conc_dye_M)

    A_min_index = df["Absorbance"].idxmin()
    
    df = df.loc[:A_min_index].copy()
    df = df.drop(13)
    df = df.drop(14)
    df = df.drop(11)

    r = df["r"]
    r_Cf = df["r_Cf"]

    K_fit, n_fit = fit_mvh(function,r, r_Cf, p0, bounds)

    param = [K_fit, n_fit]
    r_squared = calculate_R_squared(function, r, r_Cf, param)

    # Optional: visualize fit over your r-range
    r_grid = np.linspace(df["r"].min(), df["r"].max(), 200)
    y_fit = mvh_rhs(r_grid, K_fit, n_fit)

    print(f"K = {K_fit:.2e} M^-1")
    print(f"n = {n_fit:.2f}")
    print(f"R² = {r_squared:.4f}")

    ##########################        PLOT         ##############################################

    plt.figure(figsize=(8,6))
    plt.plot(df["r"], df["r_Cf"],'o', label="Data")
    plt.plot(r_grid, y_fit,'-', label=f"Fit\nK = {K_fit:.2e} M^-1\nn = {n_fit:.2f}\nR2 = {r_squared:.4f}")
    plt.xlabel("r")
    plt.ylabel("r/Cf")
    plt.title("McGhee-von Hippel fit")
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.show()
    
    return None

def compute_mols_added_DNA(df):
    df["mols_added_DNA"] = df["Volume_added"] * 1e-6 * df["Concentration_stock_DNA"] * 1e-6
    
    return df

def compute_mols_total_DNA(df):
    

    df["mols_total_DNA"] = df["mols_added_DNA"].cumsum()
    return df

def compute_total_volume(df, initial_volume_L):
    
    df["volume_final"] = df["Volume"] + initial_volume_L
    
    return df

def compute_conc_DNA(df):
    
    df["conc_DNA"] = df["mols_total_DNA"] / df["volume_final"]
    
    return df

def correct_abs(df, initial_volume_L):
    
    #Correct the absorbance
    dilution_factor = (df["Volume"] + initial_volume_L) / initial_volume_L
    df["Absorbance_pre_correction"] = df["Absorbance"]
    df["Absorbance"] = df["Absorbance_pre_correction"] * dilution_factor
    
    return df

def compute_fraction_bound(df):
    
    # Compute Fraction bound
    A_free = df["Absorbance"].iloc[0]
    A_bound = df["Absorbance"].min()
    df["Fraction_Bound"] =  (df["Absorbance"] - A_free) / (A_bound - A_free)
    
    return df

def compute_bound_ligand(df, dye_corrected):
    
    # Compute bound ligand
    df["Bound_ligand"] = dye_corrected * df["Fraction_Bound"]
    
    return df

def compute_free_ligand(df, dye_corrected):
    
    # Compute free ligand
    df["Free_ligand"] = dye_corrected -  df["Bound_ligand"]
    
    return df
    
def compute_r(df):
    
    # Compute r
    df["r"] = df["Bound_ligand"] / df["conc_DNA"]
    
    return df
    
def compute_r_Cf(df):
    # Compute r/Cf
    df["r_Cf"] = df["r"] / df["Free_ligand"]
    
    return df