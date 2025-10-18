# -*- coding: utf-8 -*-
"""
Created on Sun Oct  5 20:59:10 2025

@author: JL
"""
from mvh_analysis.mghee_von_hippel import *
from mvh_analysis.config import config_titration
import mvh_analysis.mghee_von_hippel as mvh

def main(config_path):
    
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

    print("Configuration file loaded sucessfully")
    ###########################################
    print("Preparing raw data")
    df = prepare_data(filepath, wl, conc_DNA_M, initial_volume_L, conc_dye_M)
    
    r = df["r"].values
    r_Cf = df["r_Cf"].values
    Error = False
    try:
        r = df["r"]
        r_Cf = df["r_Cf"]
        
        K_fit, n_fit = fit_mvh(function,r, r_Cf, p0, bounds)
    
        param = [K_fit, n_fit]
        r_squared = calculate_R_squared(function, r, r_Cf, param)
    
        # Optional: visualize fit over your r-range
        r_grid = np.linspace(df["r"].min(), df["r"].max(), 200)
        y_fit = mvh_rhs(r_grid, K_fit, n_fit)
        print("------------------------------")
        print("Fitting mvh curve:")
        print(f"K = {K_fit:.2e} M^-1")
        print(f"n = {n_fit:.2f}")
        print(f"R² = {r_squared:.4f}")
        print("------------------------------")
    except Exception as e:
        Error = True
        print(f"Error during analysis. Error:{e}")
    ##########################        PLOT         ##############################################

    plt.figure(figsize=(8,6))
    plt.plot(df["r"], df["r_Cf"],'o', label="Data")
    if not Error:
        plt.plot(r_grid, y_fit,'-', label=f"Fit\nK = {K_fit:.2e} M^-1\nn = {n_fit:.2f}\nR2 = {r_squared:.4f}")
    plt.xlabel("r")
    plt.ylabel("r/Cf")
    plt.title("McGhee-von Hippel Raw Data plot and fit")
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.show()
    print("Raw data plotted")
    print("------------------------------")
    print("Finished raw data plotting and analysis.")
    
    return df, cfg