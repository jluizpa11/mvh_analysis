# -*- coding: utf-8 -*-
"""
Created on Fri Oct  3 18:51:59 2025

@author: JL
"""
import pandas as pd
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import mvh_analysis.mghee_von_hippel as mvh

def prepare_data(filepath, wl, conc_DNA_M, initial_volume_L, conc_dye_M ):
    """
    Load and preprocess titration data for McGhee–von Hippel analysis.

    This function reads a CSV file containing titration data (e.g., absorbance vs. volume),
    selects the wavelength of interest, and prepares a tidy DataFrame containing
    all quantities necessary for McGhee–von Hippel binding analysis.

    The function performs:
    - Data loading
    - Wavelength selection
    - Data transformation and column labeling
    - Preparation of binding variables via `mvh_prep()`

    Parameters
    ----------
    filepath : str or Path
        Path to the CSV file containing the titration data.

    Returns
    -------
    pandas.DataFrame
        A processed DataFrame containing:
        - Volume (L)
        - Absorbance (corrected)
        - Fraction_Bound
        - r (bound sites per DNA base pair)
        - Free_ligand (M)
        - Bound_ligand (M)
        - r_Cf (binding ratio used for McGhee–von Hippel fit)
    """

    # I/O csv
    data = pd.read_csv(filepath)

    # Select wavelength 
    data_wl = select_wl(data, wl)

    # Prep wl
    df = mvh_prep(data_wl, conc_DNA_M, initial_volume_L, conc_dye_M)
    
    return df

def mvh_prep(df, conc_DNA_M, initial_volume_L, conc_dye_M):
    """
    Prepare the dataset for McGhee–von Hippel binding analysis.

    Given a titration dataset (Absorbance vs. Volume), this function calculates
    all derived quantities needed for McGhee–von Hippel model fitting:
    DNA concentration, dilution correction, fraction bound, free and bound ligand,
    and binding ratio (r/Cf).

    Parameters
    ----------
    df : pandas.DataFrame
        Input DataFrame containing at least 'Volume' and 'Absorbance' columns.
    conc_DNA_M : float
        Concentration of DNA stock solution (M).
    initial_volume_L : float
        Initial volume of the cuvette or solution before titration (L).
    conc_dye_M : float
        Initial concentration of dye in the cuvette (M).

    Returns
    -------
    pandas.DataFrame
        DataFrame with additional columns:
        - conc_DNA : DNA concentration at each titration point (M)
        - Absorbance : corrected absorbance values
        - Fraction_Bound : fraction of ligand bound to DNA
        - r : ratio of bound sites to DNA base pairs
        - Free_ligand : free dye concentration (M)
        - Bound_ligand : bound dye concentration (M)
        - r_Cf : ratio of r to free ligand concentration
    """
    
    # Convert volume DNA to [DNA]    
    df["conc_DNA"] = conc_DNA_M * df["Volume"] / (initial_volume_L + df["Volume"])
    
    #Correct the absorbance
    dilution_factor = (df["Volume"] + initial_volume_L) / initial_volume_L
    df["Absorbance_pre_correction"] = df["Absorbance"]
    df["Absorbance"] = df["Absorbance_pre_correction"] * dilution_factor
    
    #Correct the Dye Concentration
    
    dye_corrected = dilution_factor * conc_dye_M

    # Compute Fraction bound
    A_free = df["Absorbance"].iloc[0]
    A_bound = df["Absorbance"].min()
    df["Fraction_Bound"] =  (df["Absorbance"] - A_free) / (A_bound - A_free)
    
    # Compute bound ligand
    df["Bound_ligand"] = dye_corrected * df["Fraction_Bound"]
    
    # Compute free ligand
    df["Free_ligand"] = dye_corrected -  df["Bound_ligand"]
    
    # Compute r
    df["r"] = df["Bound_ligand"] / df["conc_DNA"]
    
    # Compute r/Cf
    df["r_Cf"] = df["r"] / df["Free_ligand"]
    
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["r", "Free_ligand", "r_Cf"])

    return df

def mvh_rhs(r, K, n):
    """
    McGhee–von Hippel binding isotherm (non-cooperative case).

    Computes the theoretical right-hand side of the McGhee–von Hippel equation
    for a given set of binding parameters.

    Equation:
        r/Cf = K * (1 - n*r) * ((1 - n*r)/(1 - (n - 1)*r))**(n - 1)

    Parameters
    ----------
    r : array-like
        Fraction of occupied binding sites.
    K : float
        Intrinsic binding constant (M⁻¹).
    n : float
        Number of lattice sites (base pairs) occupied per ligand molecule.

    Returns
    -------
    numpy.ndarray
        Computed r/Cf values for each r, given K and n.
    """
    # McGhee–von Hippel (non-cooperative) RHS for r/Cf
    # r/Cf = K*(1 - n*r) * ((1 - n*r)/(1 - (n - 1)*r))**(n - 1)
    # Remove rows where any critical variable is NaN or inf
    eps = 1e-12

    num = np.clip(1 - n*r, eps, None)
    den = np.clip(1 - (n - 1)*r, eps, None)

    ratio = num / den
    y = K * num * np.power(ratio, n - 1)

    # Optional: prevent overflow
    y = np.clip(y, -1e300, 1e300)

    return y
    
    # num = (1 - n*r)
    # den = (1 - (n - 1)*r)
    # # guard small denominators
    # eps = 1e-12
    # den = np.clip(den, eps, None)
    # return K * num * (num/den)**(n - 1)

def fit_mvh(function, r, r_Cf, p0, bounds):
    
    """
    Fit the McGhee–von Hippel binding model to experimental data.

    Uses nonlinear least-squares regression (SciPy's `curve_fit`)
    to optimize K and n parameters based on experimental r and r/Cf data.

    Parameters
    ----------
    function : callable
        Theoretical function to fit (e.g., `mvh_rhs`).
    r : array-like
        Experimental fractional occupancy (r).
    r_Cf : array-like
        Experimental ratio of r to free ligand concentration.
    p0 : sequence of float
        Initial guesses for parameters [K, n].
    bounds : 2-tuple of array-like
        Lower and upper bounds for [K, n].

    Returns
    -------
    tuple
        (K_fit, n_fit) : optimized values of K and n.
    """
    # Fit the curve
    
    popt, pcov =  curve_fit(function, r, r_Cf, p0=p0, bounds=bounds)
    perr = np.sqrt(np.diag(pcov))
    
    K_fit, n_fit = popt
    K_err, n_err = perr
    
    return K_fit, n_fit, K_err, n_err

def calculate_R_squared(function, x, y, param):
    
    """
    Compute the coefficient of determination (R²) for a fitted model.

    Calculates how well the model describes the experimental data
    by comparing the residual sum of squares (SS_res) to the total sum of squares (SS_tot).

    Parameters
    ----------
    function : callable
        Theoretical model used for fitting (e.g., `mvh_rhs`).
    x : array-like
        Experimental x-values (e.g., r).
    y : array-like
        Experimental y-values (e.g., r/Cf).
    param : list or tuple
        Fitted parameters (e.g., [K, n]).

    Returns
    -------
    float
        Coefficient of determination (R²), ranging from 0 to 1.
    """
    
    # Predicted values from your fitted model
    y_pred = function(x, param[0], param[1])
    
    # Compute residuals and R²
    ss_res = np.sum((y - y_pred)**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    r_squared = 1 - (ss_res / ss_tot)
    
    return r_squared

def mvh_analysis(df,cfg):
    
    print("------------------------------")
    print("Initializing McGhee - von Hippel Analysis")
    
    import numpy as np

    r = df["r"].values
    r_Cf = df["r_Cf"].values    
    r = df["r"]
    r_Cf = df["r_Cf"]
    p0 = cfg[5]
    bounds = cfg[6]
    function = mvh_rhs
    # assuming p0 = [1e6, 1.5]
    y_pred = function(r, *p0)

    try:
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
        
        print("McGhee - von Hippel Analysis complete")
        
    except Exception as e:
            
        print(f"Error during analysis. Error: {e}")
        
        plt.figure(figsize=(8,6))
        plt.plot(df["r"], df["r_Cf"],'o', label="Data")
        plt.xlabel("r")
        plt.ylabel("r/Cf")
        plt.title("McGhee-von Hippel fit")
        plt.legend()
        plt.grid()
        plt.tight_layout()
        plt.show()
        
        print("Plotted data only.")