# -*- coding: utf-8 -*-
"""
Created on Sun Oct  5 17:56:59 2025

@author: JL
"""
import tomllib
from pathlib import Path

def load_config(config_path: str | Path = "mvh_config.toml") -> dict:
    """Load and parse TOML configuration for McGhee–von Hippel analysis."""
    path = Path(config_path)
    if not path.exists():
        raise FileNotFoundError(f"Missing configuration file: {path}")
    with open(path, "rb") as f:
        return tomllib.load(f)
    
def config_titration(config_path: str | Path = "mvh_config.toml"):
    """
    Load all configuration parameters from a TOML file and return
    them as ready-to-use variables for the McGhee–von Hippel titration analysis.
    
    Converts units (µM → M, µL → L) and bundles all relevant
    experiment, fit, cleaning, and plot parameters.
    """
    cfg = load_config(config_path)

    # --- Extract sections ---
    paths = cfg.get("paths", {})
    exp = cfg.get("experiment", {})
    fit = cfg.get("fit", {})
    clean = cfg.get("cleaning", {})
    plot_cfg = cfg.get("plot", {})
    template = cfg.get("template", {})

    # --- Convert units ---
    wl = exp.get("wavelength_nm")
    conc_DNA_M = exp.get("conc_DNA_uM", 0) * 1e-6
    conc_dye_M = exp.get("conc_dye_uM", 0) * 1e-6
    initial_volume_L = exp.get("initial_volume_uL", 0) * 1e-6

    # --- Prepare fit parameters ---
    p0 = fit.get("p0", [1e6, 1.5])
    bounds = (
        fit.get("bounds_lower", [1e2, 1.0]),
        fit.get("bounds_upper", [1e8, 10.0]),
    )

    # --- Misc sections ---
    filepath = Path(paths.get("filepath", ""))
    #drop_rows = clean.get("drop_rows", [])
    #plot_params = plot_cfg
    function = fit.get("function")
    template = Path(template.get("template",""))
    
    return (
        filepath,
        wl,
        conc_DNA_M,
        conc_dye_M,
        initial_volume_L,
        p0,
        bounds,
        function,
        template
        #drop_rows,
        #plot_params,
    )


