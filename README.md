McGhee–von Hippel Binding Analysis Module

This module implements the McGhee–von Hippel model for cooperative ligand–DNA binding.
It allows you to fit experimental titration data (absorbance, fluorescence, etc.) and extract binding parameters such as:

Binding constant (K)

Site size (n)

Cooperativity parameter (ω) — optional, for cooperative models

The model is widely used to describe intercalative or groove-binding dyes interacting with DNA or polymers.

🧪 Features

Non-linear regression fitting of binding curves (r vs r_Cf)

Supports both non-cooperative and cooperative binding models

Automatic calculation of:

Best-fit parameters and their standard errors

Coefficient of determination (R²)

Optional constraints and bounds

Config file support (config.toml) for streamlined input and reproducibility

CSV import/export compatibility with UVTools

Publication-ready plots (binding isotherm and residuals)

⚙️ Configuration File

The analysis parameters are defined in a config.toml file, allowing reproducible, script-free execution.
An example config.toml:

filepath = "titration_data.csv"
wl = 506
conc_DNA_M = 1.15e-3
conc_dye_M = 44e-6
initial_volume_L = 500e-6
p0 = [1e6, 2.0]
bounds = [[1e3, 0.5], [1e8, 10.0]]
function = "mvh_rhs"


You can easily modify this file to run different analyses without editing the code.

🚀 Basic Workflow
from mvh_analysis.mghee_von_hippel import *
from mvh_analysis.config import config_titration
from pathlib import Path

# === Load Configuration ===
config_path = Path("config/config.toml")
cfg = config_titration(config_path)

# === Unpack Parameters ===
filepath = cfg[0]
wl = cfg[1]
conc_DNA_M = cfg[2]
conc_dye_M = cfg[3]
initial_volume_L = cfg[4]
p0 = cfg[5]
bounds = cfg[6]
function = getattr(mvh_analysis.mghee_von_hippel, cfg[7])

# === Run Analysis ===
df = prepare_data(filepath, wl, conc_DNA_M, initial_volume_L, conc_dye_M)
K, n = fit_mvh(function, df["r"], df["r_Cf"], p0=p0, bounds=bounds)
r2 = calculate_R_squared(function, df["r"], df["r_Cf"], [K, n])

print(f"K = {K:.2e} M⁻¹")
print(f"n = {n:.2f}")
print(f"R² = {r2:.4f}")

🧩 Function Overview
Function	Description
prepare_data(filepath, wl, conc_DNA_M, initial_volume_L, conc_dye_M)	Loads titration data, applies dilution correction, and computes derived quantities.
select_wl(data, wl)	Extracts absorbance vs. volume data for a specific wavelength from a JASCO-like dataset.
mvh_prep(df, conc_DNA_M, initial_volume_L, conc_dye_M)	Calculates DNA and ligand concentrations, fraction bound, and binding ratios.
mvh_rhs(r, K, n)	Returns theoretical McGhee–von Hippel binding curve (r/Cf) for given parameters.
fit_mvh(function, r, r_Cf, p0, bounds)	Performs non-linear regression to extract best-fit K and n using SciPy.
calculate_R_squared(function, x, y, param)	Calculates R² to assess the goodness of fit.
📁 Expected Input

Your CSV file should contain absorbance values as a function of added DNA volume:

Wavelength (nm)	0.0	1.0	2.0	3.0	...
506.0	0.602741	0.544984	0.560099	0.552854	...

The first column lists wavelengths, and each subsequent column corresponds to a titration point (volume of DNA added, in µL or L).
The script automatically selects the target wavelength and reshapes the dataset.

📊 Example Output
K = 1.26e+06 M⁻¹
n = 1.66
R² = 0.9436


The returned DataFrame includes all computed quantities (r, Free_ligand, Bound_ligand, etc.) for further analysis.

📈 Optional Plotting
import matplotlib.pyplot as plt

plt.scatter(df["r"], df["r_Cf"], label="Experimental")
plt.plot(df["r"], mvh_rhs(df["r"], K, n), "r-", label="Fit")
plt.xlabel("r")
plt.ylabel("r/Cf")
plt.legend()
plt.tight_layout()
plt.show()

📚 Reference

McGhee, J. D.; von Hippel, P. H.
Theoretical Aspects of DNA–Protein Interactions: Cooperative and Non-cooperative Binding of Large Ligands to a One-Dimensional Homogeneous Lattice.
J. Mol. Biol. 1974, 86, 469–489.