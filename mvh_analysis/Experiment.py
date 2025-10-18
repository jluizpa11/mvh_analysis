# -*- coding: utf-8 -*-
"""
Created on Sun Oct 12 21:53:02 2025

@author: JL
"""

from pathlib import Path
import pandas as pd
import numpy as np
from mvh_analysis.Titration import Titration
from mvh_analysis.mghee_von_hippel import mvh_analysis, fit_mvh, calculate_R_squared
from mvh_analysis.print_mvh import print_initializing_mvh, print_analysis_complete, print_wl_selected, print_fit_r2_with_err
from mvh_analysis.plot import plot_mvh
import matplotlib.pyplot as plt

class Experiment:
    def __init__(self, config_paths):
        self.titrations = [Titration(p) for p in config_paths]

    def initialize_all(self, wavelength=None):
        print("Initializing all titrations...")
        for t in self.titrations:
            t.initialize_mvh_analysis(wavelength)
    
    def run_all(self):
        print("Running all titrations...")
        for t in self.titrations:
            t.mvh_analysis()
            
    def summarize(self):
        """Create a summary DataFrame of all fitted results."""
        self.summary_df = pd.DataFrame([
            {
                "Name": Path(t.filepath).stem,
                "Wavelength (nm)": t.wl,
                "K (M^-1)": t.K,
                "n": t.n,
                "R²": t.r2,
            }
            for t in self.titrations
        ])
        self.summary_df["log K"] = np.log10(self.summary_df["K (M^-1)"])
        print("\n=== Summary of All Titrations ===")
        print(self.summary_df)
        return self.summary_df

    def average_results(self):
        """Compute mean and standard deviation of K, n, and R² across replicates."""
        if not hasattr(self, "summary_df"):
            self.summarize()

        df = self.summary_df

        self.avg_K = df["K (M^-1)"].mean()
        self.std_K = df["K (M^-1)"].std()
        self.avg_n = df["n"].mean()
        self.std_n = df["n"].std()
        self.avg_R2 = df["R²"].mean()

        print("\n=== Averaged Results Across Replicates ===")
        print(f"K = {self.avg_K:.2e} ± {self.std_K:.2e} M^-1")
        print(f"n = {self.avg_n:.2f} ± {self.std_n:.2f}")
        print(f"R² = {self.avg_R2:.4f}")

        return {
            "K_mean": self.avg_K,
            "K_std": self.std_K,
            "n_mean": self.avg_n,
            "n_std": self.std_n,
            "R2_mean": self.avg_R2
        }

    def plot_all(self, show_fit=True, title="All Titrations"):
        """
        Overlay all titrations' r vs r/Cf data (and fits, if available).

        Parameters
        ----------
        show_fit : bool, optional
            If True, plots both experimental data and the fitted curve.
            If False, plots only experimental data points.
        title : str, optional
            Title for the plot.
        """
        plt.figure(figsize=(8, 6))

        for t in self.titrations:
            # Check if the titration has a fitted DataFrame
            if not hasattr(t, "df"):
                print(f"Skipping {Path(t.filepath).stem}: no data found.")
                continue

            # Plot the experimental data points
            plt.plot(
                t.df["r"], t.df["r_Cf"], "o", markersize=6,
                label=f"{Path(t.filepath).stem} data"
            )

            # Optionally overlay the fitted curve
            if show_fit and hasattr(t, "K") and hasattr(t, "n"):
                r_grid = np.linspace(t.df["r"].min(), t.df["r"].max(), 200)
                y_fit = t.function(r_grid, t.K, t.n)
                plt.plot(
                    r_grid, y_fit, "-",
                    label=f"{Path(t.filepath).stem} fit\nK={t.K:.2e}, n={t.n:.2f}"
                )

        plt.xlabel("r", fontsize=12)
        plt.ylabel("r/Cf", fontsize=12)
        plt.title(title, fontsize=14)
        plt.legend(fontsize=9, frameon=False)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()

    def prepare_df_concat(self, exclude=[]):
        self.titrations_dfs = []
        
        for index, t in enumerate(self.titrations):
            if index not in exclude:
                self.titrations_dfs.append(t.df)
            
        self.df_concat = pd.concat(self.titrations_dfs, ignore_index=True)
        self.df_concat = self.df_concat.sort_values(by="Volume", ascending=True).reset_index(drop=True)
        self.function = self.titrations[0].function
        self.p0 = self.titrations[0].p0
        self.bounds = self.titrations[0].bounds
        self.wl = self.titrations[0].wl
        
    def fit(self, df):
        
        self.K, self.n, self.K_err, self.n_err = fit_mvh(self.function, df["r"], df["r_Cf"], self.p0, self.bounds)
                            
    def get_r2(self, df):
            
        self.r2 = calculate_R_squared(self.function, df["r"], df["r_Cf"], (self.K, self.n))
    
    def plot(self, df):
        
        plot_mvh(df, self.K, self.n, self.K_err, self.n_err, self.r2)
        
    def mvh_analysis_concat(self, df):
        
        print_initializing_mvh()
        print_wl_selected(self.wl)
        self.fit(df)
        self.get_r2(df)
        print_fit_r2_with_err(self.K, self.n, self.K_err, self.n_err, self.r2)
        self.plot(df)
        print_analysis_complete()