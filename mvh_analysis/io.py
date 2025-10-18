# -*- coding: utf-8 -*-
"""
Created on Fri Oct 10 17:42:33 2025

@author: JL
"""

from pathlib import Path
import pandas as pd

def load_raw_df(input_file):
    
    df = pd.read_csv(input_file)
    
    return df

def save_df_to_csv(df, filename, folder):
    
    folder.mkdir(parents=True, exist_ok=True)
    save_name = folder / filename
    df.to_csv(save_name, index=False)
    print(f"Dataframe saved to {save_name}.")


    
    