"""
Data Cleaning Script for MATLAB Format Conversion and Saving to Numpy Format

This script loads data from MATLAB `.mat` files, cleans and processes the data (including channel names, 
coordinates, and data matrices), and saves it to Numpy format (.npz).

Author: AL
Modification: HA
Date: 2024-12-03
"""

import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

import os
import numpy as np
from scipy.io import loadmat
from utils.config import TS_E_all_cond_by_block_trigs_odors_folder, subjects, trigger_to_calculate, phases

data_path = os.path.join(TS_E_all_cond_by_block_trigs_odors_folder, '_mat_by_odor')

def _clean_data(pathdata, suj, phase, trigger):
    """
    Cleans and processes the data from a MATLAB `.mat` file, extracting and organizing necessary 
    information such as channel names, coordinates, and data matrix, and saving it in a more convenient
    Numpy format.

    Args:
        pathdata (str): The directory where the MATLAB files are located.
        suj (str): Subject identifier.
        phase (str): The phase of the experiment (e.g., 'encoding', 'resting').
        trigger (str): The trigger condition to process.

    Returns:
        dict: A dictionary containing cleaned data with keys: 'x', 'xyz', 'channel', 'label', and 'sf'.
        str: The short filename (without the extension) for saving.
    """
    
    file = f'{suj}_{phase}_{trigger}.mat'
    print(f'\n{"="*50}\nProcessing file: {file}')
    
    fname = os.path.join(pathdata, file)
    try:
        mat = loadmat(fname)
    except FileNotFoundError:
        print(f"[ERROR] File not found: {file}. Skipping...\n")
        return None, None

    # Extract and clean channel names and coordinates (xyz)
    chan_xy = mat['chanmonop'][1:4]  # Coordinates for channels (x, y, z)
    chan_names = mat['chanmonop'][0]  # Channel names
    chan_label = mat['chanmonop'][4]  # Channel labels (if applicable)
    
    # Initialize the xyz coordinates matrix
    xyz = np.zeros([len(chan_names), 3], dtype=float)
    print("  Extracting channel coordinates (xyz):")
    for k in range(len(chan_names)):
        xyz[k, 0] = float(chan_xy[0, k][0])
        xyz[k, 1] = float(chan_xy[1, k][0])
        xyz[k, 2] = float(chan_xy[2, k][0])
    
    # Round coordinates to one decimal place
    xyz = np.round(xyz, decimals=1)
    print(f"    Channel Coordinates (xyz):\n    {xyz}")

    # Clean and reshape the data matrix 'x'
    x = mat['x']
    print(f"  Raw data shape: {x.shape}")
    
    if x.ndim == 2:  # In case x is a 2D matrix, reshape it to 3D
        print("    Reshaping 2D data to 3D...")
        x = x[..., np.newaxis]
    
    x = x.swapaxes(0, 1).swapaxes(1, 2)  # Reorder dimensions
    print(f"    Cleaned data shape: {x.shape}")

    # Create a dictionary with the cleaned data
    dico = dict(mat)
    dico['x'] = x[:len(chan_names), ...]  # Keep only relevant channels
    dico['xyz'] = xyz  # Cleaned channel coordinates
    dico['channel'] = chan_names  # Cleaned channel names
    dico['label'] = chan_label  # Channel labels
    dico['sf'] = dico['sf'][0][0]  # Sampling frequency (assumed to be a scalar)
    
    # Clean up the dictionary by removing unnecessary fields
    del dico['chanmonop']
    del dico['__header__']
    del dico['__globals__']
    del dico['__version__']

    # Return the cleaned data dictionary and the short filename
    return dico, fname.split('.mat')[0]


def main():
    """
    Main function that loops over subjects, phases, and triggers to process and save data.
    For each combination of subject, phase, and trigger, it calls the `_clean_data` function
    and saves the cleaned data in Numpy `.npz` format.
    """
    
    print(f"\n{'#'*50}\nStarting data processing...\n{'#'*50}")
    
    for sub in subjects:
        for phase in phases:
            for tr in trigger_to_calculate[sub]:
                print(f"\n{'-'*50}")
                print(f"Processing Subject: {sub} | Phase: {phase} | Trigger: {tr}")
                
                dico, shortname = _clean_data(data_path, sub, phase, tr)
                
                if dico is None:
                    continue
                
                path_to_save = os.path.join(data_path, 'Encoding_By_Odor', f'{shortname}.npz')
                
                np.savez(path_to_save, **dico)
                print(f"  Data saved to: {path_to_save}")
    
    print(f"\n{'#'*50}\nData processing completed!\n{'#'*50}")


if __name__ == '__main__':
    main()
