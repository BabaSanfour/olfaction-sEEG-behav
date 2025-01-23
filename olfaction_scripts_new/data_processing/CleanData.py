# -*- coding: utf-8 -*-
"""
Data Cleaning and Bipolirization 

Description:
This script loads data cleans process, and bipolrize it. The processed data is saved in 
both Numpy compressed format (.npz) and MNE-compatible FIF format for further analysis.


Author: AL
Modification: HA
Last Modified: 2024-22-03
"""

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

import re
from collections import defaultdict
import mne
import numpy as np
import pandas as pd
from utils.config import sessions, subjects, raw_data_folder, derivatives_data_folder, sf, rej_elec


def clean_data(subject):
    """
    Cleans and processes raw data for a given subject, saving the results in npz and fif formats.

    Args:
        subject (str): The identifier for the subject whose data is to be cleaned.

    Returns:
        None: Processed data is saved to a specified directory.
    """
    subject_files_path = os.path.join(raw_data_folder, subject)
    
    names_df = pd.DataFrame(np.loadtxt(os.path.join(subject_files_path, "channel_names.txt"), dtype=str))
    coord_df = pd.DataFrame(np.loadtxt(os.path.join(subject_files_path, "channel_coords.txt"), dtype=str))
    descr_df = pd.DataFrame(np.loadtxt(os.path.join(subject_files_path, "channel_descriptions.txt"), dtype=str))
    
    all_info = pd.concat([names_df, coord_df, descr_df], axis=1)    
    all_info.to_csv(os.path.join(subject_files_path, f"{subject}_all_info.txt"), header=False, index=False)
    del names_df, coord_df, descr_df

    for session in sessions:
        print(f'\n{"=" * 50}\nProcessing Session: {session}')
        data_fname = os.path.join(subject_files_path, f"{subject}_{session}_sigs.npy")
        
        try:
            data = np.load(data_fname)
        except FileNotFoundError:
            print(f"[ERROR] File not found: {data_fname}. Skipping...\n")
            continue

        coord = (
            np.array(all_info.iloc[:, 0]),
            np.array(all_info.iloc[:, 1]),
            np.array(all_info.iloc[:, 2]),
            np.array(all_info.iloc[:, 3]),
            np.array(all_info.iloc[:, 4]),
        )
        chanmonop = np.array(coord, dtype=object)

        chan_xy = chanmonop[1:4]
        chan_names = chanmonop[0]
        chan_label = chanmonop[4]

        xyz = np.zeros((len(chan_names), 3), dtype=float)
        for k in range(len(chan_names)):
            xyz[k, 0] = float(chan_xy[0][k])
            xyz[k, 1] = float(chan_xy[1][k])
            xyz[k, 2] = float(chan_xy[2][k])
        xyz = np.round(xyz, decimals=1)  # Round to one decimal place

        cleaned_data = {
            'data': data[:len(chan_names), ...],  # Retain relevant channels
            'xyz': xyz,
            'channel': chan_names,
            'label': chan_label,
            'sf': sf
        }

        save_path = os.path.join(derivatives_data_folder, f"{subject}", f"{subject}_{session}_sigs.npz")
        np.savez(save_path, **cleaned_data)
        print(f"  Data saved to: {save_path}")

        # Save using MNE
        print("  Preparing data for MNE...")
        info = mne.create_info(ch_names=list(chan_names), sfreq=sf, ch_types="seeg")  # Set channel type to sEEG
        raw = mne.io.RawArray(data[:len(chan_names)], info)

        # Add channel positions if available
        montage = mne.channels.make_dig_montage(ch_pos=dict(zip(chan_names, xyz)), coord_frame="head")
        raw.set_montage(montage)

        # Save as FIF
        mne_save_path = os.path.join(derivatives_data_folder, f"{subject}", f"{subject}_{session}_raw.fif")
        raw.save(mne_save_path, overwrite=True)
        print(f"  MNE Raw data saved to: {mne_save_path}")


def find_bipolar_pairs(channels, xyz):
    """
    Identify bipolar pairs from a list of channel names and compute their coordinates.

    Args:
        channels (list or array): List of channel names (e.g., ["a'1", "a'2", "b1"]).
        xyz (dict): Dictionary with channel names as keys and their (x, y, z) coordinates as values.

    Returns:
        pairs (list): List of bipolar pair names in the format "prefixX-prefixY".
        anodes (list): List of anode channel names (e.g., ["a'2", "b2"]).
        cathodes (list): List of cathode channel names (e.g., ["a'1", "b1"]).
        xyz_pairs (np.ndarray): Array of shape (n_pairs, 3) containing the coordinates of the new bipolar channels.
    """
    prefix_groups = defaultdict(list)
    
    # Group channels by prefix and extract numeric parts
    for idx, ch in enumerate(channels):
        match = re.match(r"([a-zA-Z']+)(\d+)", ch)
        if match:
            prefix, num = match.groups()
            prefix_groups[prefix].append((int(num), ch))  # Store numeric suffix and channel name

    pairs = []
    anodes = []
    cathodes = []
    xyz_pairs = []

    # Create bipolar pairs and compute their coordinates
    for prefix, num_ch_pairs in prefix_groups.items():
        num_ch_pairs = sorted(num_ch_pairs, key=lambda x: x[0])  # Sort by numeric suffix
        for i in range(len(num_ch_pairs) - 1):
            anode_name = num_ch_pairs[i + 1][1]
            cathode_name = num_ch_pairs[i][1]
            pair_name = f"{anode_name}-{cathode_name}"
            pairs.append(pair_name)
            anodes.append(anode_name)
            cathodes.append(cathode_name)
            
            # Compute the midpoint of the coordinates
            midpoint = (np.array(xyz[anode_name]) + np.array(xyz[cathode_name])) / 2
            xyz_pairs.append(midpoint)
    
    return pairs, anodes, cathodes, np.round(np.array(xyz_pairs), decimals=1)

def bipolarization(raw):
    """
    Apply bipolar referencing to raw MNE data.

    Args:
        raw (mne.io.Raw): Raw MNE object containing data and metadata.

    Returns:
        raw_bip_ref (mne.io.Raw): Bipolar-referenced MNE object with updated channel names and coordinates.
    """
    channels = raw.ch_names
    xyz = raw.get_montage().get_positions()["ch_pos"]  # Get channel positions
    bipolar_pairs, anodes, cathodes, xyz_bi = find_bipolar_pairs(channels, xyz)

    print(f"Found {len(bipolar_pairs)} bipolar pairs.")
    
    # Apply bipolar referencing
    raw_bip_ref = mne.set_bipolar_reference(
        raw, anode=anodes, cathode=cathodes, ch_name=bipolar_pairs, drop_refs=True
    )

    # Create and assign a new montage with updated coordinates
    montage = mne.channels.make_dig_montage(ch_pos=dict(zip(bipolar_pairs, xyz_bi)), coord_frame="head")
    raw_bip_ref.set_montage(montage)
    
    return raw_bip_ref

def bipolarize_data(subject):
    """
    Load, bipolarize, and save MNE data for a given subject.

    Args:
        subject (str): Identifier of the subject to process.

    Returns:
        None: Saves the processed data in FIF format.
    """
    subject_files_path = os.path.join(derivatives_data_folder, subject)
    
    for session in sessions:
        print(f'\n{"=" * 50}\nProcessing Session: {session}')
        raw_fif_fname = os.path.join(subject_files_path, f"{subject}_{session}_raw.fif")
        
        try:
            raw = mne.io.read_raw_fif(raw_fif_fname, preload=True)
        except FileNotFoundError:
            print(f"[ERROR] File not found: {raw_fif_fname}. Skipping...\n")
            continue

        raw_bip_ref = bipolarization(raw)
        
        mne_save_path = os.path.join(subject_files_path, f"{subject}_{session}_bipo_raw.fif")
        raw_bip_ref.save(mne_save_path, overwrite=True)
        print(f"  Bipolar-referenced data saved to: {mne_save_path}")


def drop_electrodes(subject):
    """
    Drops bad and white matter electrodes (channels) from MNE `.fif` files for a given subject.

    Args:
        subject (str): Identifier of the subject to process.

    Returns:
        None: The cleaned data is saved to the derivatives folder.
    """
    subject_files_path = os.path.join(derivatives_data_folder, subject)
    

    for session in sessions:
        print(f'\n{"=" * 50}\nProcessing Session: {session}')
        
        npz_fpath = os.path.join(subject_files_path, f"{subject}_{session}_sigs.npz")
        raw_fif_fname = os.path.join(subject_files_path, f"{subject}_{session}_bipo_raw.fif")

        try:
            mat = np.load(npz_fpath, allow_pickle=True)
            labels = mat['label']
        except FileNotFoundError:
            print(f"[ERROR] .npz file not found: {npz_fpath}. Skipping...\n")
            continue
        bad_indices = rej_elec.get(subject, [])  # Get the list of bad indices for the current subject
        # Add WM electrodes to reject
        for elec_idx, label in enumerate(labels):
            if "WM" in label or "LCR" in label or "skull" in label or "out" in label:
                bad_indices.append(elec_idx)
        print(f"Identified {len(bad_indices)} electrodes to reject: {bad_indices}")

        try:
            raw = mne.io.read_raw_fif(raw_fif_fname, preload=True)
        except FileNotFoundError:
            print(f"[ERROR] File not found: {raw_fif_fname}. Skipping...\n")
            continue

        bad_channels = [raw.ch_names[idx] for idx in bad_indices if idx < len(raw.ch_names)]
        
        raw.info['bads'] = bad_channels
        print(f"Marked bad channels for {subject}: {bad_channels}")

        raw_clean = raw.copy().drop_channels(bad_channels)
        print(f"Dropped {len(bad_channels)} channels. Remaining channels: {len(raw_clean.ch_names)}")

        new_file_path = os.path.join(subject_files_path, f"{subject}_{session}_bipo_clean_raw.fif")
        raw_clean.save(new_file_path, overwrite=True)
        print(f"Saved cleaned data to: {new_file_path}")

def filtering(subject):
    """
    Applies notch filtering to remove powerline noise and its harmonics for a given subject.

    Args:
        subject (str): Identifier of the subject whose data is to be filtered.

    Returns:
        None: The filtered data is saved to the derivatives folder.
    """
    subject_files_path = os.path.join(derivatives_data_folder, subject)

    for session in sessions:
        print(f'\n{"=" * 50}\nProcessing Session: {session}')
        raw_fif_fname = os.path.join(subject_files_path, f"{subject}_{session}_bipo_clean_raw.fif")

        try:
            raw = mne.io.read_raw_fif(raw_fif_fname, preload=True)
        except FileNotFoundError:
            print(f"[ERROR] File not found: {raw_fif_fname}. Skipping...\n")
            continue

        print(f"Applying notch filter to {raw_fif_fname}...")
        raw.notch_filter(freqs=np.arange(50, 251, 50), fir_design='firwin')

        mne_save_path = os.path.join(subject_files_path, f"{subject}_{session}_bipo_clean_filt_raw.fif")
        raw.save(mne_save_path, overwrite=True)
        print(f"  Filtered data saved to: {mne_save_path}")

def main():
    """
    Main pipeline for cleaning and bipolarizing electrophysiological data.
    """
    print(f"\n{'#' * 50}\nStarting data processing...\n{'#' * 50}")
    
    for subject in subjects:
        print(f"Processing Subject: {subject}")
        print("Cleaning Data")
        clean_data(subject)
        print("Bipolarizing Data")
        bipolarize_data(subject)
        print("Dropping Bad Electrodes")
        drop_electrodes(subject)
        print("Filtering")
        filtering(subject)

    print(f"\n{'#' * 50}\nData processing completed!\n{'#' * 50}")

if __name__ == '__main__':
    main()
