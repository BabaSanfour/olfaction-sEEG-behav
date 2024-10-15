# Olfaction sEEG Behavior Project

This repository contains the breakdown of data and scripts for the Olfaction sEEG Behavior project. The files are organized into multiple directories.

## Directory Structure

```plaintext
olfaction_scripts_old/
│
├── 0_RawData_trc_avi_behav/                         # Raw Data
│   ├── S0X_iden_ordre=FM/                           # Subject '0X' with identifier 'iden' and data presented in order 'FM'
│   │    ├── S0X_iden_Behavior/                      # Subject '0X' Behavior Data (configPos, header, ParaOlfacto, raw, res, xlsx files)
│   │    └── S0X_iden_SignalsSEEG/                   # Subject '0X' sEEG Data (TRC files)
│   └── ...                                          # 6 subjects in total
│
├── 1bis_OE_BaseSam/                                 # Behavioral results, figures, and various scripts
│   ├── JPlailly201306_seeg_ALS/                     # Old scripts
│   │    ├── behavior/
│   │    │    ├── plots/                             # Collection of plots (PNG/PDF)
│   │    │    ├── respiration_amplitude/             # Numpy arrays (.npy files)
│   │    │    ├── Compute_Plot_respi.ipynb           # Plotting notebook
│   │    │    ├── analysis_name.py                   # Python scripts for various analyses
│   │    │    └── ...                                # Additional CSV and XLS files with results
│   │    ├── conn_analysis/                          # Connectivity analysis scripts
│   │    │    └── analysis_name.py                   # Python scripts for various analyses
│   │    ├── Scripts_concat/                         # Concatenation scripts
│   │    │    ├── Concatenate_sigs_Odor_Rest.ipynb   # Notebook for concatenating signals
│   │    │    ├── Export_correct_sig.ipynb           # Notebook for exporting signals
│   │    │    └── ...                                # MATLAB files
│   │    ├── seeg_analysis/                          # sEEG analysis scripts
│   │    │    ├── bytrial_envelop_cycles/            # Figures (PNG)
│   │    │    ├── diff_cycles_maps_0-40/             # Figures (PNG)
│   │    │    ├── analysis_name.py                   # Python scripts for various analyses
│   │    │    └── ...                                # Additional PNG, PDF, and XLS files
│   │    └── database_JPlailly201306_seeg_ALS.sqlite # SQLite database (usage unclear)
│   ├── OLD/                                         # Legacy scripts (JPlailly201306_seeg_david and JPlailly201306_seeg_win)
│   └── Results_Behavior                             
│        ├── Figures                                 # Collection of figures (Behavior Data)
│        └── ...                                     # XLS and ODS files: dataframes of results
├── 2_Ripples_scripts/                               # Collection of Python scripts
│   ├── Classif_power_ML/                            # Old scripts
│   │    └── ... name.ipynb
│   ├── Vaz_ripples/                                 # Poster.pdf, detectripples.py, and log and MATLAB files
│   ├── 00_ripples_distrib.ipynb                     
│   └── ... analysis_name.py                         # Collection of Python scripts (processing and utilities)
├── 3_BrainpipeScripts/                              # Collection of Python scripts and Brainpipe source code
│   ├── 0_scripts_Etienne/                           
│   │    └── ... name_analysis.py                    # Scripts for cleaning data and computing power
│   ├── 1_Intra_scripts/                             #
│   │    ├── behaviour/
│   │    │    ├── ... analysis_name.py               # Python scripts for various analyses
│   │    │    └── ... name.ipynb                     # Notebooks for various analyses
│   │    ├── compute_canada/
│   │    │    ├── ... analysis_name.py               # Python scripts for various analyses
│   │    │    └── ... name.sh                        # Batch scripts
│   │    ├── data/                                   # Numpy files
│   │    ├── mne_expls/                              # Plotting notebook
│   │    ├── R_stats/                                # R scripts
│   │    ├── ... name.py
│   │    └── ... name.ipynb
│   ├── Notebook/                             
│   │    ├── Old_stuff/
│   │    │    └── ... name.ipynb                    # Notebooks for various analyses
│   │    ├── Other scripts/
│   │    │    └── ... name.ipynb                    # Notebooks for various analyses
│   │    ├── utils/                                  
│   │    │    └── ... name.ipynb                    # Notebooks for various analyses
│   │    └── ... name.ipynb                          # Notebooks for classification and ERP analysis
│   ├── ripples_hfo/                             
│   │    └── ... name.ipynb                          # Notebooks for HFO
│   └── src/                                         # VisPy source code
├── 4_Visbraun_scripts/                                                           
│   ├── OLD/                             
│   │    └── ... name.py                             # Plotting scripts
│   ├── src/                                         # VisPy source code
│   └── ... name.ipynb                               # Plotting scripts
├── Data_NicoF_CRNL/                                                           
│   ├── Encoding_Odor/                               # Clean data (npz files)                             
│   ├── Rappel_odor/                                 # Clean data (npz files) 
│   └── Recognition/                                 # Clean data (npz files)
├── Implantation_Patients_MNI/                       # Subjects' implantation data                                    
├── OE_Matrices_NoArt/
│   ├── _DB_2021_Janv/
│   │    ├── TS_E_all_cond_by_block_trigs_odors/
│   │    │    ├── _group_files/                      
│   │    │    │    └── ... {Ident}_odor_{i}.npy      # Unclear usage
│   │    │    ├── {ident}/                          # Subject data (available for all subjects)
│   │    │    │    ├── R1/
│   │    │    │    │    └── ... name.npy            # correct_ts_{i}_by_odor+trigs.npy files
│   │    │    │    ├── R2/                          # Same as above
│   │    │    │    ├── R3/                          # Same as above
│   │    │    │    └── ... name.txt                 # Channel descriptions files
│   │    ├── TS_R_all_cond_by_block_trigs_odors/     # Same as above
│   │    └── TS_R_all_cond_by_block_trigs_rec_time/  # Same as above
│   ├── TS_E_all_cond_by_block_trigs_odors/          # Same as above
│   │    ├── _group_files/
│   │    │    └── ... {Ident}_odor_{i}.npy           # Unclear usage
│   │    ├── _mat_by_odor/
│   │    │    └── ... {Ident}_odor_{i}.mat           # Unclear usage
│   │    ├── {ident}/                                 # Subject data (available for all subjects)
│   │    │    ├── R1/
│   │    │    │    └── ... name.npy                  # correct_ts_{i}_by_odor+trigs.npy files
│   │    │    ├── R2/                                # Same as above
│   │    │    ├── R3/                                # Same as above
│   │    │    └── ... name.txt                       # Channel descriptions files
│   ├── TS_R_all_cond_by_block_trigs_odors/          # Same as above
│   └── TS_R_all_cond_by_block_trigs_rec/            # Same as above


```

## Overview

- **Raw Data:** Contains raw behavioral and sEEG data for each subject.
- **Behavioral Analysis Scripts:** This directory includes various scripts for analyzing behavioral results, generating figures, and processing sEEG data.