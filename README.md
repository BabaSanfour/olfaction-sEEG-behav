# Olfaction sEEG Behavior Project

This repository contains the breakdown of data and scripts for the Olfaction sEEG Behavior project. The files are organized into multiple directories.

## Directory Structure

```plaintext
olfaction_scripts_old/
│
├── 0_RawData_trc_avi_behav/                         # Raw Data
│   ├── S0X_iden_ordre=FM                            # Subject '0X' with identifier 'iden' and data presented in order 'FM'
│   │    ├── S0X_iden_Behavior                       # Subject '0X' Behavior Data (configPos, header, ParaOlfacto, raw, res, xlsx files)
│   │    └── S0X_iden_SignalsSEEG                    # Subject '0X' sEEG Data (TRC files)
│   └── ...                                          # 6 Subjects in total
│
├── 1bis_OE_BaseSam/                                 # Behavioral results, figures, and various scripts
│   ├── JPlailly201306_seeg_ALS                      # Old scripts
│   │    ├── behavior
│   │    │    ├── plots                              # Collection of plots (PNG/PDF)
│   │    │    ├── respiration_amplitude              # Numpy arrays (.npy files)
│   │    │    ├── Compute_Plot_respi.ipynb           # Plotting notebook
│   │    │    ├── analysis_name.py                   # Python scripts for different analyses
│   │    │    └── ...                                # Additional CSV and XLS files with results
│   │    ├── conn_analysis                           # Connectivity analysis scripts
│   │    │    └── analysis_name.py                   # Python scripts for different analyses
│   │    ├── Scripts_concat                          # Concatenation scripts
│   │    │    ├── Concatenate_sigs_Odor_Rest.ipynb   # Notebook for concatenating signals
│   │    │    ├── Export_correct_sig.ipynb           # Notebook for exporting signals
│   │    │    └── ...                                # MATLAB files
│   │    ├── seeg_analysis                           # sEEG analysis scripts
│   │    │    ├── bytrial_envelop_cycles             # Figures (PNG)
│   │    │    ├── diff_cycles_maps_0-40              # Figures (PNG)
│   │    │    ├── analysis_name.py                   # Python scripts for different analyses
│   │    │    └── ...                                # Additional PNG, PDF, and XLS files
│   │    └── database_JPlailly201306_seeg_ALS.sqlite # SQLite database (usage unclear)
│   ├── OLD                                          # Legacy scripts (JPlailly201306_seeg_david and JPlailly201306_seeg_win)
│   └── Results_Behavior.py                          # Script for generating behavior figures (PDF, PNG)
│        ├── Figures                                 # Collection of figures (Behavior Data)
│        └── ...                                     # XLS and ODS files: dataframes of results
```

## Overview

- **Raw Data:** Contains raw behavioral and sEEG data for each subject.
- **Behavioral Analysis Scripts:** This directory includes various scripts for analyzing behavioral results, generating figures, and processing sEEG data.