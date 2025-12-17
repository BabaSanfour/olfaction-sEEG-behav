# Olfaction sEEG Behavior Analysis (Legacy Scripts)

## Overview
This directory (`olfaction_scripts_old`) contains the legacy codebase for analyzing olfactory behavior and sEEG data. It relies heavily on a custom library named `brainpipe` (currently missing/external).

## Directory Structure
*   **`3_BrainpipeScripts/`**: **CORE PIPELINE**. Contains the main scripts for data ingestion, preprocessing, power analysis, and RSA.
    *   `1_Intra_scripts/`: The most active and updated folder containing notebooks for the analysis pipeline.
    *   `0_scripts_Etienne/`: Older scripts and utilities (e.g., `01_CleanData.py`).
*   **`2_Ripples_scripts/`**: Specialized scripts for ripple detection and artifact rejection.
*   **`1bis_OE_BaseSame/`**: Older analysis folders.

## Critical Notes
*   **Missing Dependency**: The `brainpipe` library is required for most scripts.
*   **Configuration**: Please check `config.py` (if available) or look for hardcoded paths in `00_CreateStudy.py`.

## Quick Start
1.  **Environment**: Install dependencies from `requirements.txt`.
    ```bash
    pip install -r requirements.txt
    ```
2.  **Configuration**: Set up your data paths in `config.py` (or environment variables).
3.  **Data Ingestion**: Use `3_BrainpipeScripts/0_scripts_Etienne/01_CleanData.py` to convert raw Matlab files.
4.  **Analysis**: Refer to the README in `3_BrainpipeScripts/1_Intra_scripts/` for the notebook workflow.
