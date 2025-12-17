# Intra-Cranial Analysis Scripts

This directory contains the core notebooks for the sEEG pipeline, including Power Analysis, RSA, and Classification.

## Pipeline Overview

### 1. Setup & Preprocessing
*   **`00_CreateStudy.ipynb`**: Required Step 0. Initializes the "Study" object and paths.
*   **`02_Bipolarized.py`**: Applies bipolar referencing to cleaned `.npz` files.
*   **`02_Remove_bad_electrodes.ipynb`**: Quality control and electrode rejection.

### 2. Feature Extraction
*   **`03_Compute_Power.ipynb`**: Major notebook for computing Time-Frequency power (Hilbert/Wavelet).
    *   *Input*: Bipolarized `.npz` files.
    *   *Output*: Feature `.npz` files (e.g., `_power.npz`).

### 3. Representational Similarity Analysis (RSA)
*   **`RSA.py`**: Core library file containing `within_rsa` and `btw_rsa` functions.
*   **`04_RSA_Plot_Matrices_Stats.ipynb`**: Visualizes Encoding vs Retrieval similarity matrices.
*   **`05_TPSim_RSA_by_cond_sklearn_in_time.ipynb`**: Advanced Temporal Pattern Similarity analysis.

### 4. Classification (Machine Learning)
*   **`05_Classif_Power_Unbalanced_...ipynb`**: Set of notebooks for classifying conditions (Low vs High performance) using Power features.
*   **`09_Confusion_matrices...ipynb`**: Visualization of classification results.

### 5. Visualization & Figures
*   **`07_figures_results_&_count_by_rois.ipynb`**: Generates summary figures for ROIs.
*   **`13_figures_results_&_count_by_rois.ipynb`**: Updated figure generation.

## Usage
Run notebooks in the order numbered above (01 -> 02 -> 03...). Ensure `brainpipe` is importable.
