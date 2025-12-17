# Ripple Detection Scripts

This folder contains scripts specialized for detecting ripples (>80Hz oscillations) and performing artifact rejection.

## Key Files
*   **`05_detect_ripples.py`**: The main script for detection.
    *   *Method*: Likely based on Vaz et al. 2019 (Hilbert amplitude thresholding).
    *   *Input*: Raw sEEG data (`x` variable).
    *   *Output*: `rip_art.npz` files containing detected events.
*   **`utils_functions.py`**: Helper functions for convolution (`conv2`) and electrode manipulation.

## Usage
This pipeline appears independent of the main Power/RSA pipeline but shares the same validation logic (creating a Study object).
