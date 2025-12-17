
import os
from pathlib import Path

# =============================================================================
# GLOBAL CONFIGURATION
# =============================================================================

# Define the root of the project relative to this config file
# Assuming config.py is at olfaction_scripts_old/
ROOT_DIR = Path(__file__).resolve().parent

# =============================================================================
# PATH SETTINGS
# =============================================================================

# Path to the RAW data (Matlab files, etc.)
# Default: Try to find '0_RawData_trc_avi_behav' in the project, or fallback to env var
# User should export OLFACTION_RAW_DATA or edit this line.
DEFAULT_RAW_PATH = os.environ.get('OLFACTION_RAW_DATA', str(ROOT_DIR / '0_RawData_trc_avi_behav'))

# Path where the 'brainpipe' study will be initialized
# This is usually where the 'database', 'feature', etc. folders live.
# Default: The '3_BrainpipeScripts' folder.
STUDY_PATH = os.environ.get('OLFACTION_STUDY_PATH', str(ROOT_DIR / '3_BrainpipeScripts'))

# =============================================================================
# STUDY SETTINGS
# =============================================================================

STUDY_NAME = 'Olfacto'

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def get_study_path():
    """Returns the configured study path as a string."""
    return STUDY_PATH

def get_raw_data_path():
    """Returns the raw data path as a string."""
    return DEFAULT_RAW_PATH
