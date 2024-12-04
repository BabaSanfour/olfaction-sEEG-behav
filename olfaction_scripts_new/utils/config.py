import os
import sys
import numpy as np

# Set paths
user = os.path.expanduser('~')
data_folder = os.path.join(user, "scratch/data/iEEG_olfaction_data")
raw_data_folder = os.path.join(data_folder, "raw_data")
TS_E_all_cond_by_block_trigs_odors_folder = os.path.join(data_folder, "TS_E_all_cond_by_block_trigs_odors/")


# Set subjects and triggers
subjects = ['CHAF', 'VACJ', 'SEMC', 'LEFC', 'FERJ', 'PIRJ']
trigger_to_calculate = {
    'CHAF' :   ['01', '02', '03', '04', '15', '16', '17', '18'],
    'VACJ' :   ['01', '02', '03', '15', '16', '17', '18'],
    'SEMC' :   ['01', '02', '03', '16', '17', '18'],
    'LEFC' :   ['01', '02', '03', '04', '15', '16', '18'],
    'FERJ' :   ['01', '02', '03', '04', '16', '18'],
    'PIRJ' :   ['01', '02', '03', '04', '15', '16', '18'],
    }
phases = ['oders']
