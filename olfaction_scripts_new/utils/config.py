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
    'CHAF' : ['1','2','3','4','5','7','8','9'],
    'VACJ' : ['14','15','16','17','10','11','12','13'], 
    'SEMC' : ['5','7','8','9','10','11','12','13'],
    'LEFC' : ['1','2','3','4','14','15','16','17'],
    'FERJ' : ['12','13','16','17','1','2','5','7'],
    'PIRJ' : ['1','4','9','18','6','7','15','5'],
}
phases = ['odor']
