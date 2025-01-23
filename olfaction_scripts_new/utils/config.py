import os
import sys
import numpy as np

# Set paths
user = os.path.expanduser('~')
data_folder = os.path.join(user, "scratch/data/iEEG_olfaction_data")
raw_data_folder = os.path.join(data_folder, "raw_data")
derivatives_data_folder = os.path.join(data_folder, "derivatives")

# Set subjects and triggers
subjects = ['CHAF', 'FERJ', 'LEFC', 'SEMC', 'PIRJ', 'VACJ', 'MICP'] #respectively S0 S1 S2 S3 S4 S5


rej_elec = {
    'CHAF' : [27,76,98,102,103,104],
    'VACJ' : [0,66,67,68,122,128], 
    'SEMC' : [36],
    'PIRJ' : [10,14,15,25,26,41,45,47,48,53,54,55,57,58,63,64,70,71,72,79,80,81,82,83,],
    'LEFC' : [47,48,49,50,140,141,],
    'MICP' : [4,9,13,18,19,20,71,],
    'FERJ': [],
            }


trigger_to_calculate = {
    'CHAF' : ['1','2','3','4','5','7','8','9'],
    'VACJ' : ['14','15','16','17','10','11','12','13'], 
    'SEMC' : ['5','7','8','9','10','11','12','13'],
    'LEFC' : ['1','2','3','4','14','15','16','17'],
    'FERJ' : ['12','13','16','17','1','2','5','7'],
    'PIRJ' : ['1','4','9','18','6','7','15','5'],
}
sessions = ["E1", "E2"]
phases = ['odor']


sf = 512. # Sampling frequency
