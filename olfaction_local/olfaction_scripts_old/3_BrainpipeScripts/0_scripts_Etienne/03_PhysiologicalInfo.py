"""In this script, we are going to get physiological informations
from your seeg data coordinates (brodmann area/matter/gyrus...),
save everything in an excel file. Everything is going to be saved
in the physiology folder of your study.
"""
import numpy as np
import pandas as pd
from pandas import ExcelWriter
from os import path

from brainpipe.system import study
from brainpipe.preprocessing import xyz2phy
if __name__ == "__main__":
    # Get info from monop / bipo files?
    bipo = True
    if bipo:
        pattern = 'allfilter1_bipo.npz'
        savexls = 'OlfactoPhysio_bipo.xlsx'
    else:
        pattern = 'allfilter1.npz'
        savexls = 'OlfactoPhysio_monop.xlsx'
    # Define path and data files :
    st = study('Olfacto')
    files = st.search(pattern, folder='database/TS_E_all_cond_by_block_trigs_th40_art400_30_250/') # Search one file per subject

    # Create an empty "table":
    writer = ExcelWriter(path.join(st.path, 'physiology', savexls))

    for fi in files:
        suj = fi.split('_E1E2')[0]
        # Load your file :
        file2load = path.join(st.path, 'database/TS_E_all_cond_by_block_trigs_th40_art400_30_250/', fi)
        print('-> File: ', file2load)
        mat = np.load(file2load) # Load file
        # Get physiological informations :
        xyzObj = xyz2phy()
        df = xyzObj.get(mat['xyz']) 
        df['channel'] = mat['channel']
        # Add to a new sheet :
        df.to_excel(writer, sheet_name=suj)
    writer.save()
    
    

