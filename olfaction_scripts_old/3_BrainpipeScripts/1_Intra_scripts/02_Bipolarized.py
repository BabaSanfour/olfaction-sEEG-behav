"""Script to bipolarize data.
"""
import numpy as np
from os import path

from brainpipe.preprocessing import bipolarization
from brainpipe.system import study


#if __name__ == "__main__":
# Get all files include in the database (ignore files containing '_bipo'):
reps =  [   
            'Retrieval_new_odors',
            'Retrieval_new_rec',
        ]

st = study('Ripples')
for rep in reps:
    files = [k for k in st.search('_odor_', folder=('database/'+rep)) if not k.find('_bipo')+1]
    for fi in files:
        if fi.endswith('.npz'):
            # Load file :

            loadname = path.join(st.path, 'database/'+rep, fi)
            print('-> Bipolarize: '+loadname)
            mat = np.load(loadname, allow_pickle=True)
            print('all files', mat.files)

            #~ print(mat['channel'])
            # Bipolarize your data :
            x_bip, chan_bip, label_bip, xyz_bip = bipolarization(mat['x'], [mat['channel'][i][0] for i in range(len(mat['channel']))],
                                                                [mat['label'][i][0] for i in range(len(mat['label']))], xyz=mat['xyz'],)
            #~ print(mat['channel'])
            # Save bipolarized data :
            mat = dict(mat)
            mat['x'], mat['channel'], mat['xyz'], mat['label'] = x_bip, chan_bip, xyz_bip, label_bip
            savename = loadname.replace('.npz', '_bipo.npz')
            np.savez(savename, **mat)
