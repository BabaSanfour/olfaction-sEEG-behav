import os.path as op

import numpy as np
import pandas as pd

from seegpy.io import read_3dslicer_fiducial
from seegpy.pipeline import pipeline_labelling_ss
from seegpy.labelling.lab_vol import labelling_contacts_vol_fs_mgz


fs_root = '/media/karim/Datas4To/1_Analyses_Intra_EM_Odor/Implantation_Patients_MNI/db_freesurfer'
suj = 'S06_MICP'
bv_root = None
# define the path that contains electrode coordinates
xyz_path = '/media/karim/Datas4To/1_Analyses_Intra_EM_Odor/Implantation_Patients_MNI/db_slicer_labels/S06_MICP/recon.fcsv'
# define where to save the contact's labels
save_to = '/media/karim/Datas4To/1_Analyses_Intra_EM_Odor/Implantation_Patients_MNI/db_slicer_labels/S06_MICP'

# read the coordinates and the contact names
df = read_3dslicer_fiducial(xyz_path)
c_xyz = np.array(df[['x', 'y', 'z']])
c_names = np.array(df['label'])


"""
bipolar = True : info anat sur les contacts en bipolaires (coordonnée définie
comme le milieu de deux contacts successifs genre B2-B1)
bipolar = False : info anat sur les contacts en monopolaires seuls
"""
bipolar = True
"""
radius = c'est le rayon de la boule qui est utilisée pour trouver le voxel
labélisé le plus proche. Si il trouve rien, il va te mettre un `bad_label`
"""
radius = 5.
bad_label = 'none'
pipeline_labelling_ss(save_to, fs_root, bv_root, suj, c_xyz, c_names,
                      bipolar=bipolar, radius=radius, bad_label=bad_label,
                      verbose=False, testing=False)

"""
reload the exported anatomical informations for both monopolar and bipolar
contacts
"""
anat_path = op.join(save_to, f"{suj}_radius-{radius}.xlsx")
df = {}
for derivation in ['monopolar', 'bipolar']:
    df[derivation] = pd.read_excel(anat_path, sheet_name=derivation)

"""
Here we put all high resolutions files for the hippocampe. Then the rest of the
script localize each contact using the selected file and add it to the original
anatomical file. I think you should only select files that end with
"FSvoxelSpace"
"""
hires_hipp_files = [
    'rh.hippoAmygLabels-T1.v21.CA.FSvoxelSpace',
    'lh.hippoAmygLabels-T1.v21.CA.FSvoxelSpace',
    'rh.hippoAmygLabels-T1.v21.FS60.FSvoxelSpace',
    'lh.hippoAmygLabels-T1.v21.FS60.FSvoxelSpace',
    'rh.hippoAmygLabels-T1.v21.HBT.FSvoxelSpace',
    'lh.hippoAmygLabels-T1.v21.HBT.FSvoxelSpace',
    'rh.hippoAmygLabels-T1.v21.FSvoxelSpace',
    'lh.hippoAmygLabels-T1.v21.FSvoxelSpace',
]

for derivation in ['monopolar', 'bipolar']:
    # selection the derivation type and get the (xyz) coordinates
    df_der = df[derivation]
    xyz_der = np.array(df_der[['x_scanner', 'y_scanner', 'z_scanner']])

    # now, localize the sites using several files
    for file in hires_hipp_files:
        lab = labelling_contacts_vol_fs_mgz(
            fs_root, suj, xyz_der, radius=radius, file=file,
            bad_label=bad_label, verbose=None)
        df_der[file] = lab.squeeze()
    df[derivation] = df_der

"""
re-export the file
"""
save_as = op.join(save_to, f"{suj}_radius-{radius}_hip-hires.xlsx")
with pd.ExcelWriter(save_as) as writer:
    for sheet, wdf in df.items():
        wdf.to_excel(writer, sheet_name=sheet)
