import numpy as np

from seegpy.io import read_3dslicer_fiducial
from seegpy.pipeline import pipeline_labelling_ss


fs_root = '/media/karim/Datas4To/1_Analyses_Intra_EM_Odor/Implantation_Patients_MNI/db_freesurfer'
suj = 'S02_SEMC'
bv_root = None
# define the path that contains electrode coordinates
xyz_path = '/media/karim/Datas4To/1_Analyses_Intra_EM_Odor/Implantation_Patients_MNI/db_slicer_labels/S02_SEMC/recon.fcsv'
# define where to save the contact's labels
save_to = '/media/karim/Datas4To/1_Analyses_Intra_EM_Odor/Implantation_Patients_MNI/db_slicer_labels/S02_SEMC'

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