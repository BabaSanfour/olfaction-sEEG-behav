from seegpy.io import read_3dslicer_fiducial
from seegpy.pipeline import pipeline_labelling_ss
import numpy as np

"""
First, define the path to :
* Freesurfer root folder where all the subjects are stored
* If you performed your segmentation with MarsAtlas, you can also provide
  the root path to the BrainVisa folder
* The name of your subject
"""
main_path = '/media/karim/Datas4To/1_Analyses_Intra_EM_Odor/Implantation_Patients_MNI/'
fs_root = main_path + 'db_freesurfer/'
bv_root = fs_root
subj = 'S01_VACJ'
# define the path that contains electrode coordinates
xyz_path = main_path + 'db_slicer_labels/'+subj+'/recon.fcsv'
# define where to save the contact's labels
save_to = main_path + 'db_slicer_labels/'+subj+'/'

# read the coordinates and the contact names
df = read_3dslicer_fiducial(xyz_path)
c_xyz = np.array(df[['x', 'y', 'z']])
c_names = np.array(df['label'])

"""
Run the pipeline for labelling contacts :
* Here, we're going to search for voxels / vertices that are contained in a
  sphere of 5mm centered around each sEEG sites (`radius=5.`)
* Use `bipolar=True` to localize bipolar derivations
* Use `bad_label='none'` for each contact that can't be labelled
"""
pipeline_labelling_ss(save_to, fs_root, bv_root, subj, c_xyz, c_names,
                      bipolar=True, radius=5., bad_label='none', verbose=False)