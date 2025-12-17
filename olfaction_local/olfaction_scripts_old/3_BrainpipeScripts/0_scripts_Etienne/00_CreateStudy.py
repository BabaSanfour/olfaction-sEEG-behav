"""This script create a study, which mean a set of repository
which can be convenient to manage path. (might doesn't work on windows...)
"""
from brainpipe.system import study

if __name__ == "__main__":
    # Define path and name of the study :
    path = '/media/karim/Datas4To/Dropbox/Intra_EM/3_BrainpipeScripts/'
    name = 'Olfacto'
    # Create the study :
    st = study(name)
    st.add(path)
