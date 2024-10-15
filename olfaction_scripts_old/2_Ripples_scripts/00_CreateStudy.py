"""This script create a study, which mean a set of repository
which can be convenient to manage path. (might doesn't work on windows...)
"""
from brainpipe.system import study

if __name__ == "__main__":
    # Define path and name of the study :
    #path = '/media/karim/Datas4To/1_Analyses_Intra_EM_Odor/'
    path = '/media/1_Analyses_Intra_EM_Odor/'
    name = 'Ripples'
    # Create the study :
    st = study(name)
    st.add(path)
