import numpy as np
from os.path import join, exists
from os import makedirs

from brainpipe.system import study
from utils_functions import conv2

def detect_ripples(xamp_rip, xamp_lf, raw_sig, rip_th=2, rip_width=[10,100], rip_max=3,art_th=4):
    """
    Detect ripples and artefacts to select only physiological ripples
    Adapted from methods in Vaz et al. 2019,2020 and Norman et al. 2019

    Parameters
    ----------
    xamp_rip : array (npts, trials) - Hilbert amplitude [80-120 Hz]
    xamp_lf : array (npts, trials) - Hilbert amplitude [15-60 Hz]
    raw_sig : array (npts, trials) - bipolarized signal
    rip_th : threshold for detection > 2 SD
    rip_width : in sample (between 20 and 200 ms ripple duration)
    rip_max : threshold to be reached > 3 SD
    art_th : threshold for artefact > 5 SD

    Return
    ------
    rip_logic : mask array (npts, trials) - physio ripples
    art_logic : mask array (npts, trials) - artefacts (lf + IED)
    """

    npts, ntrials = raw_sig.shape

    # z-score amplitudes
    z_amp_rip = (xamp_rip - np.mean(xamp_rip)) / np.std(xamp_rip)
    z_amp_lf = (xamp_lf - np.mean(xamp_lf)) / np.std(xamp_lf)

    #detect potential ripples and artefacts
    rip_log = z_amp_rip > rip_th
    lf_log = z_amp_lf > art_th #(Yitzhak Norman, Science 2019)

    #compute eeg gradient and detect IED (Alex Vaz, Science 2019,2020)
    raw_diff = np.diff(raw_sig,1,0) #first derivative over time dimension
    raw_diff = np.insert(raw_diff,obj=-1,values=1,axis=0) #make raw and raw diff same size on time axis
    z_raw_diff = (raw_diff - np.mean(raw_diff)) / np.std(raw_diff)
    ied_log = np.abs(z_raw_diff) > art_th

    #combine art and IED
    art_log = np.logical_or(lf_log,ied_log)
    #expand windows around art by 50ms = 25samples (*2sides) #(Yitzhak Norman, Science 2019)
    art_log = conv2(np.double(art_log), np.ones((50,1)),mode='same')
    art_log = art_log > 0

    #remove ripples coinciding with artefacts
    rip_no_art = np.zeros((npts, ntrials))
    for t in range(ntrials):
        rip_no_art[:,t] = np.array([1 if (r==1 and a==0) else 0 \
                                    for r,a in zip(rip_log[:,t],art_log[:,t])])

    #Loop through trials to check parameters of detected ripples
    for t in range(ntrials):

        # loop through ripples remove ripples too short/too long/two small (3SD)
        if np.sum(rip_no_art[:,t]) > 0:
            indx_rip = np.where(rip_no_art[:,t] == 1)[0]
            rip_split = np.split(indx_rip, np.where(np.diff(indx_rip) != 1)[0]+1)

            for ripple in rip_split:
                if not rip_width[0] < len(ripple) < rip_width[1]:
                    rip_no_art[ripple,t] = 0
                if np.max(z_amp_rip[ripple,t]) < rip_max:
                    rip_no_art[ripple,t] = 0

        #once validated, merge ripples less than 16ms (8 samples) separated
        if np.sum(rip_no_art[:,t]) > 0:
            indx_rip = np.where(rip_no_art[:,t] == 1)[0]
            rip_split = np.split(indx_rip, np.where(np.diff(indx_rip) != 1)[0]+1)

            if len(rip_split) > 1:
                for i,ripple in enumerate(rip_split[:-1]):
                    if rip_split[i+1][0] - ripple[-1] < 8:
                        rip_no_art[ripple[-1]:rip_split[i+1][0],t] = 1

    return rip_no_art, art_log

if __name__ == "__main__":

    st = study('Ripples')
    reps = ['Retrieval_new_odors/','Retrieval_new_rec/','Encoding/']
    feat_save = ['f', 'time', 'xampl', 'x', 'fname', 'new_labels', 'channels', 'labels', 'xyz']

    for rep in reps:
        files = st.search('INFO_ampl.npz', folder = ('feature/'+rep))
        rep_save = join(st.path, 'ripples/'+rep)
        if not exists(rep_save):
            makedirs(rep_save)

        for fi in files:
            print('Processing --> ', rep, fi)
            loadname = join(st.path, 'feature/'+rep+fi)
            mat = np.load(loadname, allow_pickle = True) #['ripple' 'low_freq' 'HFA']
            xampl_rip, xampl_lf = mat['xampl'][0,...], mat['xampl'][1,...]
            raw_sig = mat['x']
            nelecs, npts, ntrials = xampl_rip.shape

            #compute and detect ripples on all electrodes
            rip_tot = np.zeros((nelecs, npts, ntrials))
            art_tot = np.zeros((nelecs, npts, ntrials))
            for el in range(nelecs):
                rip_tot[el], art_tot[el] = detect_ripples(xampl_rip[el],
                                                    xampl_lf[el],raw_sig[el])
            print('Total # of ripples ', np.sum(rip_tot)/10)

            #save data in ripples repertory
            dico_rip = {}
            for feat in feat_save:
                dico_rip[feat] = mat[feat]
            dico_rip['ripples'], dico_rip['artefacts'] = rip_tot, art_tot
            savename = loadname.replace('bipo_INFO_ampl.npz',
                                'rip_art.npz').replace('feature', 'ripples')
            np.savez(savename, **dico_rip)
