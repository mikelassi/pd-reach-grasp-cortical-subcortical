# ===============================
# CONFIGURATION (EDIT THIS ONLY)
# ===============================
load_path = r"C:\Percorso\Ai\Dati"   # <-- Path ai pickle files EEG_DATA_SHORT, LFP_DATA_SHORT, CLEAN_DATA

n_iterations = 1000    # Numero di iterazioni per il bootstrap

save_path = r"C:\Percorso\A\Cartella\Di\Destinazione"   # Cambiare al path di salvataggio

import os
import pickle
import mne
import numpy as np
from mne_connectivity import spectral_connectivity_epochs
import warnings

with open(os.path.join(load_path, 'EEG_DATA_SHORT.pkl'), 'rb') as f:
    EEG_DATA_SHORT = pickle.load(f)
with open(os.path.join(load_path, 'LFP_DATA_SHORT.pkl'), 'rb') as f:
    LFP_DATA_SHORT = pickle.load(f)
with open(os.path.join(load_path, 'CLEAN_DATA.pkl'), 'rb') as f:
    CLEAN_DATA = pickle.load(f)

# ----------------------------------------------------------- #
# CONFIG
# ----------------------------------------------------------- #
subject_list = ['wue02', 'wue03', 'wue05', 'wue06', 'wue07', 'wue09', 'wue10', 'wue11']

GC_BOOTSTRAP_RESULTS = {
    subject_id: {
        'LFP2EEG': {},
        'EEG2LFP': {},
        'NET': {},
        'freqs':{},
        'NET_CLUST_MEAN': {},
        'L2E_CLUST_MEAN': {},
        'E2L_CLUST_MEAN': {}
    }
    for subject_id in subject_list
}


# ----------------------------------------------------------- #
# FUNZIONE BOOTSTRAP PER UNA ITERAZIONE
# ----------------------------------------------------------- #
def bootstrap_iteration(eeg_data, lfp_data, min_len, sfreq, conn_obj,
                        GC_lags, fmin, fmax):

    n_trials = eeg_data.shape[0]
    n_cluster = eeg_data.shape[1]

    lfp_shifted = np.empty_like(lfp_data)

    for tr in range(n_trials):

        # Random shift of n_samples >= GC_lags
        shift_len = np.random.randint(GC_lags, min_len - GC_lags)
        lfp_shifted[tr, 0, :] = np.roll(lfp_data[tr, 0, :], -shift_len)

    # Epochs combinati
    combined = np.concatenate([lfp_shifted, eeg_data[:, :, :]], axis=1)
    conn_obj._data = combined

    indices_L2E = (
        [np.array([0])] * n_cluster,
        [np.array([i]) for i in range(1, n_cluster + 1)]
    )

    indices_E2L = (
        [np.array([i]) for i in range(1, n_cluster + 1)],
        [np.array([0])] * n_cluster
    )

    # GC LFP→EEG
    con_L2E = spectral_connectivity_epochs(
        conn_obj, method='gc', indices=indices_L2E, sfreq=sfreq,
        fmin=fmin, fmax=fmax, mode='multitaper', mt_adaptive=True,
        mt_low_bias=True, verbose=False, gc_n_lags=GC_lags
    )
    con_L2E_tr = spectral_connectivity_epochs(
        conn_obj, method='gc_tr', indices=indices_L2E, sfreq=sfreq,
        fmin=fmin, fmax=fmax, mode='multitaper', mt_adaptive=True,
        mt_low_bias=True, verbose=False, gc_n_lags=GC_lags
    )
    gc_L2E = con_L2E.get_data().squeeze() - con_L2E_tr.get_data().squeeze()

    # GC EEG→LFP
    con_E2L = spectral_connectivity_epochs(
        conn_obj, method='gc', indices=indices_E2L, sfreq=sfreq,
        fmin=fmin, fmax=fmax, mode='multitaper', mt_adaptive=True,
        mt_low_bias=True, verbose=False, gc_n_lags=GC_lags
    )
    con_E2L_tr = spectral_connectivity_epochs(
        conn_obj, method='gc_tr', indices=indices_E2L, sfreq=sfreq,
        fmin=fmin, fmax=fmax, mode='multitaper', mt_adaptive=True,
        mt_low_bias=True, verbose=False, gc_n_lags=GC_lags
    )
    gc_E2L = con_E2L.get_data().squeeze() - con_E2L_tr.get_data().squeeze()

    gc_net = gc_E2L - gc_L2E

    return gc_L2E, gc_E2L, gc_net, con_L2E.freqs

chan_names = EEG_DATA_SHORT['wue02'][0]['chan_names']
SIGNIFICANT_CLUSTER_CHANNELS = ['C1', 'CCP1h', 'CCP2h', 'CCP3h', 'Cz', 'FCC1h', 'FCC2h']
SIGNIFICANT_CLUSTER_CHANNELS_LEFT_HAND = ['C2', 'CCP2h', 'CCP1h','CCP4h', 'Cz', 'FCC2h', 'FCC1h']
cluster_idx = [chan_names.index(ch) for ch in SIGNIFICANT_CLUSTER_CHANNELS]
cluster_idx_left_hand = [chan_names.index(ch) for ch in SIGNIFICANT_CLUSTER_CHANNELS_LEFT_HAND]

warnings.filterwarnings("ignore", category=RuntimeWarning)

# ----------------------------------------------------------- #
# LOOP PER SOGGETTO
# ----------------------------------------------------------- #
for subject_id in subject_list:

    print(f"\n=== Bootstrap GC for subject: {subject_id} ===")

    if subject_id in ['wue05', 'wue09']:
        cluster_idx_subj = cluster_idx_left_hand
    else:
        cluster_idx_subj = cluster_idx

    subj_eeg_blocks_orig = EEG_DATA_SHORT[subject_id]
    subj_lfp_blocks_orig = LFP_DATA_SHORT[subject_id]

    sfreq = subj_eeg_blocks_orig[0]['fs']
    channel_of_interest = 2 if subject_id in ['wue02', 'wue03'] else 1

    eeg_chan_names = subj_eeg_blocks_orig[0]['chan_names']
    lfp_chan_name = subj_lfp_blocks_orig[0]['chan_names'][channel_of_interest-1]

    subj_eeg_blocks = [trial['EEG_phases'] for block in CLEAN_DATA[subject_id] for trial in block]
    subj_lfp_blocks = [trial['LFP_phases'] for block in CLEAN_DATA[subject_id] for trial in block]

    # Collect all EEG and LFP trials
    all_eeg_phases = [trial for trial in subj_eeg_blocks]
    all_lfp_phases = [trial for trial in subj_lfp_blocks]

    n_phases = len(all_eeg_phases[0])

    # -------------------------
    # LOOP OVER PHASES
    # -------------------------
    for phase_idx in range(n_phases):

        print(f"\n--- Phase {phase_idx+1} ---")

        eeg_phase_trials = [trial[phase_idx][cluster_idx_subj] for trial in all_eeg_phases]
        lfp_phase_trials = [trial[phase_idx][channel_of_interest-1] for trial in all_lfp_phases]

        # Determine minimum length for trimming
        min_len = min(
            min(t.shape[1] for t in eeg_phase_trials),
            min(t.shape[0] if t.ndim == 1 else t.shape[1] for t in lfp_phase_trials)
        )

        # Function to trim trials to min_len
        def trim(t):
            t = t[np.newaxis, :] if t.ndim == 1 else t
            return t[:, :min_len]

        eeg_data = np.stack([trim(t) for t in eeg_phase_trials])
        lfp_data = np.stack([trim(t) for t in lfp_phase_trials])

        # -------------------------
        # CREATE DUMMY EPOCHS FOR CONNECTIVITY
        # -------------------------
        combined_ch_names = [lfp_chan_name] + [eeg_chan_names[i] for i in cluster_idx_subj]
        combined_ch_types = ['seeg'] + ['eeg'] * len(cluster_idx_subj)

        info = mne.create_info(combined_ch_names, sfreq, combined_ch_types)
        dummy = np.zeros((eeg_data.shape[0], len(cluster_idx_subj)+1, min_len))
        conn_obj = mne.EpochsArray(dummy, info, verbose=False)

        # -------------------------
        # BOOTSTRAP ITERATIONS
        # -------------------------
        for it in range(n_iterations):
            
            print(f"Subject: {subject_id}, Phase: {phase_idx+1}, Iteration: {it+1}/{n_iterations}")

            T_epoch = min_len / sfreq
            fmin = max(1.0, 5 / T_epoch)

            GC_lags = 40
            fmax = 100
            gc_L2E, gc_E2L, gc_net, freqs = bootstrap_iteration(
                eeg_data, lfp_data, min_len, sfreq, conn_obj,
                GC_lags, fmin, fmax
            )
            freqs = np.asarray(freqs)
            gc_net_mean_iter = gc_net.mean(axis=0) # mean over cluster channels
            gc_L2E_mean_iter = gc_L2E.mean(axis=0)
            gc_E2L_mean_iter = gc_E2L.mean(axis=0)

            # Initialize null arrays on first iteration
            if it == 0:
                freqs_ref = freqs
                null_L2E = np.zeros((n_iterations, gc_L2E.shape[0], len(freqs_ref)))
                null_E2L = np.zeros((n_iterations, gc_E2L.shape[0], len(freqs_ref)))
                null_NET = np.zeros((n_iterations, gc_net.shape[0], len(freqs_ref)))
                null_L2E_CLUST_MEAN = np.zeros((n_iterations, len(freqs_ref)))
                null_E2L_CLUST_MEAN = np.zeros((n_iterations, len(freqs_ref)))
                null_NET_CLUST_MEAN = np.zeros((n_iterations, len(freqs_ref)))

            null_L2E[it] = gc_L2E
            null_E2L[it] = gc_E2L
            null_NET[it] = gc_net
            null_L2E_CLUST_MEAN[it] = gc_L2E_mean_iter
            null_E2L_CLUST_MEAN[it] = gc_E2L_mean_iter
            null_NET_CLUST_MEAN[it] = gc_net_mean_iter

        # -------------------------
        # SAVE RESULTS
        # -------------------------
        GC_BOOTSTRAP_RESULTS[subject_id]['LFP2EEG'][f'phase_{phase_idx+1}'] = null_L2E
        GC_BOOTSTRAP_RESULTS[subject_id]['EEG2LFP'][f'phase_{phase_idx+1}'] = null_E2L
        GC_BOOTSTRAP_RESULTS[subject_id]['NET'][f'phase_{phase_idx+1}'] = null_NET
        GC_BOOTSTRAP_RESULTS[subject_id]['freqs'][f'phase_{phase_idx+1}'] = freqs
        GC_BOOTSTRAP_RESULTS[subject_id]['L2E_CLUST_MEAN'][f'phase_{phase_idx+1}'] = null_L2E_CLUST_MEAN
        GC_BOOTSTRAP_RESULTS[subject_id]['E2L_CLUST_MEAN'][f'phase_{phase_idx+1}'] = null_E2L_CLUST_MEAN
        GC_BOOTSTRAP_RESULTS[subject_id]['NET_CLUST_MEAN'][f'phase_{phase_idx+1}'] = null_NET_CLUST_MEAN


filename_gc = "GC_BOOTSTRAP_RESULTS_final_01.pkl"
full_path_gc = os.path.join(save_path, filename_gc)

with open(full_path_gc, 'wb') as f:
    pickle.dump(GC_BOOTSTRAP_RESULTS, f)

print(f"\n GC bootstrap results saved to {full_path_gc}")

