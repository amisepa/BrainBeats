"""PARAMS - parameter defaults mirroring getparams_cmd.m / run_checks.m.

One place for every default the MATLAB sets (run_checks.m fills fs, icamethod,
ref; clean_eeg.m fills highpass/lowpass/ref/ASR/ICA thresholds; run_HEP.m fills
the HEP window and ROI; compute_hep_tf.m the TF grid).
"""
from __future__ import annotations


class Params(dict):
    """Dict with attribute access; unknown keys are added by callers like params.x = v."""

    def __getattr__(self, k):
        try:
            return self[k]
        except KeyError as e:
            raise AttributeError(k) from e

    def __setattr__(self, k, v):
        self[k] = v


def default_params(**over):
    """Defaults matching the MATLAB command-line path with analysis='hep'.

    brainbeats_process(EEG, 'analysis','hep', 'heart_signal','ECG',
                       'heart_channels',{'ECG'}, ...) equivalents.
    Values that MATLAB fills in run_checks/clean_eeg placeholders are noted.
    """
    p = Params(
        analysis='hep',
        heart_signal='ecg',          # 'ecg' | 'ppg' | 'rr' | 'off'
        heart_channels=['ECG'],      # labels of the heart channel(s)
        # --- clean_eeg (stage 0) ---
        clean_eeg=True,
        highpass=0.5,                # HEP keeps slow components; ICA 1-Hz copy
        lowpass=30,
        filttype='noncausal',        # 'noncausal' (zero-phase) | 'causal'
        ref='average',               # 'average' | 'off' | 'csd' | 'infinity'
        linenoise=0,                 # 0 = no notch; 50|60 adds a band-stop
        flatline=5, corrThresh=0.65, maxBad=0.33,  # channel cleaning
        # --- clean_eeg (stage 1) ---
        detectMethod='grubbs',       # bad epochs for HEP
        asr_cutoff=50,               # unused by HEP (find_badTrials instead)
        icamethod=2,                 # 1 picard | 2 infomax | 3 replicable infomax
        conf_thresh=0.75,            # ICLabel heart-probability removal threshold
        clean_method='asr_ica',      # 'asr_ica' | 'gedai' (gedai not ported yet)
        # --- run_HEP ---
        hep_window=(-300, 600),      # ms, or 'adaptive'
        hep_level='channels',        # 'channels' | 'ics' | 'both'
        hep_roi=['F1','F2','F3','F4','Fz','FC1','FC2','FC3','FC4','FCz',
                 'C1','C2','C3','C4','Cz'],
        hep_baseline='none',         # 'none' | 'regression'
        hep_baseline_win=(-150, -50),
        hep_tf_freqs=(4, 30),
        hep_tf_cycles=5,
        hep_tf_tstep=10,             # ms output step
        hep_surrogates=0,            # 0 = none
        hep_surrogate_mode='shuffle',  # 'shuffle' | 'rigid'
        hep_surrogate_shift=(-500, 500),  # 'rigid'
        surrogate_seed=1,
        keep_heart=False,
        # --- run_checks ---
        fs=None,                     # filled from EEG.srate
        vis_cleaning=False,
        vis_outputs=False,
        save=True,
        hrv_features=False,
        eeg_features=False,
        parpool=False,
        gong=False,
    )
    p.update(over)
    return p