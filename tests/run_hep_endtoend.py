"""Smoke/parity driver for the run_hep.py port: brainbeats_process's HEP branch
in Python. Chain: load dataset.set -> remove heart channels -> stage-0 clean
(0.5-30, CAR, clean_channels+interp) -> get the beat train (the REFERENCE
NN_times, since our PPG detection needs work; documented as the known gap) ->
run_hep(EEG, CARDIO, params, Rpeaks) with Infomax (params icamethod=2) ->
save dataset_HEPpy.set -> compare against dataset_HEP.set.

ICAMETHOD: params.icamethod=2 -> extended Infomax via eegprep; on the sample
data that took >4 core-hours without converging (log stageb2.log), so this
driver uses icamethod=1 (Picard, matching the MATLAB reference) by default
and prints that choice. Override with ICAMETHOD=2 to run Infomax.
"""
import os
os.environ.setdefault('MPLBACKEND', 'Agg')
import sys
import numpy as np
sys.path.insert(0, r'C:/Users/ccann/Documents/MATLAB/BrainBeats/python')
import eegprep
from functions.run_hep import run_hep
from functions.params import default_params

BB = r'C:/Users/ccann/Documents/MATLAB/BrainBeats'
fs = 250.0
ICAMETHOD = int(os.environ.get('ICAMETHOD', '1'))

EEG = eegprep.pop_loadset(os.path.join(BB, 'sample_data', 'dataset.set'))
CARDIO = eegprep.pop_select(dict(EEG), 'channel', ['PPG', 'ECG'])
CARDIO = CARDIO[0] if isinstance(CARDIO, tuple) else CARDIO
EEG = eegprep.pop_select(EEG, 'nochannel', ['PPG', 'ECG'])
EEG = EEG[0] if isinstance(EEG, tuple) else EEG
# The reference (hep_ppg test of run_tutorial_headless): 0.5-20 Hz CAUSAL FIR,
# ref = infinity (no reref), clean_eeg = TRUE, detectMethod 'median',
# icamethod 1 (Picard), heart PPG, hep_window adaptive.
EEG = eegprep.pop_eegfiltnew(EEG, 0.5, 20)
EEG = EEG[0] if isinstance(EEG, tuple) else EEG
out = eegprep.clean_channels(dict(EEG), corr_threshold=0.65, noise_threshold=15,
                             window_len=5, max_broken_time=0.33)
# no CAR, no re-referencing: the reference ran ref='infinity'
badEEG = out[0] if isinstance(out, tuple) else out
mask = np.asarray(badEEG['etc']['clean_channel_mask'], bool)
labels = [c['labels'] if isinstance(c, dict) else c['labels']
          for c in np.asarray(EEG['chanlocs']).ravel()]
bad_labels = [labels[i] for i in np.flatnonzero(~mask)]
print('stage-0 bad channels:', bad_labels)
tp9loc = [dict(c) for c in np.asarray(EEG['chanlocs']).ravel()
          if (c['labels'] if isinstance(c, dict) else c['labels']) in bad_labels]
badEEG['chanlocs'] = [dict(c) for c in np.asarray(badEEG['chanlocs']).ravel()]
XI = eegprep.pop_interp(badEEG, list(tp9loc), 'spherical')
EEG = XI[0] if isinstance(XI, tuple) else XI

# beat train: our PPG train (the reference's NN_times were lost with the
# overwritten dataset_HEP.set; parity numbers come from parity_stage_b3 which
# reconstructs the exact reference train from the fixture)
Rpeaks = np.load(os.path.join(BB, 'tests', 'our_train_ppg.npy'))
HEPre = eegprep.pop_loadset(os.path.join(BB, 'python', 'tests', 'data',
                                         'reference_hep_tf.npz'), allow_pickle=True)     if False else None

params = default_params(save=True, hep_surrogates=100, icamethod=ICAMETHOD,
                        heart_signal='PPG', heart_channels=['PPG'],
                        detectMethod='median', clean_eeg=True)
HEP = run_hep(EEG, CARDIO, params, Rpeaks)
print('run_hep OK: HEP', HEP['data'].shape if 'data' in HEP else HEP,
      'trials', HEP['trials'])
print('SAVED:', str(EEG.get('filename'))[:-4] + '_HEP.set')

# compare the in-memory HEP (brainbeats included) against the fixture
_fac = np.load(os.path.join(BB, 'python', 'tests', 'data', 'reference_hep_tf.npz'),
               allow_pickle=True)
hrsp_py = np.asarray(HEP['brainbeats']['hrsp']['hrsp'], float)
hepc_py = np.asarray(HEP['brainbeats']['hrsp'].get('hrpc',
           HEP['brainbeats']['hrsp'].get('hepc')), float)
hrsp_ref = np.asarray(_fac['hrsp'], float)
hepc_ref = np.asarray(_fac['hrpc'], float)
if hrsp_py.shape == hrsp_ref.shape:
    r_b = np.corrcoef(hrsp_py.ravel(), hrsp_ref.ravel())[0, 1]
    r_c = np.corrcoef(hepc_py.ravel(), hepc_ref.ravel())[0, 1]
    d_b = np.abs(hrsp_py - hrsp_ref)
    print(f'END-TO-END HRSP corr {r_b:.4f}  max|dB| {d_b.max():.3f} '
          f'median {np.median(d_b):.4f}')
    print(f'END-TO-END HRPC corr {r_c:.4f}')
else:
    print('shapes differ (different beat source):',
          hrsp_py.shape, hrsp_ref.shape)
    n = min(hrsp_py.shape[-1], hrsp_ref.shape[-1])
    m = min(hrsp_py.shape[0], hrsp_ref.shape[0])
    r_b = np.corrcoef(hrsp_py[:m, :, :n].ravel(), hrsp_ref[:m, :, :n].ravel())[0, 1]
    r_c = np.corrcoef(hepc_py[:m, :, :n].ravel(), hepc_ref[:m, :, :n].ravel())[0, 1]
    print(f'END-TO-END (first {n} TF steps): HRSP corr {r_b:.4f}  HRPC corr {r_c:.4f}')
