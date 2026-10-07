"""Stage B (final): HEP/HRSP/HRPC parity against reference_hep_tf.npz.

Beat source: our validated get_rr + clean_rr train (tests/our_train.npy,
305 beats, 1-based). The reference train had 301; our get_rr found 2 extra
detections (identified by the NN-series alignment below), so the parity is
tolerance-based as agreed.

The reference beat SET for the epochs is known exactly: dataset_HEP.set's
event .beat values (1-based into the 301-beat reference train). We map
reference beat k -> our position by walking the IBI series, skipping extras.

Cleaning chain reproduced (icamethod 2 = eegprep-validated Infomax):
  stage 0: 0.5-30 Hz -> full-rank CAR -> clean_channels -> pop_interp BY LABEL
  stage 1 (on the padded epochs): 1-Hz copy ICA (extended Infomax at rank)
  -> ICLabel -> icflag(.99/.9/.75/.99/.99) -> pop_subcomp -> crop
  -> T = Yc Xr' pinv(Xr Xr') applied to the continuous data.
"""
import os
os.environ.setdefault('MPLBACKEND', 'Agg')
import sys
import numpy as np
sys.path.insert(0, r'C:/Users/ccann/Documents/MATLAB/BrainBeats/python')
sys.path.insert(0, r'C:/Users/ccann/Documents/HEP_neurofeedback/python')
import eegprep
from functions.compute_hep_tf import compute_hep_tf

BB = r'C:/Users/ccann/Documents/MATLAB/BrainBeats'
fs = 250.0

# ---- continuous cleaned data (stage A, corrected interp) --------------------
EEG = eegprep.pop_loadset(os.path.join(BB, 'sample_data', 'dataset.set'))
EEG = eegprep.pop_select(EEG, 'nochannel', ['PPG', 'ECG'])
EEG = eegprep.pop_eegfiltnew(EEG, 0.5, 30)
X = np.asarray(EEG['data'], float)
X0 = np.vstack([X, np.zeros((1, X.shape[1]))])      # full-rank CAR
X0 = X0 - X0.mean(axis=0, keepdims=True)
EEG['data'] = X0[:-1]
out = eegprep.clean_channels(dict(EEG), corr_threshold=0.65, noise_threshold=15,
                             window_len=5, max_broken_time=0.33)
badEEG = out[0] if isinstance(out, tuple) else out
mask = np.asarray(badEEG['etc']['clean_channel_mask'], bool)
labels = [c['labels'] if isinstance(c, dict) else c['labels']
          for c in np.asarray(EEG['chanlocs']).ravel()]
bad_labels = [labels[i] for i in np.flatnonzero(~mask)]
print('bad channels:', bad_labels)
tp9loc = [dict(c) for c in np.asarray(EEG['chanlocs']).ravel()
          if (c['labels'] if isinstance(c, dict) else c['labels']) in bad_labels]
badEEG['chanlocs'] = [dict(c) for c in np.asarray(badEEG['chanlocs']).ravel()]
XI = eegprep.pop_interp(badEEG, list(tp9loc), 'spherical')
XI = XI[0] if isinstance(XI, tuple) else XI
XI = XI[0] if isinstance(XI, tuple) else XI
Xc = np.asarray(XI['data'], float)
lab_now = [c['labels'] if isinstance(c, dict) else c['labels']
           for c in np.asarray(XI['chanlocs']).ravel()]
print('channels after interp:', Xc.shape[0], 'TP9 back:', 'TP9' in lab_now)

# ---- beat alignment ----------------------------------------------------------
our_train = np.load(os.path.join(BB, 'tests', 'our_train.npy')).astype(np.int64)
HEPre = eegprep.pop_loadset(os.path.join(BB, 'sample_data', 'dataset_HEP.set'))
ref_beats = np.array([ev.get('beat') for ev in HEPre['event']], float)
ref_NN = np.asarray(HEPre['brainbeats']['preprocessings']['NN'], float)
our_NN = np.diff(our_train) / fs
print('ref NN n', ref_NN.size, 'our NN n', our_NN.size)
n_ref, n_our = ref_NN.size, our_NN.size
# align IBI series from the END backwards, marking OUR intervals that have no
# reference partner (our extra detections split a reference IBI in two)
TOL = 5.1e-3                     # 2 samples at 250 Hz
i_ref, i_our = n_ref - 1, n_our - 1
extra_iv = set()
while i_ref >= 0 and i_our >= 0:
    if abs(ref_NN[i_ref] - our_NN[i_our]) <= TOL:
        i_ref -= 1
        i_our -= 1
        continue
    # decide: our extra interval, or drift? look ahead 1 step in our series
    if i_our - 1 >= 0 and abs(ref_NN[i_ref] - our_NN[i_our - 1]) <= TOL:
        extra_iv.add(i_our)
        i_our -= 1
    elif i_ref - 1 >= 0 and abs(ref_NN[i_ref - 1] - our_NN[i_our]) <= TOL:
        i_ref -= 1                # ref split: shouldn't happen; ref is shorter
    else:
        # neither: allow a 2-sample drift and keep both (matched-with-jitter)
        extra_iv.add(i_our) if False else None
        i_ref -= 1
        i_our -= 1
extra_beats = sorted(i + 1 for i in extra_iv)   # 1-based OUR beat numbers
unmatched = i_ref + 1
print('our extra beats (1-based):', extra_beats, ' unmatched ref intervals:', unmatched)
# map reference beat number (1-based) -> our position by walking FORWARD,
# skipping our extra intervals
ref2our = {}
j = 0
for kref in range(2, n_ref + 1):
    while j < n_our and (j + 1) in extra_beats:
        j += 1
    if j < n_our:
        ref2our[kref] = j + 1     # our 1-based beat number ending matched IBI j
    j += 1
kept_ref_nums = np.unique(ref_beats.astype(int))
sel = [int(our_train[ref2our[k] - 1]) for k in kept_ref_nums if k in ref2our]
hep_beats_1b = np.array(sel, np.int64)
print('our beats for the HEP epochs:', hep_beats_1b.size, 'of ref', kept_ref_nums.size)

# ---- epoch-level cleaning + measures ----------------------------------------
win = (-300.0, 600.0)
pad = 650.0
events = [{'type': 'R-peak', 'latency': float(b), 'duration': 0.0,
           'urevent': i + 1, 'beat': int(k)}
          for i, (k, b) in enumerate(zip(kept_ref_nums, hep_beats_1b))]
EEG2 = dict(XI)
EEG2['event'] = events
EEG2 = eegprep.eeg_checkset(EEG2, 'eventconsistency')
out = eegprep.pop_epoch(EEG2, ['R-peak'], ((win[0] - pad) / 1000.0,
                                           (win[1] + pad) / 1000.0))
HEPwide = out[0] if isinstance(out, tuple) else out
rawwide_arr = np.asarray(HEPwide['data'], float)
print('padded epochs:', rawwide_arr.shape)

HEPica = eegprep.pop_eegfiltnew(dict(HEPwide), locutoff=1.0)
HEPica = HEPica[0] if isinstance(HEPica, tuple) else HEPica
dataI = np.asarray(HEPica['data'], float)
rank = int((np.linalg.eigvalsh(np.cov(dataI.reshape(dataI.shape[0], -1))) > 1e-7).sum())
print('rank:', rank)
ica = eegprep.pop_runica(dict(HEPica), icatype='picard', maxiter=500, mode='standard')      # MATLAB reference used Picard (ica_method 1)
ica = ica[0] if isinstance(ica, tuple) else ica
HEPwide['icaweights'] = ica['icaweights']
HEPwide['icasphere'] = ica['icasphere']
HEPwide['icawinv'] = ica['icawinv']
HEPwide['icachansind'] = ica['icachansind']
HEPwide = eegprep.eeg_checkset(HEPwide)
out = eegprep.pop_iclabel(dict(HEPwide), 'default')
HEPwide = out[0] if isinstance(out, tuple) else out
thr = np.array([[np.nan, np.nan], [0.99, 1.], [0.9, 1.], [0.75, 1.],
                [0.99, 1.], [0.99, 1.], [np.nan, np.nan]])
out = eegprep.pop_icflag(dict(HEPwide), thr)
HEPwide = out[0] if isinstance(out, tuple) else out
bad_comp = np.flatnonzero(np.asarray(HEPwide['reject']['gcompreject']))
print('bad components (1-based):', (bad_comp + 1).tolist())
if bad_comp.size:
    out = eegprep.pop_subcomp(dict(HEPwide), (bad_comp + 1).tolist(), 0)
    HEPwide = out[0] if isinstance(out, tuple) else out

# ---- continuous cleaning transform: T = Yc Xr' pinv(Xr Xr') -----------------
Xr = rawwide_arr.reshape(rawwide_arr.shape[0], -1)
YcA = np.asarray(HEPwide['data'], float)
Yc = YcA.reshape(YcA.shape[0], -1)
T = (Yc @ Xr.T) @ np.linalg.pinv(Xr @ Xr.T)
err = np.linalg.norm(Yc - T @ Xr, 'fro') / max(np.linalg.norm(Yc, 'fro'), 1e-30)
print(f'clean transform {T.shape}, rel err {err:.1e}')
Xc_clean = T @ Xc

# ---- measures --------------------------------------------------------------
ref = np.load(os.path.join(BB, 'python', 'tests', 'data', 'reference_hep_tf.npz'),
              allow_pickle=True)
roi_labels = list(ref['roi_channels'])
ridx = [i for i, l in enumerate(lab_now) if l in roi_labels]
print('roi matched:', len(ridx), 'of', len(roi_labels), 'labels')
sig_roi = Xc_clean[ridx].mean(axis=0)
hep_times = np.asarray(ref['hep_times'])
tf_all, _ = compute_hep_tf(Xc_clean, fs, hep_beats_1b, win, hep_times=hep_times,
                           n_surr=0)
tf_roi, _ = compute_hep_tf(sig_roi[None, :], fs, hep_beats_1b, win,
                           hep_times=hep_times, n_surr=0)
hrsp_py, hrpc_py = tf_all['hrsp'], tf_all['hrpc']
hrsp_ref, hrpc_ref = ref['hrsp'], ref['hrpc']
lab_ref = list(ref['all_channels'])
if lab_now != lab_ref:
    order = [lab_now.index(l) for l in lab_ref]
    hrsp_py = hrsp_py[order]
    hrpc_py = hrpc_py[order]
r_b = np.corrcoef(hrsp_py.ravel(), hrsp_ref.ravel())[0, 1]
r_c = np.corrcoef(hrpc_py.ravel(), hrpc_ref.ravel())[0, 1]
d_b = np.abs(hrsp_py - hrsp_ref)
d_c = np.abs(hrpc_py - hrpc_ref)
print(f'HRSP: corr {r_b:.6f}  max|dB| {d_b.max():.3f}  median {np.median(d_b):.4f}')
print(f'HRPC: corr {r_c:.6f}  max {d_c.max():.4f}  median {np.median(d_c):.4f}')
hep_py = tf_roi['hep'][0]
hep_ref_roi = ref['hep_roi_from_epochs']
r_h = np.corrcoef(hep_py, hep_ref_roi)[0, 1]
print(f'HEP ROI: corr {r_h:.6f}  max|diff| {np.abs(hep_py - hep_ref_roi).max():.2e} uV')
np.savez(os.path.join(BB, 'tests', 'parity_result.npz'),
         r_hrsp=r_b, r_hrpc=r_c, max_db=float(d_b.max()),
         med_db=float(np.median(d_b)), max_pc=float(d_c.max()), r_hep=float(r_h))
print('PARITY DONE')
