"""Stage B (definitive): parity using the REFERENCE beat train exactly.

The reference run's beat train is stored in dataset_HEP.set:
  brainbeats.preprocessings.NN_times (301 beat times in s)
  brainbeats.preprocessings.bad_heartbeats + interpolated_heartbeats (per-beat flags)
  event .beat values -> which of the 301 beats produced each of the 267 epochs
No re-detection: beats = round(NN_times * fs) + 1 (EEGLAB 1-based).

Continuous cleaned data: same chain as the reference stage-0 (0.5-30 Hz filter,
full-rank CAR, clean_channels -> interp back), then the epoch-level ICA
cleaning chain reproduced with Picard (the reference's ica_method; Infomax is
the going-forward default but weights can't match across implementations).
Xc_clean = T Xc with T = Yc Xr' pinv(Xr Xr') (epochs before/after cleaning).
Removed trials (removed_eeg_trials, 15) are applied when building Xr/Yc.

Compare HRSP/HRPC against brainbeats.hrsp; the reference HEP (all channels
and ROI) is not stored in this .set, so HEP parity is checked against the
stored epochs separately (they are the baseline-regressed HEP of the SAME
chain -- the HEP-parity test uses our epoch-level chain + regression).
"""
import os
os.environ.setdefault('MPLBACKEND', 'Agg')
import sys
import numpy as np
sys.path.insert(0, r'C:/Users/ccann/Documents/MATLAB/BrainBeats/python')
import eegprep
from functions.compute_hep_tf import compute_hep_tf

BB = r'C:/Users/ccann/Documents/MATLAB/BrainBeats'
fs = 250.0

HEPre = eegprep.pop_loadset(os.path.join(BB, 'sample_data', 'dataset_HEP.set'))
pre = HEPre['brainbeats']['preprocessings']
NN_times = np.asarray(pre['NN_times'], float).ravel()
bad_hb = np.asarray(pre['bad_heartbeats'], bool).ravel()
interp_hb = np.asarray(pre['interpolated_heartbeats'], bool).ravel()
removed_trials = np.asarray(pre['removed_eeg_trials'], int).ravel()
n_train = NN_times.size                      # 301 beats
beats_ref = np.rint(NN_times * fs).astype(np.int64) + 1   # 1-based samples
print('train:', n_train, 'bad:', bad_hb.sum(), 'interp:', interp_hb.sum(),
      'removed trials:', removed_trials.size)
ref_beats = np.array([ev.get('beat') for ev in HEPre['event']], int)   # into 1..301
print('epoch beat numbers:', ref_beats.size, ref_beats.min(), ref_beats.max())

# heart channel used by the reference: PPG
rawEEG = eegprep.pop_loadset(os.path.join(BB, 'sample_data', 'dataset.set'))
labels_raw = [c['labels'] if isinstance(c, dict) else c['labels']
              for c in np.asarray(rawEEG['chanlocs']).ravel()]
print('labels_raw has PPG:', 'PPG' in labels_raw)

# --- continuous cleaned EEG (stage 0) ---------------------------------------
EEG = eegprep.pop_select(rawEEG, 'nochannel', ['PPG', 'ECG'])
EEG = eegprep.pop_eegfiltnew(EEG, 0.5, 20)      # reference hep_ppg: lowpass 20
# ref='infinity': no CAR, no reref (reference hep_ppg params)
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
Xc = np.asarray(XI['data'], float)
lab_now = [c['labels'] if isinstance(c, dict) else c['labels']
           for c in np.asarray(XI['chanlocs']).ravel()]
print('channels:', Xc.shape, 'TP9 back:', 'TP9' in lab_now)

# --- epoch-level cleaning (the T), mirroring the reference ------------------
win = (-300.0, 600.0)
pad = 650.0
kept_train = np.flatnonzero(~bad_hb & ~interp_hb)          # 0-based into train
# beats producing the STORED epochs: .beat numbers, consecutive-ish; removed
# trials were dropped AFTER epoching (267 stored = 267 kept beats? verify)
kept_beatnums = np.array(sorted(set(int(b) for b in ref_beats)), int)  # stored epochs' beats
print('stored epochs come from', kept_beatnums.size, 'distinct beats')
hep_beats_1b = beats_ref[kept_beatnums - 1]                # 1-based sample idx
events = [{'type': 'R-peak', 'latency': float(b), 'duration': 0.0,
           'urevent': i + 1}
          for i, b in enumerate(hep_beats_1b)]
EEG2 = dict(XI)
EEG2['event'] = events
EEG2 = eegprep.eeg_checkset(EEG2, 'eventconsistency')
out = eegprep.pop_epoch(EEG2, ['R-peak'], ((win[0] - pad) / 1000.0,
                                           (win[1] + pad) / 1000.0))
rawwide = out[0] if isinstance(out, tuple) else out
rawdata = np.asarray(rawwide['data'], float)
n_ep = rawdata.shape[2]
print('padded epochs:', rawdata.shape)
# apply removed_eeg_trials? ref event .beat values already EXCLUDE removed
# epochs only if stored epochs are post-removal; try BOTH and keep the better.

# ICA: the reference's OWN solution (no Picard here)
_ica = np.load(os.path.join(BB, 'python', 'tests', 'data', 'regen_ica.npz'))
W = np.asarray(_ica['icaweights'], float)          # ncomp x nch
A = np.asarray(_ica['icawinv'], float)             # nch x ncomp
_removed = np.rint(np.asarray(_ica['removed'], float)).astype(int) - 1
print('MATLAB ICA loaded:', W.shape, 'removed comps:', (_removed + 1).tolist())
S = W @ Xc
S[_removed] = 0.0
Xc_clean = A @ S
print('component removal applied; Xc_clean', Xc_clean.shape)


# --- our ROI HEP from the cleaned continuous chain --------------------------
roi_chans = [str(c) for c in np.asarray(
    HEPre['brainbeats']['roi']['channels'], str).ravel()] \
    if hasattr(HEPre['brainbeats'].get('roi'), '__getitem__') or True else []
lab_now = [c['labels'] if isinstance(c, dict) else c['labels']
           for c in np.asarray(XI['chanlocs']).ravel()]
ridx = [lab_now.index(ch) for ch in roi_chans if ch in lab_now]
print('ROI channels found:', len(ridx), '/', len(roi_chans))
# epochs from the cleaned continuous data at the SAME beats (with padding)
EEGc = dict(XI)
EEGc['data'] = Xc_clean
EEGc['event'] = events
EEGc = eegprep.eeg_checkset(EEGc, 'eventconsistency')
outc = eegprep.pop_epoch(EEGc, ['R-peak'], ((win[0] - pad) / 1000.0,
                                            (win[1] + pad) / 1000.0))
cep = outc[0] if isinstance(outc, tuple) else outc
cdat = np.asarray(cep['data'], float)              # roi? no: (63, npnts, n)
ctimes = np.asarray(cep['times'], float).ravel()
pre_mask = ctimes < 0                              # baseline window (pre-R)
cdat = cdat[ridx]
# restrict to the HEP window (reference -300..600) before baseline regression
roi_times_ref = np.asarray(HEPre['brainbeats']['roi']['hep_times'], float).ravel()
wmask = (ctimes >= roi_times_ref[0] - 1e-6) & (ctimes <= roi_times_ref[-1] + 1e-6)
cdat = cdat[:, wmask]; ctimes_w = ctimes[wmask]
bl = cdat[:, ctimes_w < 0].mean(axis=1, keepdims=True)
cdat = cdat - bl                                   # baseline regression per channel
hep_py_roi = cdat.mean(axis=2).mean(axis=0)        # mean trials + mean ROI chs
if ctimes_w.size != roi_times_ref.size:
    hep_py_roi = np.interp(roi_times_ref, ctimes_w, hep_py_roi)
print('roi grids: ours', ctimes.size, 'ref', roi_times_ref.size)

# --- measures with the REFERENCE beats --------------------------------------
HRSPr = HEPre['brainbeats']['hrsp']
hrsp_ref = np.asarray(HRSPr['hrsp'], float)
hepc_ref = np.asarray(HRSPr.get('hepc', HRSPr.get('hrpc')), float)
power_ref = np.asarray(HRSPr['power'], float)
hep_times = np.asarray(HRSPr['hep_times'], float)
times_ref = np.asarray(HRSPr['times'], float)
freqs_ref = np.asarray(HRSPr['freqs'], float)
lab_ref = [str(c) for c in np.asarray(HRSPr['channels']).ravel()]
print('ref channel order == ours:', lab_ref == lab_now)

tf_all, _ = compute_hep_tf(Xc_clean, fs, hep_beats_1b, win,
                           hep_times=hep_times, freqs=freqs_ref, n_surr=0)
hrsp_py, hrpc_py = tf_all['hrsp'], tf_all['hrpc']
if lab_ref != lab_now:
    order = [lab_now.index(l) for l in lab_ref]
    hrsp_py = hrsp_py[order]
    hrpc_py = hrpc_py[order]
r_b = np.corrcoef(hrsp_py.ravel(), hrsp_ref.ravel())[0, 1]
r_c = np.corrcoef(hrpc_py.ravel(), hepc_ref.ravel())[0, 1]
d_b = np.abs(hrsp_py - hrsp_ref)
print(f'HRSP: corr {r_b:.6f}  max|dB| {d_b.max():.3f}  median {np.median(d_b):.4f}')
print(f'HRPC: corr {r_c:.6f}  max {np.abs(hrpc_py - hepc_ref).max():.4f}')
# ROI HEP parity (regen reference stores brainbeats.roi.hep)
_roi = HEPre['brainbeats'].get('roi') or {}
roi_hep_ref = np.asarray(_roi.get('hep'), float).ravel() \
    if _roi.get('hep') is not None else np.zeros(0)
if roi_hep_ref.size:
    r_h = np.corrcoef(hep_py_roi.ravel(), roi_hep_ref.ravel())[0, 1]
    print(f'ROI HEP: corr {r_h:.6f}  max|diff| {np.abs(hep_py_roi - roi_hep_ref).max():.5f}')
np.savez(os.path.join(BB, 'tests', 'parity_result.npz'),
         r_hrsp=r_b, r_hrpc=r_c, max_db=float(d_b.max()),
         med_db=float(np.median(d_b)))
print('STAGE-B DONE')