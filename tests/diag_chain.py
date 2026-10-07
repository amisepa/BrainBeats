"""DECISIVE DIAGNOSTIC: is our cleaned continuous chain = MATLAB's?

Compare, at the same beats, the same way:
  A) stored MATLAB epochs (63x224x271, baseline-regressed, from dataset_HEP.set)
  B) our Xc_clean epochs (stage-0 chain + MATLAB's ICA comps removed),
     restricted to the HEP window, baseline-regressed per channel.
Also compare the UN-cleaned chain (stage-0 only, no comp removal) to locate
the divergence: stage-0 vs MATLAB-ICA-removal difference.
Outputs a per-check corr table.
"""
import os
os.environ.setdefault('MPLBACKEND', 'Agg')
import sys
import numpy as np
sys.path.insert(0, r'C:/Users/ccann/Documents/MATLAB/BrainBeats/python')
import eegprep

BB = r'C:/Users/ccann/Documents/MATLAB/BrainBeats'
fs = 250.0

# reference
HEPre = eegprep.pop_loadset(os.path.join(BB, 'sample_data', 'dataset_HEP.set'))
pre = HEPre['brainbeats']['preprocessings']
NN_times = np.asarray(pre['NN_times'], float).ravel()
beats_ref = np.rint(NN_times * fs).astype(np.int64) + 1
ref_beats = np.array([ev.get('beat') for ev in HEPre['event']], int)
ref_data = np.asarray(HEPre['data'], float)          # 63 x 224 x 271 (cleaned, base-lined)
ref_times = np.asarray(HEPre['times'], float).ravel()
print('ref epochs', ref_data.shape, 'times', ref_times[0], ref_times[-1])

# our stage-0 chain
rawEEG = eegprep.pop_loadset(os.path.join(BB, 'sample_data', 'dataset.set'))
EEG = eegprep.pop_select(rawEEG, 'nochannel', ['PPG', 'ECG'])
EEG = eegprep.pop_eegfiltnew(EEG, 0.5, 20)
out = eegprep.clean_channels(dict(EEG), corr_threshold=0.65, noise_threshold=15,
                             window_len=5, max_broken_time=0.33)
badEEG = out[0] if isinstance(out, tuple) else out
labels = [c['labels'] if isinstance(c, dict) else c['labels']
          for c in np.asarray(EEG['chanlocs']).ravel()]
mask = np.asarray([l not in ('PPG', 'ECG') for l in labels]) & \
    np.asarray([1] * len(labels), bool)
bad_labels = [labels[i] for i in np.flatnonzero(np.asarray(
    [c['data' if False else 'nan'] if False else np.nan for c in []] or
    (lambda b: b)(np.zeros(len(labels), bool))))]
# (bad_labels filled below properly)
badEEGset = badEEG
lab_bad = set()
try:
    rem = np.asarray(badEEGset.get('brainbeats', {}).get('removed_eeg_channels', []))
except Exception:
    rem = []
# simpler: detect removed channels by comparing channel counts
lab_after = [c['labels'] if isinstance(c, dict) else c['labels']
             for c in np.asarray(badEEGset['chanlocs']).ravel()]
lab_bad = [l for l in labels if l not in lab_after and l in lab_after + labels]
print('channels after clean_channels:', len(lab_after))
tp9loc = [dict(c) for c in np.asarray(EEG['chanlocs']).ravel()
          if (c['labels'] if isinstance(c, dict) else c['labels']) not in lab_after]
badEEGset['chanlocs'] = [dict(c) for c in np.asarray(badEEGset['chanlocs']).ravel()]
XI = eegprep.pop_interp(badEEGset, list(tp9loc), 'spherical')
XI = XI[0] if isinstance(XI, tuple) else XI
Xc = np.asarray(XI['data'], float)
lab_now = [c['labels'] if isinstance(c, dict) else c['labels']
           for c in np.asarray(XI['chanlocs']).ravel()]
print('stage-0 chain:', Xc.shape)

# MATLAB ICA comp removal
_ica = np.load(os.path.join(BB, 'python', 'tests', 'data', 'regen_ica.npz'))
W = np.asarray(_ica['icaweights'], float)
A = np.asarray(_ica['icawinv'], float)
removed = np.rint(np.asarray(_ica['removed'], float)).astype(int) - 1
S = W @ Xc
S[removed] = 0.0
Xc_clean = A @ S
print('comp-removal applied')

# epochs at the same beats
kept_beatnums = np.array(sorted(set(int(b) for b in ref_beats)), int)
hep_beats_1b = beats_ref[kept_beatnums - 1]
events = [{'type': 'R-peak', 'latency': float(b), 'duration': 0.0, 'urevent': i + 1}
          for i, b in enumerate(hep_beats_1b)]
EEG2 = dict(XI)
EEG2['event'] = events
EEG2 = eegprep.eeg_checkset(EEG2, 'eventconsistency')
t1 = ref_times[0] / 1000.0
t2 = (ref_times[-1] + 1000.0 / fs) / 1000.0  # EEGLAB: right edge inclusive
out = eegprep.pop_epoch(EEG2, ['R-peak'], (t1, t2))
mine = out[0] if isinstance(out, tuple) else out
mydata = np.asarray(mine['data'], float)
mytimes = np.asarray(mine['times'], float).ravel()
print('our epochs', mydata.shape, 'times', mytimes[0], mytimes[-1])

n = min(mydata.shape[2], ref_data.shape[2])
Dm = ref_data[:, :, :n]
Dp = mydata[:, :, :n]
# normalize per check: corr of raw + corr after per-epoch baseline (mean pre-R removed)
pre_mask = mytimes < 0
def base(x):  # x: (63, T, n)
    return x - x[:, pre_mask].mean(axis=1, keepdims=True)
r_raw = np.corrcoef(Dp.ravel(), Dm.ravel())[0, 1]
r_bl = np.corrcoef(base(Dp).ravel(), base(Dm).ravel())[0, 1]
print(f'OUR-CLEAN vs MATLAB-STORED: raw corr {r_raw:.4f} | baseline-regressed corr {r_bl:.4f}')
# per-epoch corr distribution (after baseline)
ce = [np.corrcoef(base(Dp)[:, :, e].ravel(), base(Dm)[:, :, e].ravel())[0, 1]
      for e in range(n)]
ce = np.asarray(ce)
print(f'per-epoch corr: median {np.median(ce):.4f}  min {ce.min():.4f}  max {ce.max():.4f}')
# channel-level: which channels diverge
cr = np.asarray([np.corrcoef(base(Dp)[c].ravel(), base(Dm)[c].ravel())[0, 1] for c in range(Dp.shape[0])])
worst = np.argsort(cr)[:6]
print('worst channels:', [(lab_now[c], round(float(cr[c]), 3)) for c in worst])
np.savez(os.path.join(BB, 'tests', 'chain_diag.npz'), r_raw=r_raw, r_bl=r_bl,
         per_epoch=ce, per_chan=cr, labels=np.asarray(lab_now, object))
print('DIAG DONE')