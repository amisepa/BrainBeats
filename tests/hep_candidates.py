"""Which preprocessing does the reference HEP apply before the beat mean?

Stored-truth = data[roi].mean(0).mean(1) should equal our hand mean IF the
stored HEP is the plain mean. It isn't (corr ~0). Candidates from run_HEP.m:
  (a) per-epoch baseline_regression (Alday) before/after averaging
  (b) HEP computed on the CONTINUOUS cleaned signal (crop + T applied) whose
      samples differ from the epoch-concat by edge effects of pop_epoch pad
  (c) per-beat normalisation/division
  (d) ROI computed as mean over channels of PER-CHANNEL HEPs, each scaled
Test (d)-style alternatives and the baseline_regression quickly.
"""
import os
os.environ.setdefault('MPLBACKEND', 'Agg')
import sys
import numpy as np
sys.path.insert(0, r'C:/Users/ccann/Documents/MATLAB/BrainBeats/python')
import eegprep

HEP = eegprep.pop_loadset(r'C:/Users/ccann/Documents/MATLAB/BrainBeats/sample_data/dataset_HEP.set')
br = HEP['brainbeats']
hrsp_ref = br['hrsp']
data = np.asarray(HEP['data'], float)
labels = [c['labels'] if isinstance(c, dict) else c['labels']
          for c in np.asarray(HEP['chanlocs']).ravel()]
roi = ['F1','F2','F3','F4','Fz','FC1','FC2','FC3','FC4','FCz',
       'C1','C2','C3','C4','Cz']
ridx = [i for i, l in enumerate(labels) if l in roi]
hep_hand = data[ridx].mean(axis=0).mean(axis=1)
hep_times_ref = np.asarray(hrsp_ref['hep_times'], float)

truth = np.asarray(hrsp_ref['hep'][0]
                   if 0 else hrsp_ref['hep'], float)
print('truth shape', truth.shape, 'hand', hep_hand.shape)
t = truth[0] if truth.ndim > 1 else truth
print('corr(truth, plain-mean)',
      np.corrcoef(t, hep_hand)[0, 1])
idx0 = np.rint(hep_times_ref / 1000.0 * 250.0).astype(int)
# (e) hep of the ROI MEDIAN over channels?
med = data[ridx].mean(axis=0)
print('corr(truth, mean-epochs-then-chan-mean)', np.corrcoef(t, med.mean(axis=1))[0, 1])
# (f) maybe the stored data epochs were BASELINE-CORRECTED only (data stored is
# post-cleaning but pre-baseline) and HEP had baseline regression:
try:
    from functions.baseline_regression import baseline_regression
    ep = data[ridx]                          # (15, 224, 267)
    out = np.empty_like(ep)
    bsl = (-300.0, 0.0)
    for c in range(ep.shape[0]):
        for b in range(ep.shape[2]):
            out[c, :, b] = baseline_regression(
                ep[c, :, b], 250.0, bsl, np.asarray(HEP['times']).ravel() / 1000.0 * 0 + 0)
    hep_br = out.mean(axis=0).mean(axis=1)
    print('corr(truth, baseline-then-mean)', np.corrcoef(t, hep_br)[0, 1])
except Exception as e:
    print('baseline_regression path failed:', type(e).__name__, e)
# (g) check 'hep' in hrsp_ref keys
print('hrsp keys:', sorted(br['hrsp'].keys()) if isinstance(br['hrsp'], dict) else type(br['hrsp']))
print('preprocessings keys:', sorted(br['preprocessings'].keys())
      if isinstance(br['preprocessings'], dict) else type(br['preprocessings']))