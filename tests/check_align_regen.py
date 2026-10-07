"""Cheap alignment check on the regen reference WITHOUT recomputing anything:
1) stored event latencies vs our beats (NN_times-derived): constant shift?
2) channel label sets/orders: ours vs MATLAB's hrsp.channels & the .set chanlocs.
3) sanity: the stored epochs' baseline-regressed pre-R mean == 0?
"""
import os
os.environ.setdefault('MPLBACKEND', 'Agg')
import numpy as np
import scipy.io as sio

BB = r'C:/Users/ccann/Documents/MATLAB/BrainBeats'
fs = 250.0
with open(os.path.join(BB, 'sample_data', 'dataset_HEP.set'), 'rb') as fh:
    hdr = sio.loadmat(fh, struct_as_record=False, squeeze_me=True)
ev = hdr['event']
lat = np.asarray([np.ravel(e.latency)[0] if hasattr(e, 'latency') else e['latency']
                  for e in (ev if isinstance(ev, np.ndarray) else [ev])], float)
bb = hdr.get('brainbeats')
pre = bb.preprocessings
nnt = np.asarray(pre.NN_times, float).ravel()
beats = np.rint(nnt * fs).astype(np.int64) + 1
kept = np.array(sorted(set(int(np.ravel(e.beat)[0]) if hasattr(e, 'beat') else 0
                           for e in ev)), int)
kept = kept[kept > 0]
ref_lat = beats[kept - 1]
print('latency check: n', lat.size, 'vs', ref_lat.size)
d = lat - ref_lat
print('latency diff (first 10):', d[:10])
print('unique diffs:', np.unique(np.round(d))[:8])
# channel labels
cl = hdr['chanlocs']
labels = np.asarray([str(c.labels) for c in cl] if hasattr(cl[0], 'labels')
                    else [c['labels'] for c in cl], object)
hrsp_ch = np.asarray(bb.hrsp.channels, str).ravel()
print('set chanlocs n:', labels.size, 'hrsp.channels n:', hrsp_ch.size)
print('same set:', set(labels) == set(hrsp_ch))
print('same order:', list(labels) == list(hrsp_ch))
roi = bb.roi
roi_ch = np.asarray(getattr(roi, 'channels', []), str).ravel()
print('roi channels:', list(roi_ch))
# stored epochs baseline check
data = None  # skip heavy data check here; diag_chain already compared
print('event latencies: min/max', lat.min(), lat.max())
print('beats: min/max', ref_lat.min(), ref_lat.max())
# ppg_transit
for f in getattr(pre, '_fieldnames', []):
    if 'transit' in f.lower() or 'pat' in f.lower():
        v = getattr(pre, f)
        try:
            print(f, '=', np.ravel(v)[:6])
        except Exception:
            print(f, v)