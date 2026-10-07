"""Dig into ppg_transit (what shifts did MATLAB apply?) and whether the regen
HEP looks physiologically sane.
"""
import os
os.environ.setdefault('MPLBACKEND', 'Agg')
import numpy as np
import scipy.io as sio

BB = r'C:/Users/ccann/Documents/MATLAB/BrainBeats'
with open(os.path.join(BB, 'sample_data', 'dataset_HEP.set'), 'rb') as fh:
    hdr = sio.loadmat(fh, struct_as_record=False, squeeze_me=True)
pt = hdr['brainbeats'].preprocessings.ppg_transit
for f in pt._fieldnames:
    v = getattr(pt, f)
    try:
        arr = np.asarray(v, float).ravel()
        print(f, arr.shape, arr[:8])
    except Exception:
        print(f, type(v))
    # nested?
    if hasattr(v, '_fieldnames'):
        print('  nested fields:', v._fieldnames)
# ROI HEP sanity: peak/trough timing
roi = hdr['brainbeats'].roi
h = np.asarray(getattr(roi, 'hep'), float).ravel()
t = np.asarray(getattr(roi, 'hep_times'), float).ravel()
ipk = int(np.argmax(np.abs(h)))
print('roi.hep |max| at t =', t[ipk], 's rel; sign', np.sign(h[ipk]))
print('roi.hep first 6 vals:', np.round(h[:6], 4))
print('roi.hep std:', h.std())
# old fixture (267-beat reference) for comparison
fit = np.load(os.path.join(BB, 'python', 'tests', 'data', 'reference_hep_tf.npz'),
              allow_pickle=True)
h2 = np.asarray(fit['hep_roi_from_epochs'], float)
t2 = np.asarray(fit['hep_times'], float)
ipk2 = int(np.argmax(np.abs(h2)))
print('OLD fixture hep: |max| at t =', t2[ipk2], '; std', h2.std())