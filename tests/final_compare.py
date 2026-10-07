"""Final end-to-end comparison against the fixture reference (the MATLAB
63x224x267 dataset_HEP.set that existed before the overwrite; its hrsp struct
was saved in python/tests/data/reference_hep_tf.npz).

Our end-to-end output is the current sample_data/dataset_HEP.set (Python run).
Grid: ours 226 pts (-300..597?), reference 224 (-300..592) -- the HEP grids
differ only at the tail (our win passed 600: -300..600 -> 226 pts, reference
used 592). For HRSP/HRPC the TF grids (75 pts) should match exactly.
"""
import os
os.environ.setdefault('MPLBACKEND', 'Agg')
import sys
import numpy as np
sys.path.insert(0, r'C:/Users/ccann/Documents/MATLAB/BrainBeats/python')
import eegprep

ref = np.load(r'C:/Users/ccann/Documents/MATLAB/BrainBeats/python/tests/data/reference_hep_tf.npz',
              allow_pickle=True)
PY = eegprep.pop_loadset(
    r'C:/Users/ccann/Documents/MATLAB/BrainBeats/sample_data/dataset_HEP.set')
hp = PY['brainbeats']['hrsp']
hrsp_py = np.asarray(hp['hrsp'], float)
hepc_py = np.asarray(hp['hrpc', ][0] if False else hp.get('hrpc', hp.get('hepc')), float)
hrsp_ref = ref['hrsp']
hepc_ref = ref['hrpc']
print('python-run hrsp:', hrsp_py.shape, ' reference:', hrsp_ref.shape)
times_py = np.asarray(hp['times'], float)
print('tf grid match:', times_py.shape == ref['times'].shape,
      'max |t diff|', np.abs(times_py - ref['times']).max() if times_py.shape == ref['times'].shape else 'NA')
# channel alignment
labels_py = [c['labels'] if isinstance(c, dict) else c['labels']
             for c in np.asarray(PY['chanlocs']).ravel()]
lab_ref = [str(c) for c in np.asarray(ref['all_channels']).ravel()]
if labels_py != lab_ref:
    missing = set(lab_ref) - set(labels_py)
    print('channels missing in our output:', missing or 'none')
    order = [labels_py.index(l) for l in lab_ref if l in labels_py]
    lab_used = [l for l in lab_ref if l in labels_py]
    hrsp_py = hrsp_py[order]
    hepc_py = hepc_py[order]
    rmask = [i for i, l in enumerate(lab_ref) if l in labels_py]
    hrsp_ref, hepc_ref = hrsp_ref[rmask], hepc_ref[rmask]
r_b = np.corrcoef(hrsp_py.ravel(), hrsp_ref.ravel())[0, 1]
r_c = np.corrcoef(hepc_py.ravel(), hepc_ref.ravel())[0, 1]
d_b = np.abs(hrsp_py - hrsp_ref)
print(f'END-TO-END HRSP corr {r_b:.4f}  max|dB| {d_b.max():.3f}  median {np.median(d_b):.4f}')
print(f'END-TO-END HRPC corr {r_c:.4f}')
# ROI HEP of our saved epochs vs the reference epoch mean (both regressed?)
data = np.asarray(PY['data'], float)
roi = ['F1','F2','F3','F4','Fz','FC1','FC2','FC3','FC4','FCz',
       'C1','C2','C3','C4','Cz']
ridx = [i for i, l in enumerate(labels_py) if l in roi]
hep_py = data[ridx].mean(axis=0).mean(axis=1)
hep_ref = ref['hep_roi_from_epochs']
n = min(hep_py.size, hep_ref.size)
r_h = np.corrcoef(hep_py[:n], hep_ref[:n])[0, 1]
print(f'ROI HEP (epoch mean) corr {r_h:.4f} over first {n} points')