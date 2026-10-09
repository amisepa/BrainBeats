"""Parity stage B: HEP/HRSP/HRPC from the stage-A cleaned continuous data.

R-peak train: re-detected with the validated get_rr (Pan-Tompkins port) on the
RAW ECG, then clean_rr -- the exact beat source the MATLAB run used. Only the
267 beats present in dataset_HEP.set should matter for the TF measures; beats
the python pipeline keeps in addition are harmless (compute_hep_tf re-filters
by valid_beats), but a COUNT mismatch here points at a real beat-source
divergence, so it is reported.
"""
import os
os.environ.setdefault('MPLBACKEND', 'Agg')
import sys
import numpy as np
sys.path.insert(0, r'C:/Users/ccann/Documents/MATLAB/BrainBeats/python')
sys.path.insert(0, r'C:/Users/ccann/Documents/HEP_neurofeedback/python')
import eegprep
from functions.compute_hep_tf import compute_hep_tf
from get_rr_v25 import get_rr
from rr_cleaning_v25 import clean_rr, qrs_bandpass

BB = r'C:/Users/ccann/Documents/MATLAB/BrainBeats'
STAGE_A = os.path.join(BB, 'tests', 'parity_data.npz')

d = np.load(STAGE_A, allow_pickle=True)
Xc = d['data']                       # cleaned continuous EEG (62 ch, CAR+interp)
lat = d['lat']                       # 267 reference 1-based latencies
fs = 250.0

# --- beat train (the MATLAB path: get_RR on the raw ECG, then clean_rr) ----
EEGraw = eegprep.pop_loadset(os.path.join(BB, 'sample_data', 'dataset.set'))
labels = [c['labels'] if isinstance(c, dict) else c['labels']
          for c in np.asarray(EEGraw['chanlocs']).ravel()]
i_ecg = labels.index('ECG')
ecg = np.asarray(EEGraw['data'])[i_ecg].astype(float)
peaks_full, info = get_rr(ecg, fs, drop_first=False)   # full train incl. first beat
rr = np.diff(peaks_full) / fs
peaks_oneb = peaks_full[1:]                    # get_RR drops the first beat (its contract)
peak_amp = qrs_bandpass(ecg, fs)[peaks_oneb - 1]     # Step 2 reconstruction (NEVER inert)
nn, peaks_clean, idx_bad, info_rr = clean_rr(rr, peaks_oneb, fs, peak_amp=peak_amp,
                                ecg_signal=ecg, interpolate_missing=False)
peaks_clean = np.asarray(peaks_clean, float)
finite = np.isfinite(peaks_clean)
peaks_clean = peaks_clean[finite].astype(np.int64)     # 1-based indices
print(f'beat train: {peaks_full.size} raw -> {int(np.isfinite(peaks_clean).sum())} after clean_rr.')
print('reference kept beats:', lat.size)
# The reference run kept beats subject to the same window/IBI filters; compare
# overlap: ref latencies (1-based) vs our cleaned train
lset = np.asarray([int(round(float(l))) for l in lat])
our = peaks_clean
overlap = np.isin(lset, our)
print(f'reference latencies found in the python train: {overlap.sum()}/{lset.size}')
extra = int((~np.isin(our, lset)).sum())
print(f'python train beats not in the reference epochs (removed there by '
      f'cleaning/epoching): {extra}')

# --- measures on the ROI mean and all channels ----------------------------
ref = np.load(os.path.join(BB, 'python', 'tests', 'data', 'reference_hep_tf.npz'),
              allow_pickle=True)
roi_labels = list(ref['roi_channels'])
all_labels = list(ref['all_channels'])
lab_now = [c['labels'] if isinstance(c, dict) else c['labels']
           for c in np.asarray(EEGraw['chanlocs']).ravel()]
lab_now = [l for l in lab_now if l not in ('PPG', 'ECG')]
ridx = [i for i, l in enumerate(lab_now) if l in roi_labels]
print('roi chans in cleaned data:', [lab_now[i] for i in ridx])
sig_roi = Xc[ridx].mean(axis=0)
hep_beats_1b = lset                      # the SAME beats as the MATLAB epochs
win = (-300.0, 600.0)
hep_times = np.asarray(ref['hep_times'])

tf_roi, surr_roi = compute_hep_tf(sig_roi[None, :], fs, hep_beats_1b, win,
                                  hep_times=hep_times, n_surr=0)
tf_all, surr_all = compute_hep_tf(Xc, fs, hep_beats_1b, win,
                                  hep_times=hep_times, n_surr=0)

hrsp_py = tf_all['hrsp']
hrpc_py = tf_all['hrpc']
power_py = tf_all['power']
hrsp_ref = ref['hrsp']
hrpc_ref = ref['hrpc']
power_ref = ref['power']
print('shapes py/ref:', hrsp_py.shape, hrsp_ref.shape)

# channel order: the reference channels are the cleaned epochs' channel order
# (TP9 interpolated at its original position); our cleaned data dropped and
# re-inserted TP9 via pop_interp, so the order should match all_labels.
assert lab_now == list(all_labels), 'channel order mismatch'

r = np.corrcoef(hrsp_py.ravel(), hrsp_ref.ravel())[0, 1]
r_c = np.corrcoef(hrpc_py.ravel(), hrpc_ref.ravel())[0, 1]
r_p = np.corrcoef(power_py.ravel(), power_ref.ravel())[0, 1]
d_b = np.abs(hrsp_py - hrsp_ref)
d_c = np.abs(hrpc_py - hrpc_ref)
print(f'HRSP: corr {r:.6f}, max|dB diff| {d_b.max():.4f}, median {np.median(d_b):.4f}')
print(f'HRPC: corr {r_c:.6f}, max|diff| {d_c.max():.4f}, median {np.median(d_c):.4f}')
print(f'power: corr {r_p:.6f}')
# HEP ROI: compare against the epoch-mean reference
hep_roi_ref = ref['hep_roi_from_epochs']
hep_py = tf_roi['hep'][0]
r_h = np.corrcoef(hep_py, hep_roi_ref)[0, 1]
print(f'HEP ROI: corr {r_h:.6f}, max|diff| {np.abs(hep_py - hep_roi_ref).max():.2e} uV')

out = dict(r_hrsp=r, r_hrpc=r_c, r_power=r_p, r_hep=r_h,
           max_db=float(d_b.max()), med_db=float(np.median(d_b)),
           max_pc=float(d_c.max()), max_hep=float(np.abs(hep_py - hep_roi_ref).max()))
np.savez(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'parity_result.npz'), **out)
print('PARITY RESULT:', out)