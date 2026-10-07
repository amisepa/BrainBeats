"""run_hep - Epoch EEG around heartbeats for HEP / HRSP / HRPC analyses.

Python port of functions/run_HEP.m (v1.6). Steps and their MATLAB lines:

  1.  Epoch window: params['hep_window'] in ms, or 'adaptive' (run_HEP.m:87-106).
  2.  Reject beats followed by the next within epoch_end + 50 ms (qrsMargin),
      then IBI outliers (Grubbs kept-only, isoutlier).     (108-129)
  3.  Append 'R-peak' events (1-based latencies, .beat numbering), sort via
      eeg_checkset eventconsistency.                        (131-141)
  4.  Epoch at {'R-peak'} with 650 ms padding, 'epochinfo','yes'.  (176-180)
  5.  Cleaning (clean_eeg stage-1 on the epochs): bad epochs (find_badTrials,
      Grubbs default), ICA (if none: Picard at the effective rank; Infomax
      default), fitted on a 1-Hz high-passed copy when highpass < 1, ICLabel,
      pop_icflag thresholds [.99 .9 conf .99 .99] (heart conf = .75 default),
      pop_subcomp.                                          (182-212, clean_eeg)
  6.  Crop the padded epochs to the analysis window (pop_select 'point'). (214-216)
  7.  hepBeats = original-train beats of the surviving epochs (241-243);
      boundaries from 'boundary' events.                    (244-247)
  8.  Continuous data cleaned with the SAME linear cleaning, estimated from
      the epochs: T = Yc Xr' pinv(Xr Xr'), error check, Xc = T EEG.data. (249-267)
  9.  hep_level 'ics'/'both': hep_all_ics (A = W*Xc[icachansind], polarity
      signed by the IC map's max; ICLabel posteriors if present). (269-273, 433-462)
 10.  'ref','csd': surface Laplacian (apply_csd) - NOT ported yet.  (275-282)
 11.  hep_baseline 'regression': baseline_regression on the epochs, slopes +
      baselines stored.                                     (284-297)
 12.  ROI (default frontocentral 15, fallback Fz then Cz then chan 1). (299-309)
 13.  Channel measures: compute_hep_tf on the ROI mean -> brainbeats.roi, on
      ALL channels -> brainbeats.hrsp (rmfield 'hep'), surrogate ->
      brainbeats.surrogate (+ channels/times/freqs/hep.real).  (311-336)
 14.  keep_heart: heart channels re-epoched in their own units, appended.  (339-353)
 15.  Save <name>_HEP.set next to the input.               (424-428)

NOT ported: the plotting block (356-423; the MATLAB GUI owns the figures; the
Python package gets its own figure module later).

Divergences, named:
  * Surrogate RNG (numpy PCG64 vs MATLAB twister): same distribution, different
    draws for the same seed; deterministic outputs are faithful.
  * MATLAB's pop_epoch 'epochinfo' event bookkeeping is reproduced through
    HEP['epoch'] entries carrying {'event': [index], ...}; the eegprep
    pop_epoch returns (EEG, indices), used directly.
  * find_badTrials' 45-Hz low-pass port uses scipy firwin/filtfilt (the MATLAB
    design_fir/filtfilt_fast pair); the outlier test is the same Grubbs.

Copyright (C) Cedric Cannard, 2023 (MATLAB); Python port 2026-10.
"""

from __future__ import annotations

import os

import numpy as np

from .params import default_params
from .compute_hep_tf import compute_hep_tf
from .baseline_regression import baseline_regression
from .matlab_utils import matlab_round, percentile_matlab, grubbs_outliers

__all__ = ['run_hep']


def run_hep(EEG, CARDIO=None, params=None, Rpeaks=None):
    """HEP/HRSP/HRPC of a continuous EEG dataset time-locked to heartbeats.

    EEG    - continuous eegprep dataset (pop_loadset dict), EEG channels only
             (the caller removes the heart channels, like brainbeats_process)
    CARDIO - optional dataset holding the heart channel(s) (keep_heart)
    params - mapping (default_params()); missing keys get MATLAB defaults
    Rpeaks - 1-based heartbeat sample indices into EEG.data (clean_rr output)
    """
    import eegprep
    params = default_params(**(dict(params) if params else {}))
    params['analysis'] = 'hep'
    params['fs'] = EEG['srate']
    fs = float(EEG['srate'])

    # --- 1. epoch window (87-106) ------------------------------------------
    Rpeaks = np.asarray(Rpeaks, dtype=float).ravel()
    IBI = np.concatenate([np.diff(Rpeaks), [np.nan]]) / fs * 1000.0   # ms
    qrs_margin = 50.0
    if isinstance(params.get('hep_window'), str):
        if params['hep_window'] != 'adaptive':
            raise ValueError("'hep_window' must be [start end] in ms or 'adaptive'.")
        win_end = np.floor((percentile_matlab(IBI, 5) - qrs_margin) / 10.0) * 10.0
        win_end = float(min(max(win_end, 400), 1000))
        epoch_win = [-300.0, win_end]
        print(f"Adaptive HEP window: {epoch_win[0]:g} to {epoch_win[1]:g} ms "
              f"(5th percentile of the inter-beat intervals: "
              f"{percentile_matlab(IBI, 5):.0f} ms).")
    else:
        epoch_win = [float(v) for v in params['hep_window']]
    params['hep_window'] = epoch_win

    # --- 2. beat rejection (108-129) ----------------------------------------
    min_ibi = epoch_win[1] + qrs_margin
    keep = ~(IBI < min_ibi)
    if (~keep).any():
        print(f'Removing {int((~keep).sum())}/{IBI.size} heartbeats followed by '
              f'the next one within {min_ibi:g} ms (epoch end + {qrs_margin:g} ms).')
    kept = np.flatnonzero(keep & np.isfinite(IBI))
    outl = kept[grubbs_outliers(IBI[kept])]
    if outl.size:
        print(f'Removing {outl.size} outlier trials with the following '
              'interbeat intervals (IBI):')
        for o in outl:
            print(f'   - IBI = {IBI[o]:g} (ms)')
        keep[outl] = False
    print(f'{int(keep.sum())}/{Rpeaks.size} heartbeats kept for HEP epochs.')
    all_beats = Rpeaks.astype(np.int64)          # whole train, for surrogates
    Rpeaks = Rpeaks[keep].astype(np.int64)
    IBI = IBI[keep]

    # --- 3. events + sort (131-141) -----------------------------------------
    n_ev = len(EEG['event'])
    events = list(EEG['event'])
    for i, b in enumerate(Rpeaks):
        events.append({'type': 'R-peak', 'latency': float(b), 'duration': 0.0,
                       'urevent': n_ev + i + 1, 'beat': i + 1})
    EEG['event'] = events
    EEG = eegprep.eeg_checkset(EEG, 'eventconsistency')

    # --- 4. padded epochs (176-180) -------------------------------------------
    pad = 650.0
    out = eegprep.pop_epoch(EEG, ['R-peak'],
                            ((epoch_win[0] - pad) / 1000.0,
                             (epoch_win[1] + pad) / 1000.0))
    HEPwide = out[0] if isinstance(out, tuple) else out
    raw_wide = {k: v for k, v in HEPwide.items()}

    # --- 5. cleaning on the epochs (182-212; clean_eeg stage 1) ---------------
    removed_trials = np.empty(0, dtype=int)
    removed_components = np.empty(0, dtype=int)
    do_chans = params.get('hep_level', 'channels') in ('channels', 'both')
    do_ics = params.get('hep_level', 'channels') in ('ics', 'both')
    use_csd = bool(params.get('clean_eeg')) and params.get('ref') == 'csd'
    if params.get('clean_eeg'):
        HEPwide, removed_trials, removed_components = _clean_eeg_epochs(
            HEPwide, params, fs)
        HEPwide.setdefault('brainbeats', {}).setdefault(
            'preprocessings', {})
        HEPwide['brainbeats']['preprocessings']['removed_eeg_trials'] = removed_trials
        HEPwide['brainbeats']['preprocessings']['removed_eeg_components'] = removed_components
    elif do_ics and not (lambda w: np.asarray(w).size if w is not None else 0)(HEPwide.get('icaweights')):
        # no cleaning: ICA (Picard) + nothing removed, for the IC measures only
        src = HEPwide
        data = np.asarray(src['data'], float)
        rank = int((np.linalg.eigvalsh(
            np.cov(data.reshape(data.shape[0], -1))) > 1e-7).sum())
        kw = (dict(icatype='picard', maxiter=400, mode='standard') if icamethod == 1 else
               dict(icatype='runica', extended=1, pca=rank))
        ica = eegprep.pop_runica(dict(src), **kw)
        ica = ica[0] if isinstance(ica, tuple) else ica
        for f in ('icaweights', 'icasphere', 'icawinv', 'icachansind'):
            HEPwide[f] = ica[f]
        HEPwide = eegprep.eeg_checkset(HEPwide)

    # --- 6. crop (214-216) ------------------------------------------------------
    times_w = np.asarray(HEPwide['times']).ravel()
    i_win = np.flatnonzero((times_w >= epoch_win[0]) & (times_w < epoch_win[1]))
    HEP = _pop_select(HEPwide, 'point', [int(i_win[0]) + 1, int(i_win[-1]) + 1])
    HEP = eegprep.eeg_checkset(HEP)
    HEP.setdefault('brainbeats', {})['preprocessings'] = {
        'hep_window': epoch_win,
        'removed_eeg_trials': removed_trials,
        'removed_eeg_components': removed_components}

    # --- 7. beats of the surviving epochs + boundaries (241-247) --------------
    # HEP.epoch(e).event is a 1-based index (or list) into HEP.event; eegprep
    # event dicts carry 'beat' -- use it, else the event's urevent position.
    beat_nums = []
    for ep in HEP['epoch']:
        evs = ep.get('event') or []
        evs = list(evs) if isinstance(evs, (list, np.ndarray)) else [evs]
        if not evs:
            continue                       # epoch without a beat event
        ev_i = int(evs[0]) - 1 if not isinstance(evs[0], dict) else None
        if ev_i is None:
            beat_nums.append(int(evs[0].get('beat', 0)))
            continue
        ev = HEP['event'][ev_i]
        beat_nums.append(int(ev.get('beat', ev_i + 1 - n_ev)))
    if not beat_nums:
        raise RuntimeError('run_hep: no epochs carry a beat event.')
    beat_idx = np.asarray(beat_nums, dtype=int)
    hep_beats = Rpeaks[beat_idx - 1]
    HEP['brainbeats']['preprocessings']['hep_beats'] = hep_beats

    bnd = []
    for ev in EEG['event']:
        if str(ev.get('type', '')).lower() == 'boundary':
            bnd.append(float(ev['latency']))
    bnd = np.asarray(bnd, float) if bnd else np.empty(0)

    # --- 8. continuous data with the same linear cleaning (249-267) -----------
    Xc = np.asarray(EEG['data'], dtype=float)
    if (do_chans or do_ics) and params.get('clean_eeg'):
        # HEPwide may have lost epochs at the data edges (pop_epoch boundary
        # warnings), which are NOT in removed_trials: align by count from the
        # end (kept trials keep the original order).
        n_kept = int(np.asarray(HEPwide['data']).shape[2])
        kept_trials = np.setdiff1d(np.arange(1, raw_wide['trials'] + 1),
                                   np.asarray(removed_trials, dtype=int))
        if kept_trials.size > n_kept:
            print(f'WARNING: {kept_trials.size - n_kept} edge epoch(s) lost at '
                  'the recording boundaries; aligning the T estimate.')
            kept_trials = kept_trials[kept_trials.size - n_kept:]
        Xr = np.asarray(raw_wide['data'], float)[:, :, kept_trials - 1]
        Xr = Xr.reshape(Xr.shape[0], -1)
        Yc = np.asarray(HEPwide['data'], float)
        Yc = Yc.reshape(Yc.shape[0], -1)
        T = (Yc @ Xr.T) @ np.linalg.pinv(Xr @ Xr.T)
        err = np.linalg.norm(Yc - T @ Xr, 'fro') / np.linalg.norm(Yc, 'fro')
        if err > 1e-3:
            print(f'WARNING: the EEG cleaning of the epochs is not reproduced '
                  f'exactly on the continuous data (relative error {err:.1e}).')
        Xc = T @ Xc
        params['_clean_transform'] = T          # kept for apply_csd later

    # --- 9. IC measures (269-273) ------------------------------------------------
    if do_ics:
        HEP['brainbeats']['ics'] = hep_all_ics(
            HEPwide, Xc, hep_beats, epoch_win, fs,
            freqs=_tf_freqs(params), n_surr=params.get('hep_surrogates', 0),
            surr_mode=params.get('hep_surrogate_mode', 'shuffle'),
            all_beats=all_beats, boundaries=bnd,
            hep_times=np.asarray(HEP['times']).ravel(),
            seed=params.get('surrogate_seed', 1))

    # --- 10. CSD (275-282) --------------------------------------------------------
    if use_csd:
        raise NotImplementedError("ref='csd' not ported yet (apply_csd.m); "
                                  "use ref='average'.")

    # --- 11. baseline regression (284-297) --------------------------------------
    if str(params.get('hep_baseline', 'none')).lower() == 'regression':
        bw = params.get('hep_baseline_win') or (-150, -50)
        bw = [float(v) for v in bw]
        print(f'Regression-based baseline correction ({bw[0]:g} to {bw[1]:g} ms)...')
        data, beta, bl = baseline_regression(np.asarray(HEP['data'], float),
                                             np.asarray(HEP['times']).ravel(), bw)
        HEP['data'] = data
        HEP['brainbeats']['preprocessings']['baseline_regression'] = {
            'window': bw,
            'channels': [c['labels'] if isinstance(c, dict) else c['labels']
                         for c in np.asarray(HEP['chanlocs']).ravel()],
            'beta': beta, 'baseline': bl}

    # --- 12. ROI (299-309) ---------------------------------------------------------
    labels = [c['labels'] if isinstance(c, dict) else c['labels']
              for c in np.asarray(HEP['chanlocs']).ravel()]
    roi_labels = params.get('hep_roi') or         ['F1','F2','F3','F4','Fz','FC1','FC2','FC3','FC4','FCz',
         'C1','C2','C3','C4','Cz']
    low = [str(l).lower() for l in labels]
    roi = [i for i, l in enumerate(low) if l in [str(r).lower() for r in roi_labels]]
    if not roi:
        for cand in ('fz', 'cz'):
            roi = [i for i, l in enumerate(low) if l == cand]
            if roi:
                break
        if not roi:
            roi = [0]
        print(f'WARNING: none of the ROI channels is in the data: using '
              f'{labels[roi[0]]}.')
    HEP['brainbeats']['preprocessings']['hep_roi'] = [labels[i] for i in roi]

    # --- 13. channel measures (311-336) -------------------------------------------
    tf_opts = dict(freqs=_tf_freqs(params),
                   n_surr=params.get('hep_surrogates', 0),
                   surr_mode=params.get('hep_surrogate_mode', 'shuffle'),
                   all_beats=all_beats, boundaries=bnd,
                   hep_times=np.asarray(HEP['times']).ravel(),
                   seed=params.get('surrogate_seed', 1))
    if do_chans:
        sig = Xc[np.asarray(roi, int)].mean(axis=0)
        n_surr = params.get('hep_surrogates', 0)
        if n_surr and n_surr > 0:
            print(f'Surrogate heartbeat control ({n_surr:g} surrogates, '
                  f'{params.get("hep_surrogate_mode", "shuffle")})...')
        tf_main, surr_main = compute_hep_tf(
            sig, fs, hep_beats, epoch_win, cycles=params.get('hep_tf_cycles', 5),
            tstep=params.get('hep_tf_tstep', 10), **tf_opts)
        roi_struct = {'channels': [labels[i] for i in roi],
                      'hep': tf_main['hep'],
                      'hep_times': tf_main['hep_times'],
                      'tf': {k: v for k, v in tf_main.items() if k != 'hep'},
                      'surrogate': surr_main}
        HEP['brainbeats']['roi'] = roi_struct

        print(f'Computing HRSP and HRPC on {Xc.shape[0]} channels...')
        tf, surr = compute_hep_tf(Xc, fs, hep_beats, epoch_win,
                                  cycles=params.get('hep_tf_cycles', 5),
                                  tstep=params.get('hep_tf_tstep', 10), **tf_opts)
        tf['channels'] = labels
        HEP['brainbeats']['hrsp'] = {k: v for k, v in tf.items() if k != 'hep'}
        if surr is not None:
            surr['channels'] = tf['channels']
            surr['times'] = tf['times']
            surr['hep_times'] = tf['hep_times']
            surr['freqs'] = tf['freqs']
            surr.setdefault('hep', {})['real'] = tf['hep']
            HEP['brainbeats']['surrogate'] = surr
            pf = surr['hep']['p_fdr']
            print(f'Surrogate control: {int((pf < .05).sum())}/{pf.size} HEP '
                  'points (channels x latencies) differ from the surrogates '
                  '(FDR-corrected p < .05).')

    # --- 14. keep the heart channel (339-353) --------------------------------------
    n_eeg = np.asarray(HEP['data']).shape[0]
    if params.get('keep_heart') and CARDIO is not None:
        hep_idx = np.rint(np.asarray(HEP['times']).ravel() / 1000.0 * fs).astype(int)
        cdata = np.asarray(CARDIO['data'], float)
        Hd = np.stack([cdata[:, (hep_beats[b] - 1) + hep_idx]
                       for b in range(HEP['trials'])], axis=-1)
        HEP['data'] = np.concatenate([np.asarray(HEP['data']), Hd], axis=0)
        HEP['nbchan'] = HEP['data'].shape[0]
        chans = list(np.asarray(HEP['chanlocs']).ravel())
        chans += [{'labels': ch} for ch in params['heart_channels']]
        HEP['chanlocs'] = np.array(chans, dtype=object)
        HEP = eegprep.eeg_checkset(HEP)

    # --- 15. save (424-428) ----------------------------------------------------------
    if params.get('save'):
        fdir = EEG.get('filepath') or '.'
        fname = EEG.get('filename') or 'dataset.set'
        newname = os.path.join(fdir, str(fname)[:-4] + '_HEP.set')
        eegprep.pop_saveset(HEP, newname)
        # eegprep 0.3.0 pop_saveset drops custom (non-EEGLAB) fields such as
        # EEG.brainbeats; re-add them to the saved header (MAT v5) so the
        # output is a faithful BrainBeats .set.
        try:
            import scipy.io as sio
            with open(newname, 'rb') as fh:
                hdr = sio.loadmat(fh, struct_as_record=False, squeeze_me=True)
            # eegprep writes EEG fields at the top level of the .set but
            # DROPS custom fields like .brainbeats; re-attach them via a v5
            # roundtrip and keep a sidecar copy in case the roundtrip fails.
            import scipy.io as sio

            def _safe(o):
                if o is None:
                    return np.zeros(0)
                if isinstance(o, dict):
                    return {str(k): _safe(v) for k, v in o.items()}
                if isinstance(o, (list, tuple)):
                    return np.asarray([np.zeros(0) if v is None else v
                                       for v in o], dtype=object)
                return o

            bb = _safe(HEP.get('brainbeats'))
            sidecar = newname[:-4] + '_brainbeats.mat'
            sio.savemat(sidecar, {'brainbeats': bb}, oned_as='row')
            with open(newname, 'rb') as fh:
                hdr = sio.loadmat(fh, struct_as_record=True, squeeze_me=False)
            hdr = {k: v for k, v in hdr.items() if not k.startswith('__')}
            hdr['brainbeats'] = bb
            with open(newname, 'wb') as fh:
                sio.savemat(fh, hdr, oned_as='row')
        except Exception as exc:
            print(f'WARNING: could not re-attach the brainbeats struct to the '
                  f'saved file ({exc}).')
    return HEP


# ---------------------------------------------------------------------------
def hep_all_ics(HEPwide, Xc, hep_beats, epoch_win, fs, *, freqs, n_surr,
                surr_mode, all_beats, boundaries, hep_times, seed=1):
    """HEP/HRSP/HRPC + surrogates of every IC (run_HEP.m:433-462)."""
    import eegprep
    icaweights = HEPwide.get('icaweights')
    if icaweights is None or np.asarray(icaweights).size == 0:
        print("WARNING: No ICA decomposition: the IC measures need one.")
        return {}
    icachansind = np.asarray(HEPwide['icachansind']).ravel().astype(int)
    W = np.asarray(icaweights, float) @ np.asarray(HEPwide['icasphere'], float)
    maps = np.asarray(HEPwide['icawinv'], float)
    probs = None
    try:
        probs = HEPwide['etc']['ic_classification']['ICLabel']['classifications']
    except (KeyError, TypeError):
        probs = None
    if probs is None or np.asarray(probs).shape[0] != W.shape[0]:
        try:
            tmp = eegprep.pop_iclabel(dict(HEPwide), 'default')
            tmp = tmp[0] if isinstance(tmp, tuple) else tmp
            probs = tmp['etc']['ic_classification']['ICLabel']['classifications']
        except Exception:
            probs = None
    i_max = np.argmax(np.abs(maps), axis=0)
    sgn = np.sign(maps[i_max, np.arange(maps.shape[1])])
    A = (W @ Xc[icachansind - 1]) * sgn[None, :]
    print(f'Computing HEP, HRSP and HRPC of {A.shape[0]} independent components...')
    tf, surr = compute_hep_tf(A, fs, hep_beats, epoch_win, freqs=freqs,
                              n_surr=n_surr, surr_mode=surr_mode,
                              all_beats=all_beats, boundaries=boundaries,
                              hep_times=hep_times, seed=seed)
    tf_r = {k: v for k, v in tf.items() if k != 'hep'}
    return {'ic': np.arange(1, A.shape[0] + 1)[None, :].T,
            'maps': maps * sgn[None, :], 'polarity': sgn,
            'chanlocs': np.asarray(HEPwide['chanlocs'])[icachansind - 1],
            'iclabel': probs,
            'iclabel_classes': ['Brain', 'Muscle', 'Eye', 'Heart', 'Line Noise',
                                'Channel Noise', 'Other'],
            'hep': tf['hep'], 'hep_times': tf['hep_times'], 'tf': tf_r,
            'surrogate': surr}


def _tf_freqs(params):
    f = params.get('hep_tf_freqs')
    if f is None:
        return None      # compute_hep_tf default 4:30
    return np.arange(float(f[0]), float(f[-1]) + 0.5, 1.0)


# ---------------------------------------------------------------------------
# eegprep wrappers
# ---------------------------------------------------------------------------
def _pop_select(EEG, mode, value):
    import eegprep
    out = eegprep.pop_select(EEG, mode, value)
    return out[0] if isinstance(out, tuple) else out


def _clean_eeg_epochs(HEPwide, params, fs):
    """clean_eeg stage-1 on the padded epochs (HEP route)."""
    import eegprep
    hp = params.get('highpass', 0.5)
    icamethod = params.get('icamethod', 2)

    # (a) bad epochs (find_badTrials.m)
    bad = _find_bad_trials(HEPwide, params.get('detectMethod', 'grubbs'))
    if bad.size:
        out = eegprep.pop_rejepoch(dict(HEPwide), bad + 1, 0)
        HEPwide = out[0] if isinstance(out, tuple) else out
    print(f'find_badTrials: {bad.size} bad epochs removed.')

    # (b) ICA if none (run_HEP 184-191 + clean_eeg 360-393)
    has_ica = (lambda w: np.asarray(w).size if w is not None else 0)(HEPwide.get('icaweights'))
    if not has_ica:
        src = HEPwide
        if hp < 1:
            out = eegprep.pop_eegfiltnew(dict(HEPwide), locutoff=1.0,
                                         minphase=params.get('filttype') == 'causal')
            src = out[0] if isinstance(out, tuple) else out
            print('ICA fitted on a 1-Hz high-passed copy of the data.')
        data = np.asarray(src['data'], float)
        rank = int((np.linalg.eigvalsh(
            np.cov(data.reshape(data.shape[0], -1))) > 1e-7).sum())
        if icamethod == 1:
            kw = (dict(icatype='picard', maxiter=400, mode='standard') if icamethod == 1 else
               dict(icatype='runica', extended=1, pca=rank))
        elif icamethod == 2:
            kw = dict(icatype='runica', extended=1, pca=rank)
        else:
            kw = dict(icatype='runica', extended=1, pca=rank,
                      lrate=1e-5, maxsteps=2000)
        ica = eegprep.pop_runica(dict(src), **kw)
        ica = ica[0] if isinstance(ica, tuple) else ica
        for f in ('icaweights', 'icasphere', 'icawinv', 'icachansind'):
            HEPwide[f] = ica[f]
        HEPwide = eegprep.eeg_checkset(HEPwide)

    # (c) ICLabel + flag + remove (clean_eeg 404-451)
    out = eegprep.pop_iclabel(dict(HEPwide), 'default')
    HEPwide = out[0] if isinstance(out, tuple) else out
    conf = params.get('conf_thresh', 0.75)
    if params.get('rm_heart_ics') is False:
        conf = np.nan
    thr = np.array([[np.nan, np.nan], [0.99, 1.], [0.9, 1.], [conf, 1.],
                    [0.99, 1.], [0.99, 1.], [np.nan, np.nan]])
    out = eegprep.pop_icflag(dict(HEPwide), thr)
    HEPwide = out[0] if isinstance(out, tuple) else out
    bad_comp = np.flatnonzero(np.asarray(HEPwide['reject']['gcompreject']))
    if bad_comp.size:
        print(f'Removing {bad_comp.size} bad component(s).')
        out = eegprep.pop_subcomp(dict(HEPwide), (bad_comp + 1).tolist(), 0)
        HEPwide = out[0] if isinstance(out, tuple) else out
    return HEPwide, bad, bad_comp


def _find_bad_trials(EEG, method='grubbs'):
    """Port of find_badTrials.m: amplitude (RMS over channels of per-channel
    RMS) + high-frequency noise (RMS of the MAD of the data minus its 45-Hz
    low-passed version) outliers."""
    from scipy.signal import firwin, filtfilt
    from .matlab_utils import grubbs_outliers
    data = np.asarray(EEG['data'], float)
    n_trials = data.shape[2]
    # design_fir(100, [2*[0 45 50]/fs, 1], [1 1 0 0]) = LOW-PASS 0-45 Hz, order 100
    b = firwin(101, 45.0 / (EEG['srate'] / 2.0), pass_zero=True)
    sig_amp = np.empty(n_trials)
    sig_snr = np.empty(n_trials)
    for i in range(n_trials):
        ep = data[:, :, i]
        sig_amp[i] = np.sqrt((np.sqrt((ep ** 2).mean(axis=1)) ** 2).mean())
        lp = filtfilt(b, 1.0, ep, axis=1, padlen=3 * (b.size - 1))
        sig_snr[i] = np.sqrt((np.abs(ep - lp).mean(axis=1) ** 2).mean())
    bad = grubbs_outliers(sig_amp, 0.05) | grubbs_outliers(sig_snr, 0.05)
    return np.flatnonzero(bad)
