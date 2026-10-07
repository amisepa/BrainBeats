"""compute_hep_tf - HEP, HRSP, HRPC and the surrogate heartbeat control.

Python port of BrainBeats functions/compute_hep_tf.m (v1.6, in order):

  1. Output grids (MATLAB 102-109):
       tIdx   = round(win1/1000*fs) : k : round(win2/1000*fs)-1   (k ~ tstep ms)
       hepIdx = round(win1/1000*fs) : round(win2/1000*fs)-1       (every sample)
     The END SAMPLE IS EXCLUSIVE (MATLAB ":-1"); hep_times (e.g. HEP.times)
     overrides the HEP grid when given.
  2. Wavelet support: half = ceil(3*max_sigma*fs); beats are kept when
     beat + hepIdx(1) - half >= 1 and beat + hepIdx(end) + half <= nPts and no
     boundary event lies inside the window (valid_beats, MATLAB 236-244).
     Fewer than 10 kept beats raises.
  3. Surrogate trains (130-163):
       'rigid': shifts are uniform integers over shift ms (sample-rounded),
                |shift| >= min(0.25*median_IBI, half the max |shift range|);
                each surrogate = beats + shift, valid_beats-filtered.
       'shuffle': start = allB(1) + floor(rand*ibi(1)); the WHOLE IBI sequence
                (of allBeats) is permuted and cumsum'd; trimmed to nB valid
                beats (sorted).
  4. HEP: mean of the epochs of the same beats (mean_epochs, 246-253).
  5. HRSP: zero-phase Morlet coefficients per (channel, frequency):
       w = 2*exp(-(tk/fs)^2/(2 sig^2))*exp(1i*2pi f tk/fs)/sum|w|,
       tk = -ceil(3 sig fs) : ceil(3 sig fs) samples, sig = cycles/(2 pi f);
       nfft = 2^nextpow2(nPts + 2*half + 1); co = ifft(fft(X).*W);
       power = mean over beats of |co(beats+tIdx)|^2  (tf_stats, 255-263)
       hrpc  = (|sum(z)|^2 - N)/(N(N-1)), z = co/|co| (Vinck 2010, unbiased)
       hrsp  = 10*log10(power / mean-over-beats power)  (dB, 223)
  6. Surrogate stats (265-282): z vs the null mean/sd, p two-sided (HEP, HRSP)
     or one-sided (HRPC), empirical p (+1)/(nSurr+1), Benjamini-Hochberg FDR
     over all points, null_lo/hi at 2.5/97.5 percentiles.

NOTE vs heart_core_v25.morlet_cwt (the HEP_neurofeedback kernel): L1-normalised
without the factor 2 and a slightly different support (floor+1 vs ceil samples
per side). Power scale differs by x4 and the dB measure cancels it; HRSP/
HRPC/HEP are all unaffected. The BrainBeats-native kernel is used here for
parity with compute_hep_tf.m, including raw tf.power (uV^2).

Divergence, named: the surrogate RNG (numpy PCG64 vs MATLAB twister) gives
different surrogate draws for the same seed. Real outputs are deterministic.

Copyright (C) Cedric Cannard, 2026 (MATLAB); Python port 2026-10.
"""
from __future__ import annotations

import numpy as np

from .matlab_utils import matlab_round, fdr_bh

__all__ = ['compute_hep_tf', 'surrogate_trains', 'morlet_kernels_bb',
           'surrogate_stats', 'valid_beats']


# ---------------------------------------------------------------------------
# kernels (compute_hep_tf.m:188-199 -- BrainBeats' own construction)
# ---------------------------------------------------------------------------
def morlet_kernels_bb(freqs, fs, cycles=5.0, nfft=None):
    """Zero-phase Morlet kernels (x2 amplitude-normalised), one per frequency.

    Returns (W, half): W is (nF x nfft) complex FFTs, half the SHARED wavelet
    support in samples (ceil(3*sigma_max*fs)). With nfft=None, time-domain
    kernels are returned instead (list of complex arrays).
    """
    freqs = np.atleast_1d(np.asarray(freqs, float))
    sig = cycles / (2 * np.pi * freqs)              # temporal SD, seconds
    half = int(np.ceil(3 * sig.max() * fs))
    kernels = []
    for f, s in zip(freqs, sig):
        tk = np.arange(-int(np.ceil(3 * s * fs)), int(np.ceil(3 * s * fs)) + 1)
        w = np.exp(-(tk / fs) ** 2 / (2 * s ** 2)) * \
            np.exp(1j * 2 * np.pi * f * tk / fs)
        w = 2 * w / np.sum(np.abs(w))               # sinusoid amp A -> |coef| = A
        kernels.append(w)
    if nfft is None:
        return kernels, half
    W = np.zeros((freqs.size, nfft), dtype=complex)
    for i, w in enumerate(kernels):
        wp = np.zeros(nfft, dtype=complex)
        wp[:w.size] = w
        wp = np.roll(wp, -(w.size - 1) // 2)        # kernel centre -> sample 0
        W[i] = np.fft.fft(wp)
    return W, half


def _nextpow2(n):
    return int(2 ** np.ceil(np.log2(n)))


# ---------------------------------------------------------------------------
# beat validity (valid_beats, 236-244)
# ---------------------------------------------------------------------------
def valid_beats(b, lims, half, n_pts, bnd):
    """b, lims, bnd are MATLAB 1-based, n_pts the number of samples."""
    b = np.asarray(b, dtype=np.int64)
    a = b + lims[0] - half          # 1-based first sample of the window
    z = b + lims[1] + half          # 1-based last sample
    v = (a >= 1) & (z <= n_pts)
    if bnd is not None and len(np.atleast_1d(bnd)):
        for x in np.asarray(bnd, float).ravel():
            v &= ~((x > a) & (x < z))
    return v


# ---------------------------------------------------------------------------
# mean epochs (mean_epochs, 246-253)
# ---------------------------------------------------------------------------
def _mean_epochs(X, beats, idx):
    """X 0-based (nChan x nPts); beats/idx 1-based. Beating NaN handling:
    MATLAB indexes blindly (`X(:, beats(i)+idx)`) and errors outside; the port
    zero-pads nothing -- callers pre-filter with valid_beats, so the windows
    always fit. Boundary crossing is pre-filtered too."""
    X = np.atleast_2d(np.asarray(X, float))
    pos = (np.asarray(beats, dtype=np.int64)[:, None] +
           np.asarray(idx, dtype=np.int64)[None, :]) - 1      # 0-based
    if pos.min() < 0 or pos.max() >= X.shape[1]:
        raise IndexError('beat window outside the data (valid_beats slipped)')
    return X[:, pos].mean(axis=1)               # nChan x nTimes (mean over beats)


# ---------------------------------------------------------------------------
# tf_stats (255-263)
# ---------------------------------------------------------------------------
def _tf_stats(co, beats, t_idx, n_pts):
    """co: (n_samples,) complex; power and PPC over beats at tIdx."""
    pos = (np.asarray(beats, dtype=np.int64)[:, None] +
           np.asarray(t_idx, dtype=np.int64)[None, :]) - 1    # 0-based
    ok = (pos >= 0) & (pos < n_pts)
    C = co[np.clip(pos, 0, n_pts - 1)]
    C = np.where(ok, C, np.nan + 0j)
    N = int(np.asarray(beats).size)
    with np.errstate(invalid='ignore', divide='ignore'):
        pw = np.nanmean(np.abs(C) ** 2, axis=0)
        z = np.where(np.isfinite(C.real), C / np.abs(C), np.nan + 0j)
        s = np.nansum(np.nan_to_num(z), axis=0)
        pc = (np.abs(s) ** 2 - N) / (N * (N - 1)) if N > 1 else np.full(t_idx.size, np.nan)
    return pw, pc


# ---------------------------------------------------------------------------
# surrogate stats (265-282)
# ---------------------------------------------------------------------------
def surrogate_stats(real, null, tail='both'):
    """null (..., nSurr); real broadcastable against null[..., 0]."""
    null = np.asarray(null, float)
    real = np.asarray(real, float)
    n_surr = null.shape[-1]
    null_mean = null.mean(axis=-1)
    null_sd = null.std(axis=-1, ddof=1)
    with np.errstate(invalid='ignore', divide='ignore'):
        z = (real - null_mean) / null_sd
        dev = np.abs(null - null_mean[..., None])
    from scipy.special import erfc
    if tail == 'both':
        p = erfc(np.abs(z) / np.sqrt(2))
        obs = np.abs(real - null_mean)[..., None]
        p_emp = (np.sum(dev >= obs, axis=-1) + 1) / (n_surr + 1)
    else:
        p = 0.5 * erfc(z / np.sqrt(2))
        p_emp = (np.sum(null >= real[..., None], axis=-1) + 1) / (n_surr + 1)
    return {'null_mean': null_mean, 'null_sd': null_sd,
            'null_lo': np.percentile(null, 2.5, axis=-1),
            'null_hi': np.percentile(null, 97.5, axis=-1),
            'z': z, 'p': p, 'p_emp': p_emp, 'p_fdr': fdr_bh(p)}


# ---------------------------------------------------------------------------
# surrogate trains (130-163)
# ---------------------------------------------------------------------------
def surrogate_trains(beats, *, n_surr, mode, fs, all_beats=None,
                     lims=None, half=0, n_pts=None, bnd=None,
                     shift=(-500.0, 500.0), seed=1, n_target=None):
    """beats: 1-based valid real beats. Returns (trains, shifts):
    a list of 1-based int arrays (already valid_beats-filtered) and, for
    'rigid', the shifts in samples (None for 'shuffle')."""
    rng = np.random.default_rng(seed)
    beats = np.asarray(beats, dtype=np.int64)
    out = []
    if str(mode).lower() == 'rigid':
        lo = int(matlab_round(shift[0] / 1000.0 * fs))
        hi = int(matlab_round(shift[1] / 1000.0 * fs))
        min_abs = int(matlab_round(0.25 * np.median(np.diff(beats))))
        min_abs = min(min_abs, int(np.floor(0.5 * max(abs(lo), abs(hi)))))
        shifts = np.empty(int(n_surr), dtype=np.int64)
        for s in range(int(n_surr)):
            d = 0
            while abs(d) < min_abs:
                d = lo + int(np.floor(rng.random() * (hi - lo + 1)))
            shifts[s] = d
        for s in range(int(n_surr)):
            b = beats + shifts[s]
            out.append(b[valid_beats(b, lims, half, n_pts, bnd)])
        return out, shifts
    # 'shuffle' (default)
    all_b = np.sort(np.round(np.asarray(
        all_beats if all_beats is not None and len(np.atleast_1d(all_beats))
        else beats, dtype=float))).astype(np.int64)
    ibi = np.diff(all_b)
    if n_target is None:
        n_target = beats.size
    for s in range(int(n_surr)):
        start = all_b[0] + int(np.floor(rng.random() * ibi[0]))
        perm = rng.permutation(ibi.size)
        b = np.concatenate([[start], start + np.cumsum(ibi[perm])])
        b = np.asarray(b, dtype=np.int64)
        b = b[valid_beats(b, lims, half, n_pts, bnd)]
        if b.size > n_target:
            keep = rng.permutation(b.size)[:n_target]
            b = np.sort(b[keep])
        out.append(b)
    return out, None


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def compute_hep_tf(X, fs, beats, win, *, freqs=None, cycles=5.0, tstep=10.0,
                   tf=True, n_surr=0, surr_mode='shuffle', all_beats=None,
                   shift=(-500.0, 500.0), boundaries=None, seed=1,
                   hep_times=None, keep_null=False):
    """[tf, surr] = compute_hep_tf(X, fs, beats, win, opts).

    X     (nChan x nPts) continuous (already-cleaned) data
    beats 1-based heartbeat sample indices
    win   (start, end) ms
    hep_times - the epochs' own time points in ms (overrides the HEP grid)
    """
    X = np.asarray(X, dtype=float)
    if X.ndim == 1:
        X = X[None, :]
    n_chan, n_pts = X.shape
    if freqs is None:
        freqs = np.arange(4, 31, dtype=float)          # default 4:30
    else:
        freqs = np.atleast_1d(np.asarray(freqs, float))
        if freqs.size == 2:                            # keep MATLAB 4:30 semantics
            pass
    n_f = freqs.size
    win = (float(win[0]), float(win[1]))
    beats_in = np.asarray(beats, np.int64).ravel()

    # --- 1. output grids (102-109) ------------------------------------------
    k = max(1, int(matlab_round(tstep / 1000.0 * fs)))
    i1 = int(matlab_round(win[0] / 1000.0 * fs))
    i2 = int(matlab_round(win[1] / 1000.0 * fs)) - 1   # end EXCLUSIVE (MATLAB :-1)
    t_idx = np.arange(i1, i2 + 1, k)
    times = t_idx / fs * 1000.0
    if hep_times is not None and len(np.atleast_1d(hep_times)):
        hep_idx = np.rint(np.asarray(hep_times, float).ravel() /
                          1000.0 * fs).astype(int)
    else:
        hep_idx = np.arange(i1, i2 + 1)
    lims = (int(hep_idx[0]), int(hep_idx[-1]))

    # --- 2. beat validity ----------------------------------------------------
    sig = cycles / (2 * np.pi * freqs)
    half = int(np.ceil(3 * sig.max() * fs)) if tf else 0
    bnd = (np.asarray(boundaries, float).ravel()
           if boundaries is not None and len(np.atleast_1d(boundaries))
           else np.empty(0))
    v = valid_beats(beats_in, lims, half, n_pts, bnd)
    if (~v).any():
        print(f'compute_hep_tf: {int((~v).sum())}/{beats_in.size} heartbeats '
              'too close to the edges or to a discontinuity were left out.')
    beats = beats_in[v]
    n_b = int(beats.size)
    if n_b < 10:
        raise ValueError('compute_hep_tf: fewer than 10 heartbeats.')

    # --- 3. surrogate trains ---------------------------------------------------
    shifts = None
    surr_beats = []
    if n_surr and n_surr > 0:
        surr_beats, shifts = surrogate_trains(
            beats, n_surr=n_surr, mode=surr_mode, fs=fs, all_beats=all_beats,
            lims=lims, half=half, n_pts=n_pts, bnd=bnd, shift=shift, seed=seed)

    # --- 4. HEP ---------------------------------------------------------------
    tf_out = {'hep': _mean_epochs(X, beats, hep_idx),
              'hep_times': hep_idx / fs * 1000.0}
    surr = None
    if n_surr and n_surr > 0:
        surr = {'nSurr': int(n_surr), 'mode': str(surr_mode),
                'shifts': None if shifts is None else shifts / fs * 1000.0}
        null_hep = np.stack([_mean_epochs(X, b, hep_idx) for b in surr_beats],
                            axis=-1)
        surr['hep'] = surrogate_stats(tf_out['hep'], null_hep, 'both')
        if keep_null:
            surr['hep']['null'] = null_hep

    tf_out.update(times=times, freqs=freqs, cycles=cycles, nBeats=n_b)
    if not tf:
        return tf_out, surr

    # --- 5. HRSP + HRPC --------------------------------------------------------
    nfft = _nextpow2(n_pts + 2 * half + 1)
    W, _ = morlet_kernels_bb(freqs, fs, cycles, nfft=nfft)
    power = np.full((n_chan, n_f, t_idx.size), np.nan)
    hrpc = np.full_like(power, np.nan)
    null_p = np.zeros((n_chan, n_f, t_idx.size, int(n_surr))) if n_surr else None
    null_c = np.zeros_like(null_p) if n_surr else None
    for i_ch in range(n_chan):
        fx = np.fft.fft(X[i_ch], nfft)
        for i_f in range(n_f):
            co = np.fft.ifft(fx * W[i_f])[:n_pts]
            pw, pc = _tf_stats(co, beats, t_idx + 1, n_pts)   # 1-based tIdx
            power[i_ch, i_f] = pw
            hrpc[i_ch, i_f] = pc
            if n_surr:
                for s, sb in enumerate(surr_beats):
                    pw_s, pc_s = _tf_stats(co, sb, t_idx + 1, n_pts)
                    null_p[i_ch, i_f, :, s] = pw_s
                    null_c[i_ch, i_f, :, s] = pc_s
    tf_out['power'] = power
    with np.errstate(divide='ignore', invalid='ignore'):
        tf_out['hrsp'] = 10.0 * np.log10(power / power.mean(axis=2, keepdims=True))
    tf_out['hrpc'] = hrpc

    # --- 6. surrogate stats for HRSP/HRPC --------------------------------------
    if n_surr:
        null_p = null_p.astype(float)
        with np.errstate(divide='ignore', invalid='ignore'):
            null_h = 10.0 * np.log10(null_p / null_p.mean(axis=2, keepdims=True))
        surr['hrsp'] = surrogate_stats(tf_out['hrsp'], null_h, 'both')
        surr['hrpc'] = surrogate_stats(tf_out['hrpc'], null_c.astype(float), 'right')
    return tf_out, surr
