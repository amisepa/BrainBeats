"""csd_transform - Current source density (surface Laplacian) of EEG data,
with spherical splines (Perrin et al., 1989; Kayser & Tenke, 2006).

Python port of BrainBeats functions/csd_transform.m (+ GetGH.m,
current_source_density.m, read_chanlocs.m). Reference-free: channel positions
from the 10-05 template (standard_1005.ced) by label; returns the linear
channels x channels transform so the same C can be applied to other data of
the same montage (csd = C * data).

Defaults: mcont = 4, smoothl = 1e-5, headrad = 10 cm.

Copyright (C) Cedric Cannard, 2025 (MATLAB); Python port 2026-10.
"""
from __future__ import annotations

import os

import numpy as np
from scipy.special import eval_legendre

__all__ = ['read_ced_chanlocs', 'get_gh', 'csd_matrix', 'csd_transform_data']

_SPLINE_ORDER = list(range(1, 51))    # N = 50 Legendre iterations (GetGH.m:111)


# ---------------------------------------------------------------------------
# read_chanlocs.m, 'ced' case (24-45)
# ---------------------------------------------------------------------------
def read_ced_chanlocs(path):
    """labels, theta, phi from a .ced file (sph_theta + 90 wrap, sph_phi)."""
    rows = []
    with open(path, 'r', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            tok = line.split()
            if len(tok) >= 10:
                try:                      # first token numeric index
                    float(tok[0])
                except ValueError:
                    continue              # header line
                rows.append(tok[1:10])
    labels = [r[0] for r in rows]
    sph_theta = np.asarray([float(r[6]) for r in rows], float)
    sph_phi = np.asarray([float(r[7]) for r in rows], float)
    theta = sph_theta + 90.0
    theta = np.where(theta > 180.0, theta - 360.0, theta)
    phi = sph_phi
    return labels, theta, phi


# ---------------------------------------------------------------------------
# GetGH.m:75-131
# ---------------------------------------------------------------------------
def get_gh(theta, phi, m=4):
    """G, H matrices (nElec x nElec) for spherical spline flexibility m (2-10)."""
    theta = np.asarray(theta, float) % 360.0
    phi = np.asarray(phi, float)
    if not (2 <= int(m) <= 10):
        raise ValueError(f'Invalid spline flexibility m = {m} (use 2..10)')
    th, ph = np.deg2rad(theta), np.deg2rad(phi)
    x, y, z = (np.sin(ph) * np.cos(th), np.sin(ph) * np.sin(th), np.cos(ph))
    # MATLAB sph2cart(az, el, r): x = r cos(el) cos(az) -- az=theta, el=phi
    x = np.cos(ph) * np.cos(th)
    y = np.cos(ph) * np.sin(th)
    z = np.sin(ph)
    d2 = ((x[:, None] - x[None, :]) ** 2 + (y[:, None] - y[None, :]) ** 2 +
          (z[:, None] - z[None, :]) ** 2)
    ef = 1.0 - d2 / 2.0                       # cosine distances
    if ef.max() > 1.0:
        ef /= ef.max()
    if ef.min() < -1.0:
        ef /= abs(ef.min())
    mm = float(m)
    g = np.zeros_like(ef)
    h = np.zeros_like(ef)
    for n in _SPLINE_ORDER:
        p = eval_legendre(n, ef)
        nn = float(n)
        g += (2.0 * nn + 1.0) * p / (nn * nn + nn) ** mm
        h += (-(2.0 * nn + 1.0)) * p / (nn * nn + nn) ** (mm - 1.0)
    return g / (4.0 * np.pi), -h / (4.0 * np.pi)


# ---------------------------------------------------------------------------
# current_source_density.m:45-88 (vectorized already)
# ---------------------------------------------------------------------------
def csd_matrix(n_elec, G, H, smoothl=1e-5, headrad=10.0):
    """The linear CSD transform C (channels x channels): csd = C * data."""
    Gl = G + smoothl * np.eye(n_elec)
    Gi = np.linalg.inv(Gl)
    tc = Gi.sum(axis=1)                       # nElec x 1
    sgi = tc.sum()
    # C applied to a (nElec x nPnts) data matrix Z (current_source_density.m
    # with the per-sample electrode centering folded in exactly):
    #   C * Z = H * (Gi*Z - TC * (TC'/sgi * Z)) / head^2        (head^2 = headScale)
    #     = (H * (Gi - TC*TC'/sgi)) * Z / head^2
    c = H @ (Gi - np.outer(tc, tc) / sgi) / (headrad ** 2)
    return c


# ---------------------------------------------------------------------------
# csd_transform.m:42-58 wrapper
# ---------------------------------------------------------------------------
def csd_transform_data(data, labels, chanlocfile=None,
                       mcont=4, smoothl=1e-5, headrad=10.0):
    """Apply the CSD transform to data (nChan x ... , any trailing shape).

    labels: the data's channel labels, all found in the 10-05 template.
    Returns (csd_data, C). EEG channels only (the template holds no ECG/PPG).
    """
    if chanlocfile is None:
        chanlocfile = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                   'standard_1005.ced')
    if not os.path.isfile(chanlocfile):
        raise FileNotFoundError(f'Channel location file not found: {chanlocfile}')
    tpl_labels, tpl_theta, tpl_phi = read_ced_chanlocs(chanlocfile)
    tpl = {str(l).lower(): i for i, l in enumerate(tpl_labels)}
    idx = np.empty(len(labels), dtype=int)
    missing = []
    for i, lab in enumerate(labels):
        j = tpl.get(str(lab).lower())
        if j is None:
            missing.append(lab)
        else:
            idx[i] = j
    if missing:
        raise ValueError('csd_transform: channel labels not in the 10-05 '
                         f'template: {", ".join(map(str, missing))}')
    G, H = get_gh(tpl_theta[idx], tpl_phi[idx], mcont)
    C = csd_matrix(len(labels), G, H, smoothl, headrad)
    d = np.asarray(data, float)
    shape = d.shape
    out = C @ d.reshape(d.shape[0], -1)
    return out.reshape(shape), C