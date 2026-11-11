% BB_FUZZY_FRACTAL - Fuzzy entropy and fractal dimension of each row of a
% signal matrix, computed by the ASCENT plugin (version 2.1 or later).
%
% Fuzzy entropy: ASCENT's COMPUTE_FUZZEN with its defaults, 'local' mode (the
% mean of each embedding vector is removed, as in the original FuzzyEn of
% Chen et al., 2007), exponential kernel, each signal z-scored first.
% Fractal dimension: ASCENT's COMPUTE_FRACDIM (box counting, Theil-Sen fit,
% burst rejection).
%
% Up to BrainBeats 1.6, fuzzy entropy was computed without removing the
% vector means: values are higher now (median +0.18 on the NN series of 66
% resting recordings, rank correlation with the old values 0.86; +0.16 on
% 48 channels of 6 cleaned EEG recordings, rank correlation 0.91), so they
% must not be pooled with results of older versions. A warning says so at
% each call (warning ID 'BrainBeats:fuzzyEntropy'), with a dialog once per
% MATLAB session. Fractal dimension hardly changed (median difference 0 on
% HRV, +0.006 on EEG, at most 0.055).
%
% Usage:
%   [fe, fd, vers] = bb_fuzzy_fractal(signals, m, r, n, tau, parallel)
%
% Inputs:
%   signals  - one signal per row (channels x samples)
%   m, r, n, tau - fuzzy entropy embedding dimension, tolerance (fraction of
%              the SD), fuzzy power and time lag
%   parallel - parfor over the rows on the open pool (default false)
%
% Outputs:
%   fe, fd   - fuzzy entropy and fractal dimension (rows x 1), rounded to 3
%              decimals; NaN when ASCENT 2.1 or later is not installed
%   vers     - ASCENT version used ('' if unknown or not installed)
%
% Copyright (C) - Cedric Cannard, 2026

function [fe, fd, vers] = bb_fuzzy_fractal(signals, m, r, n, tau, parallel)

persistent dialogShown
if nargin < 6, parallel = false; end
nRows = size(signals,1);
fe = nan(nRows,1); fd = nan(nRows,1);

vers = '';
if exist('eegplugin_ascent','file')
    try vers = eegplugin_ascent; catch, end
end
v = sscanf(vers, '%d.%d')';
if ~exist('compute_FuzzEn','file') || ~exist('compute_FracDim','file') || (numel(v) == 2 && v * [100; 1] < 201)
    msg = sprintf(['Fuzzy entropy and fractal dimension need the ASCENT plugin, version 2.1 or later ' ...
        '(installed: %s). Install or update it in EEGLAB (File > Manage extensions) or from ' ...
        'https://github.com/amisepa/Ascent. They are set to NaN.'], ifelse(isempty(vers), 'none', vers));
    warning('BrainBeats:ascent', '%s', msg);
    bb_warndlg(msg);
    return
end

msg = ['Fuzzy entropy is now computed by the ASCENT plugin as the original FuzzyEn (each embedding ' ...
    'vector minus its mean): values are higher than with BrainBeats 1.6 and earlier, so do not ' ...
    'pool them with results of older versions.'];
warning('BrainBeats:fuzzyEntropy', '%s Turn this message off with: warning(''off'',''BrainBeats:fuzzyEntropy'')', msg);
if isempty(dialogShown)
    bb_warndlg(msg);
    dialogShown = true;
end

% ASCENT's defaults given explicitly, so that a later change of its defaults
% does not change BrainBeats' values
opts = {'Parallel', logical(parallel), 'Progress', nRows > 1};
fe = round(compute_FuzzEn(double(signals), 'm', m, 'r', r, 'n', n, 'tau', tau, ...
    'Mode', 'local', 'Kernel', 'exponential', 'Normalize', true, opts{:}), 3);
fd = round(compute_FracDim(double(signals), 'RejectBursts', true, 'RobustFit', 'theilsen', ...
    'ScaleTrimIQR', true, opts{:}), 3);
fe = fe(:); fd = fd(:);
end

function o = ifelse(c, a, b)
if c, o = a; else, o = b; end
end
