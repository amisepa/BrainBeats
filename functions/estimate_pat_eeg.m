% ESTIMATE_PAT_EEG - Pulse arrival time (PAT) of a PPG estimated from the
% cardiac field artifact (CFA) of the EEG, when no ECG was recorded.
%
% The QRS complex of each heartbeat reaches the scalp electrodes with no
% delay (cardiac field artifact), while the PPG pulse arrives ~200-450 ms
% later. Averaging the EEG around the PPG beats (median across beats, average
% reference, 5-30 Hz) therefore shows the QRS field before the pulse: the
% latency of the global field power (GFP) peak between -650 and -100 ms gives
% the PAT. On the sample dataset: 420 ms (all 291 beats) and 432 ms (first 60
% beats), for 424 ms measured with the ECG (the QRS field peaks 4 ms after the
% R-peak).
%
% Usage:
%   [pat, info] = estimate_pat_eeg(X, fs, ppgbeats)
%
% Inputs:
%   X        - EEG data (channels x samples), heart channels excluded
%   fs       - sampling rate (Hz)
%   ppgbeats - PPG beat sample indices (pulse onsets or peaks)
%
% Outputs:
%   pat  - estimated PAT (ms); NaN if the CFA peak is not clear enough
%   info - struct: .pat (ms), .ratio (GFP peak / median GFP of the window;
%          the estimate is kept if >= 2), .nBeats, .gfp, .times (ms)
%
% Copyright (C) - Cedric Cannard, 2026

function [pat, info] = estimate_pat_eeg(X, fs, ppgbeats)

minRatio = 2;          % GFP peak / median GFP needed to trust the peak
searchWin = [-650 -100];

X = double(X);
X = X - mean(X, 1);                                 % average reference
[b, a] = butter(2, [5 30]/(fs/2));
X = filtfilt(b, a, X')';

lags = round(-0.7*fs):round(0.2*fs);
t = lags / fs * 1000;
ppgbeats = round(ppgbeats(:));
ppgbeats = ppgbeats(ppgbeats + lags(1) >= 1 & ppgbeats + lags(end) <= size(X,2));
nB = numel(ppgbeats);
info = struct('pat', NaN, 'ratio', NaN, 'nBeats', nB, 'gfp', [], 'times', t);
if nB < 30
    pat = NaN;
    return
end

% Median across beats (robust to large artifacts), then GFP
E = zeros(size(X,1), numel(lags), nB);
for i = 1:nB
    E(:,:,i) = X(:, ppgbeats(i) + lags);
end
gfp = std(median(E, 3), 0, 1);

w = t >= searchWin(1) & t <= searchWin(2);
tw = t(w);
[pk, iMax] = max(gfp(w));
info.gfp = gfp;
info.ratio = pk / median(gfp);
if info.ratio >= minRatio
    pat = -tw(iMax);
else
    pat = NaN;
end
info.pat = pat;
