% ESTIMATE_PAT - Pulse arrival time (PAT): delay between each ECG R-peak and
% the following PPG pulse fiducial (pulse-wave onset or peak).
%
% A PPG pulse reaches the sensor after the electrical heartbeat: PAT is the
% pre-ejection period plus the pulse transit time to the sensor (typically
% ~200-300 ms to the finger for the pulse onset). To time-lock HEPs to
% heartbeats from a PPG, the PPG beats must be shifted back by the PAT
% (BRAINBEATS_PROCESS option 'ppg_transit').
%
% Usage:
%   [patMed, pat, info] = estimate_pat(rpeaks, ppgpeaks, fs)
%   [patMed, pat, info] = estimate_pat(rpeaks, ppgpeaks, fs, patRange)
%
% Inputs:
%   rpeaks   - ECG R-peak sample indices
%   ppgpeaks - PPG fiducial sample indices (same sampling rate and time axis)
%   fs       - sampling rate (Hz)
%   patRange - plausible PAT range in ms (default [50 600]); each PPG beat
%              is paired with the closest preceding R-peak within this range
%
% Outputs:
%   patMed - median PAT (ms)
%   pat    - PAT of each paired PPG beat (ms, NaN if no R-peak in range)
%   info   - struct: .median, .iqr, .sd (ms), .nPaired, .nPPG
%
% Copyright (C) - Cedric Cannard, 2026

function [patMed, pat, info] = estimate_pat(rpeaks, ppgpeaks, fs, patRange)

if nargin < 4 || isempty(patRange), patRange = [50 600]; end
rpeaks = sort(rpeaks(:));
ppgpeaks = ppgpeaks(:);

pat = nan(size(ppgpeaks));
for i = 1:numel(ppgpeaks)
    iR = find(rpeaks < ppgpeaks(i), 1, 'last');   % closest preceding R-peak
    if ~isempty(iR)
        d = (ppgpeaks(i) - rpeaks(iR)) / fs * 1000;
        if d >= patRange(1) && d <= patRange(2)
            pat(i) = d;
        end
    end
end

patMed = median(pat, 'omitnan');
info.median = patMed;
info.iqr = iqr(pat(~isnan(pat)));
info.sd = std(pat, 'omitnan');
info.nPaired = sum(~isnan(pat));
info.nPPG = numel(ppgpeaks);
if info.nPaired < 0.5*info.nPPG
    warning('estimate_pat: only %g of %g PPG beats could be paired with an R-peak: check that both signals are aligned.', ...
        info.nPaired, info.nPPG)
end
