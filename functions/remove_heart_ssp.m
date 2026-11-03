% REMOVE_HEART_SSP - Remove the cardiac field artifact from continuous EEG by
% signal-space projection (SSP), as MNE's compute_proj_ecg: the first n
% spatial components of the heartbeat-locked average EEG (-200 to 400 ms
% around the R-peaks) are projected out of the data.
%
% The projection removes those scalp topographies at all times, so any
% heartbeat-locked brain activity with a similar topography is removed too.
%
% Usage:
%   [EEG, info] = remove_heart_ssp(EEG, beats, params)
%
% Inputs:
%   EEG    - continuous EEG (EEG channels only, filtered and re-referenced)
%   beats  - R-peak sample indices
%   params - BrainBeats parameters; params.ssp_n: number of components
%            (default 2, as MNE's n_eeg)
%
% Outputs:
%   EEG  - EEG with the components projected out (rank lowered by n)
%   info - struct: n_components, patterns (channels x n), explained (% of
%          the heartbeat-locked average), beats_used, and the heart-locked
%          EEG amplitude before and after (see CFA_REPORT)
%
% Copyright (C) - Cedric Cannard, 2026

function [EEG, info] = remove_heart_ssp(EEG, beats, params)

n = 2;
if isfield(params,'ssp_n') && ~isempty(params.ssp_n), n = params.ssp_n; end
fs = EEG.srate;
X = double(EEG.data);
w = round(-0.2*fs):round(0.4*fs);
beats = beats(beats + w(1) >= 1 & beats + w(end) <= size(X,2));

% Beats whose window holds no large EEG artifact (GFP under median + 6
% robust SD, as for the ECG regression)
g = std(X - mean(X,1), 0, 1);
ok = g < median(g) + 6 * 1.4826 * mad(g, 1);
beats = beats(arrayfun(@(b) all(ok(b + w)), beats));
if numel(beats) < 30
    error('BrainBeats:input', 'SSP: fewer than 30 heartbeats without large artifacts.')
end

% Heartbeat-locked average, then its first n spatial components
avg = zeros(size(X,1), numel(w));
for i = 1:numel(beats)
    avg = avg + X(:, beats(i) + w);
end
avg = avg / numel(beats);
avg = avg - mean(avg(:, w/fs < -0.1), 2);
[U, Sv] = svd(avg, 'econ');
sv = diag(Sv);
U = U(:, 1:n);
Y = X - U * (U' * X);

info = cfa_report(cfa_amplitude(X, beats, fs), cfa_amplitude(Y, beats, fs), sprintf('SSP, %d components', n));
info.n_components = n;
info.patterns = U;
info.explained = 100 * sv(1:n).^2 / sum(sv.^2);
info.beats_used = numel(beats);
EEG.data = cast(Y, 'like', EEG.data);
