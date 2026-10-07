% REMOVE_HEART_REGRESSION - Remove the cardiac field artifact from continuous
% EEG by regressing the ECG out of each EEG channel.
%
% The cardiac field artifact is the volume-conducted field of the heart: at
% each electrode, a scaled copy of the cardiac signal. The ECG channel(s),
% filtered like the EEG and shifted by -20 to +20 ms (one copy per sample),
% are fitted to each EEG channel by least squares over the samples without
% large EEG artifacts, and the fitted part is subtracted from all samples. The few short lags absorb the
% differences of waveform between the ECG lead and the field at the scalp.
% No ICA is needed, and the heart signal must be an ECG.
%
% The lags are fixed: longer lags remove more artifact but also heartbeat-
% locked brain responses, which are partly predictable from the ECG. On the
% sample dataset, the heart-locked EEG amplitude from -50 to 100 ms goes
% from 2.11 to 0.59 uV (-72%), with 0.8% of the EEG variance removed. In a
% test on the same data (291 heartbeats), 91% of a synthetic 1-uV response
% added 300 ms after each R-peak was kept with +/-20 ms (99% with no lag,
% for -44% of artifact; 59% with +/-100 ms, for -76%).
%
% Usage:
%   [EEG, info] = remove_heart_regression(EEG, CARDIO, beats, params)
%
% Inputs:
%   EEG    - continuous EEGLAB dataset (EEG channels only, filtered)
%   CARDIO - EEGLAB dataset with the ECG channel(s), same sampling rate and
%            samples as EEG (raw: it is filtered here like the EEG)
%   beats  - R-peak sample indices, for the amplitude report
%   params - BrainBeats parameters (reads highpass, lowpass, filttype,
%            vis_cleaning)
%
% Outputs:
%   EEG  - EEG with the ECG-related part removed
%   info - struct: lags_ms, samples_fitted (%), cfa_before, cfa_after (heart-locked EEG amplitude
%          from -50 to 100 ms, uV; see CFA_AMPLITUDE), variance_removed (%)
%
% Copyright (C) - Cedric Cannard, 2026

function [EEG, info] = remove_heart_regression(EEG, CARDIO, beats, params)

maxLag = 20;   % ms

% ECG filtered like the EEG (same band and filter type, so same delays)
hp = 0.5; lp = 30; causal = false;
if isfield(params,'highpass') && ~isempty(params.highpass), hp = params.highpass; end
if isfield(params,'lowpass') && ~isempty(params.lowpass), lp = params.lowpass; end
if isfield(params,'filttype'), causal = strcmpi(params.filttype,'causal'); end
ECG = pop_eegfiltnew(CARDIO, 'locutoff', hp, 'minphase', causal);
ECG = pop_eegfiltnew(ECG, 'hicutoff', lp, 'minphase', causal);

n = min(EEG.pnts, ECG.pnts);
X = double(EEG.data(:,1:n));
ecg = double(ECG.data(:,1:n));

% Lagged copies of the ECG channel(s) as regressors
L = round(maxLag/1000 * EEG.srate);
R = zeros(size(ecg,1)*(2*L+1), n);
k = 0;
for iChan = 1:size(ecg,1)
    for lag = -L:L
        k = k + 1;
        R(k,:) = circshift(ecg(iChan,:), lag);
    end
end
R = R - mean(R, 2);

% Least-squares fit on the samples free of large EEG artifacts (global field
% power within 6 robust SD of its median): a few seconds of large artifacts
% otherwise dominate the fit (on the sample data, the electrode artifacts of
% the first seconds made the heart-locked amplitude 20 times larger instead
% of smaller). The fitted part is then subtracted from all samples.
gfp = std(X, 0, 1);
ok = gfp < median(gfp) + 6 * 1.4826 * mad(gfp, 1);
fprintf('ECG regression fitted on %.1f%% of the samples (the others hold large EEG artifacts). \n', 100*mean(ok));
B = (X(:,ok) * R(:,ok)') * pinv(R(:,ok) * R(:,ok)');
Y = X - B * R;

% Report on the artifact-free samples too: heartbeats whose window (-200 to
% 400 ms) holds no large EEG artifact
w = round(-0.2*EEG.srate):round(0.4*EEG.srate);
beats = beats(beats + w(1) >= 1 & beats + w(end) <= n);
beats = beats(arrayfun(@(b) all(ok(b + w)), beats));
info.lags_ms = (-L:L) / EEG.srate * 1000;
info.samples_fitted = 100 * mean(ok);   % %
info.cfa_before = cfa_amplitude(X, beats, EEG.srate);
info.cfa_after = cfa_amplitude(Y, beats, EEG.srate);
info.variance_removed = 100 * (1 - sum(sum(Y(:,ok).^2)) / sum(sum(X(:,ok).^2)));
fprintf('ECG regression (lags -%g to %g ms): heart-locked EEG amplitude (-50 to 100 ms) %.2f uV before, %.2f uV after (%.0f%% reduction); %.1f%% of the EEG variance removed. \n', ...
    maxLag, maxLag, info.cfa_before, info.cfa_after, 100*(1 - info.cfa_after/info.cfa_before), info.variance_removed);

oriEEG = EEG;
EEG.data(:,1:n) = Y;
if isfield(params,'vis_cleaning') && params.vis_cleaning
    vis_artifacts(EEG, oriEEG, 'ShowSetname', false);
    set(gcf,'Toolbar','none','Menu','none','Name','EEG before (red) and after (blue) ECG regression','NumberTitle','Off')
    finish_figure(gcf)
end
