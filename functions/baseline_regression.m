% BASELINE_REGRESSION - Regression-based baseline correction of epoched data
% (Alday, 2019).
%
% Instead of subtracting each trial's baseline (which assumes a slope of 1
% between baseline and post-stimulus activity), the baseline mean is used as a
% trial-level regressor. At each channel and time point, the slope between the
% trial baselines and the data is estimated by least squares, and each trial is
% corrected by that slope times its centered baseline:
%
%   corrected(t) = data(t) - beta(t) * (baseline - mean(baseline))
%
% The corrected epochs can then be analyzed like any other epochs. Because the
% baseline is centered, the average across all trials is unchanged: what is
% removed is the trial-to-trial variance explained by the baseline, and, when
% groups are given, the baseline differences between them (as in an ANCOVA).
% For heartbeat-evoked potentials this matters because the pre-R-peak window
% contains activity from the previous cardiac cycle, so a plain subtraction can
% move pre-R-peak differences into the post-R-peak window. Steinfath et al.
% (2026) recommend avoiding baseline correction for HEPs, or using this
% regression form, which does not introduce baseline differences.
%
% Usage:
%   [data, beta, bl] = baseline_regression(data, times, blWin)
%   [data, beta, bl] = baseline_regression(data, times, blWin, groups)
%
% Inputs:
%   data   - epoched data (channels x times x trials), e.g. EEG.data
%   times  - epoch time vector (ms), e.g. EEG.times
%   blWin  - baseline window [start end] (ms), e.g. [-300 -100]
%   groups - (optional) condition label for each trial (numeric or cell of
%            char). When the corrected epochs will be used to compare
%            conditions, give them here: the slope is then estimated within
%            conditions, so condition effects are not absorbed by the baseline.
%
% Outputs:
%   data   - baseline-corrected epochs (same size as input)
%   beta   - baseline slope for each channel and time point (channels x times)
%   bl     - baseline mean of each trial (channels x trials). The correction
%            can be undone with: data + beta .* permute(bl - mean(bl,2), [1 3 2])
%
% Example:
%   [EEG.data, beta] = baseline_regression(EEG.data, EEG.times, [-300 -100]);
%
% Reference:
%   Alday, P. M. (2019). How much baseline correction do we need in ERP
%   research? Extended GLM model can replace baseline correction while lifting
%   its limits. Psychophysiology, 56(12), e13451.
%   Steinfath et al. (2026). Heartbeat-evoked responses in M/EEG: A
%   systematic review of methods with suggestions. Psychophysiology, 63(4), e70297.
%
% Copyright (C) - Cedric Cannard, 2026

function [data, beta, bl] = baseline_regression(data, times, blWin, groups)

[nChan, nTimes, nTrials] = size(data);
blIdx = times >= blWin(1) & times <= blWin(2);
if ~any(blIdx)
    error('baseline_regression: no time point in the baseline window [%g %g] ms.', blWin)
end
if nTrials < 3
    error('baseline_regression: at least 3 trials are needed.')
end

% Design: one intercept per group (or a single intercept) + centered baseline
if nargin < 4 || isempty(groups)
    G = ones(nTrials,1);
else
    [~, ~, g] = unique(groups(:));
    G = double(g == 1:max(g));
end

bl = squeeze(mean(data(:,blIdx,:), 2));      % channels x trials
if nChan == 1, bl = bl(:)'; end
blc = bl - mean(bl, 2);                       % centered on the mean over all trials

beta = nan(nChan, nTimes);
for iChan = 1:nChan
    Y = squeeze(data(iChan,:,:))';            % trials x times
    if nTimes == 1, Y = Y(:); end
    X = [G blc(iChan,:)'];
    B = X \ Y;                                % least squares, all time points at once
    beta(iChan,:) = B(end,:);
    data(iChan,:,:) = permute(Y - blc(iChan,:)' * B(end,:), [3 2 1]);
end
