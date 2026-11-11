function test_thresholds
% TEST_THRESHOLDS - Thresholds that must not depend on the data's units
% (absolute ones failed in volts before v1.7). No EEGLAB needed.
%   1. DATA_RANK: an average-referenced, then component-removed dataset has
%      the expected rank in uV, in V and in single precision.
%   2. GET_RR flatline check: a synthetic ECG passes with an offset and in
%      volts; a constant (with an offset) and zeros are refused.
%   3. GET_SQI_PPG: the same beat classes for a synthetic PPG in any unit.
%
% Copyright (C) - Cedric Cannard, 2026

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'functions'));
rng(1);

% 1. rank
X = randn(32, 5000) * 20;                           % uV
X = X - mean(X, 1);                                 % average reference: rank 31
[U, ~, ~] = svd(X, 'econ');
X = X - U(:,1) * (U(:,1)' * X);                     % one component removed: 30
r = [data_rank(X) data_rank(1e-6*X) data_rank(single(X)) data_rank(single(1e-6*X))];
assert(all(r == 30), 'data_rank: %s (expected 30)', mat2str(r))
fprintf('1 OK: rank 30 in uV, V, single\n');

% 2. ECG flatline
fs = 250; n = 120*fs; t = (0:n-1)/fs;
ecg = zeros(1, n);
for b = 1:round(0.8*fs):n-fs
    ecg = ecg + 1000*exp(-((1:n) - b - fs/2).^2 / (2*(0.01*fs)^2));   % QRS-like, 1 mV
end
ecg = ecg + 10*randn(1, n);
p = struct('fs', fs, 'heart_signal', 'ecg');
for s = {ecg + 5e5, 1e-6*ecg}
    [~, ~, pk] = evalc_rr(s{1}, t*1000, p);
    assert(abs(numel(pk) - 149) <= 2, 'get_RR: %d beats (expected ~149)', numel(pk))
end
for s = {5e5 + zeros(1, n), zeros(1, n)}
    failed = false;
    try
        evalc_rr(s{1}, t*1000, p);
    catch ME
        failed = contains(ME.message, 'flat line');
    end
    assert(failed, 'get_RR: a flat line was not refused')
end
fprintf('2 OK: ECG with an offset and in volts detected; flat lines refused\n');

% 3. PPG classes in any unit
ppg = zeros(1, n); onsets = 1 + (fs:round(0.8*fs):n-fs);
for b = onsets
    k = b:min(n, b + round(0.8*fs) - 1); tt = (k - b) / fs;
    ppg(k) = exp(-(tt - 0.15).^2/(2*0.05^2)) + 0.4*exp(-(tt - 0.4).^2/(2*0.06^2));
end
ppg = 2000 * (ppg + 0.01*randn(1, n));
[~, ~, a1] = get_sqi_ppg(onsets, ppg(:), fs);
[~, ~, a2] = get_sqi_ppg(onsets, 1e-3*ppg(:), fs);
assert(isequal(a1, a2), 'get_sqi_ppg: beat classes change with the units')
fprintf('3 OK: PPG classes identical at x1 and x1e-3 (%.0f%% ''Q'')\n', 100*mean(strcmp(a1,'Q')));
fprintf('All threshold checks passed.\n')
end

function varargout = evalc_rr(x, t, p)
[~, varargout{1:3}] = evalc('get_RR(x, t, p)');
end
