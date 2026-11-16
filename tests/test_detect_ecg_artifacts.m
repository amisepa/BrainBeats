function test_detect_ecg_artifacts
% TEST_DETECT_ECG_ARTIFACTS - detect_ecg_artifacts on a synthetic ECG with
% known artifacts. No EEGLAB needed.
%   1. 'window' flags a muscle burst, a slow movement excursion and two
%      bursts 2.4 s apart (merged into one segment by MinGapSec: a 1-s
%      window already bridges gaps under ~2 s), and little else.
%   2. 'moving' with AmpThresh flags the excursion only at its extremes,
%      the reason for 'window'.
%   3. 'window' gives the same mask in uV and in V.
%
% Copyright (C) - Cedric Cannard, 2026

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'functions'));
rng(3);

fs = 250; n = 600*fs; t = (0:n-1)/fs;
ecg = zeros(1, n);
ibi = round(fs * (0.85 + 0.05*randn(1, 800)));      % ~70 bpm
b = round(fs/2);
k = 1;
while b < n - fs
    ecg = ecg + 1000*exp(-((1:n) - b).^2 / (2*(0.008*fs)^2)) ...   % R wave, 1 mV
              + 200*exp(-((1:n) - b - round(0.25*fs)).^2 / (2*(0.04*fs)^2));  % T wave
    b = b + ibi(k); k = k + 1;
end
ecg = ecg + 15*randn(1, n) + 100*sin(2*pi*0.25*t);   % noise, respiration drift

% artifacts
[bh, ah] = butter(4, [35 90]/(fs/2));
emg = filtfilt(bh, ah, randn(1, n));
emg = emg / std(emg) * 300;
burst = false(1, n); burst(round(120*fs):round(122*fs)) = true;          % 2-s burst
pair  = false(1, n); pair(round(400*fs):round(400.8*fs)) = true;         % two bursts
pair(round(403.2*fs):round(404.0*fs)) = true;                            % 2.4 s apart
ecg(burst | pair) = ecg(burst | pair) + emg(burst | pair);
exc = round(300*fs):round(305*fs);                                       % 5-s excursion
ecg(exc) = ecg(exc) + 5000*sin(pi*(0:numel(exc)-1)/(numel(exc)-1)).^2;

truth = burst | pair; truth(exc) = true;

% 1. window
[m, seg] = evalc_detect(ecg, fs, 'Method','window');
hit = @(a, b) mean(m(round(a*fs):round(b*fs)));
assert(hit(120.2, 121.8) == 1, 'window: burst not flagged')
assert(hit(300.5, 304.5) == 1, 'window: excursion not flagged whole (%.2f)', hit(300.5, 304.5))
assert(hit(400.2, 403.8) == 1, 'window: the two close bursts are not one segment')
inpair = seg(:,1) <= 400*fs + 50 & seg(:,2) >= 404*fs - 50;
assert(sum(inpair) == 1, 'window: close bursts not merged')
[~, segB] = evalc_detect(ecg, fs, 'Method','window', 'MinGapSec', 0.1);
assert(~any(segB(:,1) <= 400*fs + 50 & segB(:,2) >= 404*fs - 50), ...
    'window: the pair is one segment even without merging (merge check is void)')
pad = conv(double(truth), ones(1, 4*fs), 'same') > 0;    % 2 s around each artifact
fa = mean(m(~pad));
assert(fa < 0.01, 'window: %.2f%% flagged away from the artifacts', 100*fa)
fprintf('1 OK: window flags burst, excursion, merged pair; %.2f%% elsewhere\n', 100*fa);

% 2. moving + AmpThresh: the excursion is flagged only in part
m2 = evalc_detect(ecg, fs, 'Method','moving', 'AmpThresh', 6);
cov2 = mean(m2(round(300.5*fs):round(304.5*fs)));
assert(cov2 < 1, 'moving: excursion flagged whole (test no longer shows the difference)')
fprintf('2 OK: moving flags %.0f%% of the excursion, window 100%%\n', 100*cov2);

% 3. units
m3 = evalc_detect(1e-6*ecg, fs, 'Method','window');
assert(isequal(m, m3), 'window: uV and V masks differ')
fprintf('3 OK: same mask in uV and V\n');
end

function [m, seg] = evalc_detect(x, fs, varargin)
[~, m, seg] = evalc('detect_ecg_artifacts(x, ''SampleRate'', fs, varargin{:})');
end
