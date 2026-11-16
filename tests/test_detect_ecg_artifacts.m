function test_detect_ecg_artifacts
% TEST_DETECT_ECG_ARTIFACTS - detect_ecg_artifacts on a synthetic ECG with
% known artifacts. No EEGLAB needed.
%   1. 'window' flags a muscle burst, a movement (2-s, 1.5-Hz swing) and two
%      bursts 2.4 s apart (merged into one segment by MinGapSec: a 1-s
%      window already bridges gaps under ~2 s), and little else.
%   2. 'moving' with AmpThresh flags the movement in part at best (here not at all),
%      the reason for 'window'.
%   3. 'window' gives the same mask in uV and in V.
%   4. 'window' does not flag a slow baseline drift under clean R waves
%      (its 1-Hz high-pass), and does without the high-pass.
%   5. 'window' does not flag a 2-min stretch of 3x taller heartbeats (the
%      movement rule's local reference), and does against the recording's
%      median range.
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
beats = ecg;                                         % heartbeats alone (test 5)
ecg = ecg + 15*randn(1, n) + 100*sin(2*pi*0.25*t);   % noise, respiration drift

% artifacts
[bh, ah] = butter(4, [35 90]/(fs/2));
emg = filtfilt(bh, ah, randn(1, n));
emg = emg / std(emg) * 300;
burst = false(1, n); burst(round(120*fs):round(122*fs)) = true;          % 2-s burst
pair  = false(1, n); pair(round(400*fs):round(400.8*fs)) = true;         % two bursts
pair(round(403.2*fs):round(404.0*fs)) = true;                            % 2.4 s apart
ecg(burst | pair) = ecg(burst | pair) + emg(burst | pair);
exc = round(300*fs):round(302*fs);                                       % 2-s movement
ecg(exc) = ecg(exc) + 2000*sin(2*pi*1.5*(0:numel(exc)-1)/fs);            % 1.5-Hz swing

truth = burst | pair; truth(exc) = true;
drift = round(460*fs):round(480*fs);                                     % 20-s baseline drift
ecg(drift) = ecg(drift) + 8000*sin(pi*(0:numel(drift)-1)/(numel(drift)-1)).^2;

% 1. window
[m, seg] = evalc_detect(ecg, fs, 'Method','window');
hit = @(a, b) mean(m(round(a*fs):round(b*fs)));
assert(hit(120.2, 121.8) == 1, 'window: burst not flagged')
assert(hit(300.2, 301.8) == 1, 'window: excursion not flagged whole (%.2f)', hit(300.2, 301.8))
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

% 2. moving + AmpThresh: the movement is flagged in part at best
m2 = evalc_detect(ecg, fs, 'Method','moving', 'AmpThresh', 6);
cov2 = mean(m2(round(300.2*fs):round(301.8*fs)));
assert(cov2 < 1, 'moving: excursion flagged whole (test no longer shows the difference)')
fprintf('2 OK: moving flags %.0f%% of the excursion, window 100%%\n', 100*cov2);

% 3. units
m3 = evalc_detect(1e-6*ecg, fs, 'Method','window');
assert(isequal(m, m3), 'window: uV and V masks differ')
fprintf('3 OK: same mask in uV and V\n');

% 4. drift: not an artifact
inner = round(462*fs):round(478*fs);
assert(mean(m(inner)) < 0.01, 'window: %.0f%% of the drift flagged', 100*mean(m(inner)))
m4 = evalc_detect(ecg, fs, 'Method','window', 'HighPass', 0);
assert(mean(m4(inner)) > 0.2, 'window without high-pass: drift not flagged (check is void)')
fprintf('4 OK: drift flagged %.0f%% (%.0f%% without the high-pass)\n', 100*mean(m(inner)), 100*mean(m4(inner)));

% 5. a stretch of taller heartbeats (e.g. after a change of posture): not an artifact
tall = round(150*fs):round(270*fs);
ecg5 = ecg; ecg5(tall) = ecg5(tall) + 2 * beats(tall);   % heartbeats 3x, same noise
mid  = round(160*fs):round(260*fs);
m5 = evalc_detect(ecg5, fs, 'Method','window', 'HFthresh', Inf);   % the range rule alone
assert(mean(m5(mid)) < 0.05, 'window: %.0f%% of the taller heartbeats flagged', 100*mean(m5(mid)))
m5g = evalc_detect(ecg5, fs, 'Method','window', 'HFthresh', Inf, 'AmpLocalSec', 1e6);   % one reference for all
assert(mean(m5g(mid)) > 0.5, 'window with the recording median: taller heartbeats not flagged (check is void)')
fprintf('5 OK: taller heartbeats flagged %.0f%% (%.0f%% against the recording median)\n', 100*mean(m5(mid)), 100*mean(m5g(mid)));
end

function [m, seg] = evalc_detect(x, fs, varargin)
[~, m, seg] = evalc('detect_ecg_artifacts(x, ''SampleRate'', fs, varargin{:})');
end
