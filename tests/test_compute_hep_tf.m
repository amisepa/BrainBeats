function test_compute_hep_tf
% TEST_COMPUTE_HEP_TF - Checks of compute_hep_tf on synthetic data (no EEGLAB
% needed). Errors on the first failed check.
%   1. No group delay: a 10-Hz burst centered 300 ms after each heartbeat
%      peaks at 300 ms in the HRSP at every frequency of the burst.
%   2. HEPC is high for a phase-locked burst and ~0 for random phases (same
%      power), and the HRSP sees both.
%   3. Surrogate control: an evoked deflection at 400 ms is significant
%      (FDR-corrected), and pure noise gives no significant point.
%   4. Heartbeats whose window crosses a discontinuity are left out.

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'functions'));
fs = 250; nPts = 250*240; rng(3);
beats = (2*fs:round(0.8*fs):nPts-2*fs)' + randi([-10 10], numel(2*fs:round(0.8*fs):nPts-2*fs), 1);
t = (-0.3*fs:0.6*fs)/fs;
env = exp(-(t-0.3).^2/(2*0.05^2));                 % burst envelope, peak at 300 ms

% 1-2. phase-locked burst on channel 1, random-phase burst on channel 2
X = 0.5*randn(2,nPts);
for b = beats'
    idx = b + (-0.3*fs:0.6*fs);
    X(1,idx) = X(1,idx) + env .* cos(2*pi*10*t);
    X(2,idx) = X(2,idx) + env .* cos(2*pi*10*t + 2*pi*rand);
end
tf = compute_hep_tf(X, fs, beats, [-300 600], struct('freqs',6:14));
f10 = tf.freqs == 10;
[~,iMax] = max(squeeze(tf.hrsp(1,f10,:)));
assert(abs(tf.times(iMax) - 300) <= 10, 'HRSP peak at %g ms (expected 300 ms)', tf.times(iMax))
[~,iMax2] = max(squeeze(tf.hrsp(2,f10,:)));
assert(abs(tf.times(iMax2) - 300) <= 20, 'Random-phase HRSP peak at %g ms (expected 300 ms)', tf.times(iMax2))
i300 = abs(tf.times - 300) <= 10;
c1 = mean(tf.hepc(1,f10,i300)); c2 = mean(tf.hepc(2,f10,i300));
assert(c1 > 0.5 && abs(c2) < 0.05, 'HEPC phase-locked %.2f (expected > .5), random %.2f (expected ~0)', c1, c2)
for f = [8 12]
    [~,iM] = max(squeeze(tf.hrsp(1,tf.freqs==f,:)));
    assert(abs(tf.times(iM) - 300) <= 20, 'HRSP peak at %g ms at %g Hz (group delay?)', tf.times(iM), f)
end
fprintf('1-2 OK: HRSP peak at %g ms; HEPC %.2f (locked) vs %.3f (random)\n', tf.times(iMax), c1, c2);

% 3. surrogate control: evoked deflection at 400 ms on channel 1, noise on 2
X = randn(2,nPts);
for b = beats'
    X(1,b + round(0.4*fs) + (-5:5)) = X(1,b + round(0.4*fs) + (-5:5)) + 1.5*hann(11)';
end
[tf, surr] = compute_hep_tf(X, fs, beats, [-300 600], struct('freqs',6:10,'nSurr',100,'tf',false));
i400 = abs(tf.hep_times - 400) < 5;
assert(all(surr.hep.p_fdr(1,i400) < .05), 'Evoked deflection not significant after correction')
assert(~any(surr.hep.p_fdr(2,:) < .05), 'Pure noise gave an FDR-significant point')
assert(all(abs(surr.shifts) >= 0.25*median(diff(beats))/fs*1000 - 1e-9), 'Surrogate shift too small')
fprintf('3 OK: evoked peak z = %.1f (FDR p = %.1e, empirical p = %.3f); noise: min FDR p = %.2f\n', ...
    min(surr.hep.z(1,i400)), max(surr.hep.p_fdr(1,i400)), max(surr.hep.p_emp(1,i400)), min(surr.hep.p_fdr(2,:)));

% 4. discontinuity: beats around it are left out
tf2 = compute_hep_tf(X, fs, beats, [-300 600], struct('tf',false,'boundaries',beats(10)+10));
assert(tf2.nBeats < numel(beats), 'Heartbeat crossing a discontinuity was kept')
fprintf('4 OK: %d/%d heartbeats kept with one discontinuity\n', tf2.nBeats, numel(beats));
fprintf('All compute_hep_tf checks passed.\n')
