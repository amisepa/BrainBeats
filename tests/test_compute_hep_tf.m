function test_compute_hep_tf
% TEST_COMPUTE_HEP_TF - Checks of compute_hep_tf on synthetic data (no EEGLAB
% needed). Errors on the first failed check.
%   1. No group delay: a 10-Hz burst centered 300 ms after each heartbeat
%      peaks at 300 ms in the HRSP at every frequency of the burst.
%   2. HRPC is high for a phase-locked burst and ~0 for random phases (same
%      power), and the HRSP sees both.
%   3. Surrogate control: an evoked deflection at 400 ms is significant
%      (FDR-corrected), and pure noise gives no significant point. The null
%      band of a noise channel is at the noise level of an average of nB
%      windows.
%   4. Heartbeats whose window crosses a discontinuity are left out.
%   5. HRPC p-values: pure noise gives (almost) no FDR-significant point, and
%      a phase-locked burst is detected.
%   6. HEP family-wise statistics (TFCE and max statistic, 1000 surrogates):
%      the evoked deflection of check 3 is significant, the noise channel is
%      not. (False-positive rate and power on simulated noise: README.)
%   7. Induced power (HEP subtracted): the phase-locked burst of check 1
%      leaves no induced power, the random-phase burst keeps it.
%   8. A QRS-like transient: the total HRSP shows a significant "decrease"
%      away from it (its power in the window mean), the induced HRSP does not.
%   6b. TFCE and max statistic equal a direct computation; a weak negative
%      deflection is detected with the right sign.
%   9. Pre-R-peak band power per heartbeat recovers the amplitude of an
%      oscillation that changes from beat to beat; a phase-locked burst adds
%      no induced power.
%   10. Family-wise HRSP/HRPC statistics on white noise (and the FDR rate).
%   11. One output latency; fewer HEP surrogates than HRSP ones; the default
%      statistic is the max statistic (TFCE only on request).
%   12. Excluded spans: surrogates avoid them and keep as many beats as the HEP.

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
c1 = mean(tf.hrpc(1,f10,i300)); c2 = mean(tf.hrpc(2,f10,i300));
assert(c1 > 0.5 && abs(c2) < 0.05, 'HRPC phase-locked %.2f (expected > .5), random %.2f (expected ~0)', c1, c2)
for f = [8 12]
    [~,iM] = max(squeeze(tf.hrsp(1,tf.freqs==f,:)));
    assert(abs(tf.times(iM) - 300) <= 20, 'HRSP peak at %g ms at %g Hz (group delay?)', tf.times(iM), f)
end
fprintf('1-2 OK: HRSP peak at %g ms; HRPC %.2f (locked) vs %.3f (random)\n', tf.times(iMax), c1, c2);

% 7. induced power (HEP subtracted from each heartbeat): the phase-locked
% burst leaves nothing, the random-phase burst is unchanged (from -200 ms:
% before, the wavelets reach the previous beat's burst, at -500 ms +/- the
% inter-beat jitter, which is not locked to this beat)
ind1 = squeeze(tf.hrsp_induced(1,f10,tf.times >= -200)); ind2 = squeeze(tf.hrsp_induced(2,f10,:));
tot2 = squeeze(tf.hrsp(2,f10,:));
assert(max(abs(ind1)) < 1, 'Induced HRSP of the phase-locked burst: %.1f dB (expected ~0)', max(abs(ind1)))
[~, iMax3] = max(ind2);
assert(abs(tf.times(iMax3) - 300) <= 20 && abs(max(ind2) - max(tot2)) < 0.5, ...
    'Induced HRSP of the random-phase burst: peak %.1f dB at %g ms (total %.1f dB)', max(ind2), tf.times(iMax3), max(tot2))
fprintf('7 OK: induced HRSP of the phase-locked burst max |%.2f| dB; random-phase burst %.1f dB (total %.1f)\n', ...
    max(abs(ind1)), max(ind2), max(tot2));

% 3. surrogate control: evoked deflection at 400 ms on channel 1, noise on 2
X = randn(2,nPts);
for b = beats'
    X(1,b + round(0.4*fs) + (-5:5)) = X(1,b + round(0.4*fs) + (-5:5)) + 1.5*hann(11)';
end
[tf, surr] = compute_hep_tf(X, fs, beats, [-300 600], struct('freqs',6:10,'nSurr',100,'tf',false,'stat','tfce'));
i400 = abs(tf.hep_times - 400) < 5;
assert(all(surr.hep.p_fdr(1,i400) < .05), 'Evoked deflection not significant after correction')
assert(~any(surr.hep.p_fdr(2,:) < .05), 'Pure noise gave an FDR-significant point')
sdNoise = median(surr.hep.null_sd(2,:)); expected = 1/sqrt(tf.nBeats);
assert(abs(sdNoise/expected - 1) < 0.25, 'Null SD %.3f, expected ~%.3f', sdNoise, expected)
fprintf('3 OK: evoked peak z = %.1f (FDR p = %.1e, empirical p = %.3f); noise: min FDR p = %.2f\n', ...
    min(surr.hep.z(1,i400)), max(surr.hep.p_fdr(1,i400)), max(surr.hep.p_emp(1,i400)), min(surr.hep.p_fdr(2,:)));
assert(all(surr.hep.p_tfce(1,i400) < .05) && all(surr.hep.p_tmax(1,i400) < .05), ...
    'Evoked deflection not significant with TFCE or the max statistic')
assert(~any(surr.hep.p_tfce(2,:) < .05 | surr.hep.p_tmax(2,:) < .05), ...
    'Pure noise significant with TFCE or the max statistic')
fprintf('6 OK: TFCE and max statistic (%d surrogate HEPs): evoked deflection p = %.3f and %.3f; noise channel: none\n', ...
    surr.nSurrHep, max(surr.hep.p_tfce(1,i400)), max(surr.hep.p_tmax(1,i400)));

% 6b. family-wise statistics against a direct computation (TFCE: E = 0.5,
% H = 2, dh = 0.1 along time, both signs; max statistic), with a weak
% negative deflection at 200 ms that TFCE must detect
X = randn(3,nPts);
for b = beats'
    X(2,b + round(0.2*fs) + (-12:12)) = X(2,b + round(0.2*fs) + (-12:12)) - 0.25*hann(25)';
end
[tf, surr] = compute_hep_tf(X, fs, beats, [-300 600], struct('tf',false,'nSurr',200,'keepNull',true,'stat','tfce'));
H = surr.hep; z = H.z;
T = zeros(size(z));
for sgn = [1 -1]
    for h = 0.1:0.1:max(sgn*z(:))
        for c = 1:size(z,1)
            m = sgn*z(c,:) >= h;
            d = diff([0 m 0]); on = find(d == 1); off = find(d == -1) - 1;
            for k = 1:numel(on)
                T(c,on(k):off(k)) = T(c,on(k):off(k)) + sgn * h^2 * 0.1 * (off(k)-on(k)+1)^0.5;
            end
        end
    end
end
assert(max(abs(T(:) - H.tfce(:))) < 1e-9 * max(abs(T(:))), 'TFCE differs from the direct computation')
nS = size(H.null,3); zn = zeros(size(H.null));
for j = 1:nS                                  % leave-one-out standardization
    o = H.null(:,:,[1:j-1 j+1:nS]);
    zn(:,:,j) = (H.null(:,:,j) - mean(o,3)) ./ std(o,0,3);
end
mx = squeeze(max(abs(zn), [], [1 2]));
pt = (sum(mx' >= abs(z(:)), 2) + 1) / (numel(mx) + 1);
assert(max(abs(pt - H.p_tmax(:))) < 1e-12, 'Max-statistic p-values differ from the direct computation')
i200 = abs(tf.hep_times - 200) <= 20;
assert(any(H.p_tfce(2,i200) < .05) && all(H.tfce(2,i200 & H.p_tfce(2,:) < .05) < 0), ...
    'Weak negative deflection not detected (or with the wrong sign) by TFCE')
fprintf('6b OK: TFCE and max statistic match the direct computation; negative deflection (min z %.1f) detected, TFCE p = %.3f\n', ...
    min(z(2,i200)), min(H.p_tfce(2,i200)));

% 4. discontinuity: beats around it are left out
tf2 = compute_hep_tf(X, fs, beats, [-300 600], struct('tf',false,'boundaries',beats(10)+10));
assert(tf2.nBeats < numel(beats), 'Heartbeat crossing a discontinuity was kept')
fprintf('4 OK: %d/%d heartbeats kept with one discontinuity\n', tf2.nBeats, numel(beats));

% 5. HRPC calibration (a Gaussian tail gave 2.5-4% false points here)
o = struct('freqs',6:2:14,'nSurr',100);
hits = zeros(1,5);
for r = 1:5
    rng(10+r); o.seed = r;
    [~, surr] = compute_hep_tf(randn(1,nPts), fs, beats, [-300 600], o);
    hits(r) = mean(surr.hrpc.p_fdr(:) < .05);
end
assert(mean(hits) < .005, 'HRPC: %.2f%% FDR-significant points on pure noise', 100*mean(hits))
X = randn(1,nPts);
for b = beats'
    idx = b + (-0.3*fs:0.6*fs);
    X(idx) = X(idx) + 0.5*env .* cos(2*pi*10*t);
end
[tf, surr] = compute_hep_tf(X, fs, beats, [-300 600], o);
i10 = tf.freqs == 10; i300 = abs(tf.times - 300) <= 10;
assert(all(surr.hrpc.p_fdr(1,i10,i300) < .05), 'HRPC: phase-locked burst not detected')
fprintf('5 OK: HRPC FDR-significant points on noise %.2f%%; phase-locked burst detected\n', 100*mean(hits));

% 8. a QRS-like transient (20 ms, 5 x the noise SD, amplitude varying 15%
% across beats as with breathing) on noise: the total HRSP shows a "decrease"
% away from it that the surrogates flag; the induced HRSP does not
X = randn(1,nPts);
q = exp(-((-12:12)/fs).^2/(2*0.004^2));            % ~20 ms wide
for b = beats'
    X(b + (-12:12)) = X(b + (-12:12)) + 5*(1 + 0.15*randn)*q;
end
[tf, surr] = compute_hep_tf(X, fs, beats, [-300 600], struct('freqs',4:2:30,'nSurr',100));
% points farther than 2 wavelet SDs from this QRS and from the next one
% (whose latency varies with the inter-beat interval)
sd2 = 2 * tf.cycles ./ (2*pi*tf.freqs(:)) * 1000;   % ms, freqs x 1
away = tf.times > sd2 & tf.times < min(diff(beats))/fs*1000 - sd2;
pT = squeeze(surr.hrsp.p_fdr(1,:,:)); pI = squeeze(surr.hrsp_induced.p_fdr(1,:,:));
hT = squeeze(tf.hrsp(1,:,:)); hI = squeeze(tf.hrsp_induced(1,:,:));
fT = mean(pT(away) < .05); fI = mean(pI(away) < .05);
fprintf('8: QRS transient, HRSP away from the QRSs (%d points): total %.2f dB (%.0f%% FDR points), induced %.2f dB (%.1f%%)\n', ...
    sum(away(:)), median(hT(away)), 100*fT, median(hI(away)), 100*fI);
assert(fT > .1, 'Total HRSP: the QRS "decrease" was expected (%.0f%% FDR points)', 100*fT)
assert(fI < .005, 'Induced HRSP: %.1f%% FDR points away from the QRS', 100*fI)
fprintf('8 OK\n');

% 9. pre-R-peak power per heartbeat: a 10-Hz oscillation whose amplitude
% changes every 10 heartbeats (1 or 3) is recovered in the alpha band (9 x
% more power at 3; the band mean of the wavelet responses to a 10-Hz
% sinusoid of amplitude A is ~0.65 A^2 for 8-12 Hz), and not in beta
amp = 1 + 2*mod(floor((0:numel(beats))/10), 2);   % per inter-beat interval
X = 0.2*randn(1,nPts); tt = (0:nPts-1)/fs;
edgesB = [1; beats; nPts];
for i = 1:numel(edgesB)-1
    seg = edgesB(i):edgesB(i+1)-1;
    X(seg) = X(seg) + amp(i) * sin(2*pi*10*tt(seg) + 2*pi*rand);
end
pre = struct('win',[-300 -50], 'bands',[8 13; 13 30], 'names',{{'alpha' 'beta'}});
tf = compute_hep_tf(X, fs, beats, [-300 600], struct('freqs',4:30,'prestim',pre));
[~, ib] = ismember(tf.beats(:)', beats);
a2 = amp(ib).^2;                                   % amplitude^2 before each beat
inner = mod(ib-1, 10) >= 2 & mod(ib-1, 10) <= 8;  % not next to an amplitude change
pA = squeeze(tf.prestim.power(1,1,:))'; pB = squeeze(tf.prestim.power(1,2,:))';
k9 = median(pA(inner & a2 == 9)) / median(pA(inner & a2 == 1));
rB = median(pB ./ pA);
assert(abs(k9/9 - 1) < .25 && rB < .1 && corr(pA(:), a2(:)) > .9, ...
    'Pre-R alpha power ratio %.1f between amplitudes 3 and 1 (expected 9), beta/alpha %.2f (expected ~0)', k9, rB)
fprintf('9 OK: pre-R alpha power ratio %.1f between amplitudes 3 and 1 (expected 9), %.2f x A^2; beta/alpha %.3f\n', ...
    k9, median(pA(inner) ./ a2(inner)), rB);
% a 10-Hz burst phase-locked to the heartbeats before the R-peak adds no
% induced power (it is part of the HEP)
Xl = X;
for b = beats'
    idx = b + (-0.3*fs:0);
    Xl(idx) = Xl(idx) + 3 * hann(numel(idx))' .* cos(2*pi*10*(-0.3*fs:0)/fs);
end
tfl = compute_hep_tf(Xl, fs, beats, [-300 600], struct('freqs',4:30,'prestim',pre));
pL = squeeze(tfl.prestim.power(1,1,:))';
assert(abs(median(pL ./ pA) - 1) < .1, 'Phase-locked burst changed the induced pre-R power (x %.2f)', median(pL ./ pA))
fprintf('9b OK: phase-locked pre-R burst: induced alpha power x %.3f\n', median(pL ./ pA));

% 10. family-wise statistics of HRSP and HRPC on white noise (6 channels,
% 27 frequencies, 100 surrogates, 3 runs): FDR is liberal for HRSP, the
% family-wise statistics are not
anyF = zeros(3,3); anyW = zeros(3,3);
for r = 1:3
    rng(100+r);
    [~, s] = compute_hep_tf(randn(6,nPts), fs, beats, [-300 600], struct('nSurr',100,'seed',r,'stat','tfce'));
    M = {s.hrsp, s.hrsp_induced, s.hrpc};
    for k = 1:3
        anyF(r,k) = any(M{k}.p_fdr(:) < .05);
        anyW(r,k) = any(M{k}.p_tfce(:) < .05 | M{k}.p_tmax(:) < .05);
    end
end
assert(sum(anyW(:)) <= 1, 'Family-wise HRSP/HRPC statistics: significant points on noise in %d of 9 cases', sum(anyW(:)))
fprintf('10 OK: noise, runs with any significant point (HRSP / induced HRSP / HRPC): FDR %s, family-wise %s\n', ...
    mat2str(sum(anyF)), mat2str(sum(anyW)));

% 11. edge cases: one output latency; fewer HEP surrogates than HRSP ones
tf = compute_hep_tf(randn(1,nPts), fs, beats, [0 10], struct('freqs',8:10,'tstep',10));
assert(numel(tf.times) == 1 && all(isfinite(tf.hrsp(:))), 'One output latency failed')
[~, s] = compute_hep_tf(randn(1,nPts), fs, beats, [-300 600], struct('freqs',8:10,'nSurr',30,'nSurrHep',10));
assert(s.nSurrHep >= 30, 'nSurrHep below nSurr')
% default statistic: max statistic, TFCE not computed
assert(strcmp(s.stat,'tmax') && isfield(s.hrpc,'p_tmax') && ~isfield(s.hep,'p_tfce') && ~isfield(s.hrsp,'p_tfce'), ...
    'Default statistic: TFCE computed or max statistic missing')
fprintf('11 OK: one output latency; nSurrHep raised to %d; default max statistic, no TFCE\n', s.nSurrHep);

% 12. excluded spans (epochs rejected by the cleaning): large artifacts in 13
% epoch windows; excluding them, the null is at the noise level of an
% average of as many windows as the HEP (surrogates refilled to nB beats)
x = randn(1,nPts); w = round([-0.3 0.6]*fs);
bad = beats(round(linspace(20, numel(beats)-20, 13)));
for b = bad', x(b + (w(1):w(2))) = x(b + (w(1):w(2))) + 30*randn(1, diff(w)+1); end
kb = setdiff(beats, bad);
o = struct('tf',false,'nSurr',100,'allBeats',beats);
[tfA, sA] = compute_hep_tf(x, fs, kb, [-300 600], o);
o.exclude = bad + w;
[~, sB] = compute_hep_tf(x, fs, kb, [-300 600], o);
expected = 1/sqrt(tfA.nBeats);
assert(abs(median(sB.hep.null_sd)/expected - 1) < 0.1 && median(sA.hep.null_sd) > 3*expected, ...
    'Exclusion: null SD %.3f (expected %.3f; without exclusion %.3f)', median(sB.hep.null_sd), expected, median(sA.hep.null_sd))
fprintf('12 OK: null SD %.3f with the rejected windows excluded (noise level %.3f), %.3f without\n', ...
    median(sB.hep.null_sd), expected, median(sA.hep.null_sd));
fprintf('All compute_hep_tf checks passed.\n')
