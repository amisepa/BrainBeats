function [fp, hit] = sim_hep_fwer(nRep, nChan, noiseType, amp)
% SIM_HEP_FWER - False positives and power of the HEP surrogate statistics
% of COMPUTE_HEP_TF (per-point FDR, max statistic, TFCE) on simulated data.
%
% 300 s at 250 Hz, heartbeats every 0.85 +/- 0.05 s, white noise or pink
% noise (1/f, 0.5-30 Hz), 100 surrogates (1000 for the HEP). amp > 0 adds a
% response (Hann window of 200 ms centered at 400 ms, amp x the noise SD) to
% channel 1 at every heartbeat; points outside 300-500 ms of channel 1 and on
% the other channels are then false positives. No EEGLAB needed.
%
% Usage:
%   sim_hep_fwer(200, 1, 'white', 0)     % false-positive rate, 1 channel
%   sim_hep_fwer(100, 16, 'pink', 0.15)  % power and false positives, 16 channels
%
% Outputs:
%   fp  - share of runs with any false positive (FDR, tmax, TFCE)
%   hit - share of runs where the response is detected (amp > 0)
%
% Results (2026-11-03): no response, any false positive (FDR / tmax / TFCE):
% white 1 ch .050/.055/.050 (200 runs), pink 1 ch .030/.040/.035 (200), pink
% 16 ch .010/.010/.030 (100). Response 0.15 SD, pink: 1 ch detected
% .43/.46/.51, false positives elsewhere .14/.06/.07; 16 ch detected
% .13/.13/.20, false positives .03/.01/.03 (100 runs each).
% Rerun 2026-11-04 with each surrogate standardized by the others
% (leave-one-out): white 1 ch .050/.055/.050; response 0.15 SD, pink 1 ch:
% detected .43/.45/.51, false positives elsewhere .14/.04/.07.
% Stronger responses (2026-11-04, 1 channel, 20 runs; detected in all):
% false positives elsewhere at 0.5 SD white 1.00/.40/.60, pink .85/.30/.35;
% at 1 SD 1.00/1.00/1.00 (white and pink). The surrogates sample all
% cardiac phases, so their mean at every latency is the average over the
% cardiac cycle; a large heartbeat-locked response moves the rest of the
% cycle the other way, and the latencies without a response then differ
% from that average too. Family-wise control holds when nothing is
% heartbeat-locked; with a response (or the QRS), significance away from
% it means "differs from the cycle average", not "locked activity there".
%
% Copyright (C) - Cedric Cannard, 2026

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'functions'));
fs = 250; n = 300*fs;
hit = zeros(nRep,3); fp = zeros(nRep,3); t0 = tic;
for r = 1:nRep
    rng(r)
    ibi = round((0.85 + 0.05*randn(400,1))*fs);
    beats = cumsum([fs; ibi]); beats = beats(beats < n - fs);
    X = randn(nChan, n);
    if strcmp(noiseType,'pink')
        f = (0:n-1)*fs/n; f = min(f, fs - f);
        g = zeros(1,n); keep = f >= 0.5 & f <= 30; g(keep) = 1 ./ sqrt(f(keep));
        X = real(ifft(fft(X, [], 2) .* g, [], 2));
        X = X ./ std(X, 0, 2);
    end
    if amp > 0
        w = hann(round(0.2*fs))'; c = round(0.4*fs); idx = c - floor(numel(w)/2) + (0:numel(w)-1);
        for b = beats(:)'
            X(1, b + idx) = X(1, b + idx) + amp*w;
        end
    end
    [tf, S] = compute_hep_tf(X, fs, beats, [-300 600], struct('tf',false,'nSurr',100,'allBeats',beats,'seed',r,'stat','tfce'));
    resp = tf.hep_times >= 300 & tf.hep_times <= 500;
    P = {S.hep.p_fdr, S.hep.p_tmax, S.hep.p_tfce};
    for k = 1:3
        sig = P{k} < .05;
        if amp > 0
            hit(r,k) = any(sig(1,resp));
            sig(1,resp) = false;
        end
        fp(r,k) = any(sig(:));
    end
end
fp = mean(fp); hit = mean(hit);
fprintf('%s noise, %d channel(s), response %.2f SD, %d runs (%.0f s): any false positive FDR %.3f, tmax %.3f, TFCE %.3f\n', ...
    noiseType, nChan, amp, nRep, toc(t0), fp);
if amp > 0
    fprintf('  response detected: FDR %.3f, tmax %.3f, TFCE %.3f\n', hit);
end
end
