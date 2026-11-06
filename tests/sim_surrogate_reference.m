function R = sim_surrogate_reference(nRep, cfaAmp, hepAmp)
% SIM_SURROGATE_REFERENCE - Which reference should the HEP surrogate test use?
%
% The 'random' surrogates of COMPUTE_HEP_TF compare each latency with the
% cardiac-cycle average. After a high-pass filter the data have no DC, so a
% heartbeat-locked waveform with a nonzero area (QRS, T wave) shifts every
% latency without locked activity by minus its mean: those latencies then
% differ from the surrogates. This script compares four references on
% simulated data where the latencies without locked activity are known:
%   cycle      - none (current: the surrogates' mean, the cycle average)
%   pre150     - each HEP (real and surrogate) minus its mean at -150 to -50 ms
%                (the default 'hep_baseline_win')
%   pre300     - the same at -300 to -200 ms
%   noQRS      - minus its mean over the epoch outside the QRS (-50 to 100 ms)
%
% 300 s at 250 Hz, heartbeats every 0.85 +/- 0.05 s, pink noise (1/f, 0.5-30
% Hz, SD 1), a cardiac field artifact shaped like an ECG (P -170 ms, Q, R
% (amplitude cfaAmp x the noise SD), S, T at 250 ms; Gaussian waves) and an
% optional response (Hann window of 200 ms at 400 ms, hepAmp x the noise SD),
% high-pass 0.5 Hz and low-pass 30 Hz (zero-phase), 1000 surrogate HEPs.
% False positive: a point with p < .05 where neither the cardiac field
% (below 1% of the R wave) nor the response is. No EEGLAB needed.
%
% Usage:
%   R = sim_surrogate_reference(100, 0.5, 0)      % 100 runs, R = 0.5 SD, no HEP
%   R = sim_surrogate_reference(100, 0.5, 0.15)
%
% Output: R - table, one row per reference: any false positive (share of
% runs) and response detected (share of runs) for TFCE and the max statistic.
%
% Results: see the end of this file.
%
% Copyright (C) - Cedric Cannard, 2026

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'functions'));
fs = 250; n = 300*fs; win = [-300 600];
t = (round(win(1)/1000*fs):round(win(2)/1000*fs)-1) / fs * 1000;     % HEP latencies (ms)
g = @(mu, sd) exp(-(t - mu).^2 / (2*sd^2));
cfa = cfaAmp * (0.08*g(-170,25) - 0.1*g(-25,6) + g(0,8) - 0.3*g(25,8) + 0.25*g(250,50));
resp = zeros(size(t)); inR = t >= 300 & t < 500;
resp(inR) = hepAmp * hann(sum(inR))';
quiet = abs(cfa) < .01*max(cfaAmp, eps) & ~inR;
refs = {'cycle', []; 'pre150', t >= -150 & t <= -50; 'pre300', t >= -300 & t <= -200; ...
    'noQRS', t < -50 | t > 100};
nRef = size(refs, 1);
fp = zeros(nRep, nRef, 2); hit = fp;
f = (0:n-1)*fs/n; f = min(f, fs - f);
band = min(1, max(0, (f - .25)/.25)) .* min(1, max(0, (37 - f)/7));   % 0.5-30 Hz, zero-phase
t0 = tic;
for r = 1:nRep
    rng(r)
    ibi = round((0.85 + 0.05*randn(400,1))*fs);
    beats = cumsum([fs; ibi]); beats = beats(beats < n - fs);
    x = randn(1, n);
    gp = zeros(1,n); k = f >= 0.5 & f <= 30; gp(k) = 1 ./ sqrt(f(k));
    x = real(ifft(fft(x) .* gp)); x = x / std(x);
    w = cfa + resp; idx = round(t/1000*fs);
    for b = beats(:)'
        x(b + idx) = x(b + idx) + w;
    end
    x = real(ifft(fft(x) .* band));                    % no DC, like filtered EEG
    [tf, S] = compute_hep_tf(x, fs, beats, win, struct('tf',false,'nSurr',100,'nSurrHep',1000, ...
        'allBeats',beats,'seed',r,'keepNull',true));
    for k = 1:nRef
        h = tf.hep; N = S.hep.null;                    % 1 x times, 1 x times x surrogates
        if ~isempty(refs{k,2})
            h = h - mean(h(refs{k,2}));
            N = N - mean(N(:,refs{k,2},:), 2);
        end
        [pT, pM] = fwer(h, N);
        for s = 1:2
            P = {pT, pM}; sig = P{s} < .05;
            fp(r,k,s) = any(sig(quiet));
            hit(r,k,s) = any(sig(inR));
        end
    end
end
R = table(refs(:,1), mean(fp(:,:,1))', mean(fp(:,:,2))', mean(hit(:,:,1))', mean(hit(:,:,2))', ...
    'VariableNames', {'reference' 'fp_tfce' 'fp_tmax' 'hit_tfce' 'hit_tmax'});
fprintf('R wave %.2f SD, response %.2f SD, %d runs (%.0f s), %.0f%% of the latencies quiet:\n', ...
    cfaAmp, hepAmp, nRep, toc(t0), 100*mean(quiet));
disp(R)
end

function [pT, pM] = fwer(h, N)
% TFCE and max-statistic p-values of h (1 x times) against the surrogates N
% (1 x times x nS), each surrogate standardized by the others (as in
% COMPUTE_HEP_TF)
nS = size(N, 3);
mu = mean(N, 3); sd = std(N, 0, 3);
z = (h - mu) ./ sd;
d = reshape(N - mu, [], nS)';                         % surrogates x times
SS = (nS-1) * sd.^2;
zn = (d * nS/(nS-1)) ./ sqrt(max(SS - d.^2 * nS/(nS-1), 0) / (nS-2));
maxZ = max(abs(zn), [], 2);
maxT = max(abs(tfce(zn)), [], 2);
pM = (sum(maxZ' >= abs(z(:)), 2) + 1)' / (nS + 1);
pT = (sum(maxT' >= abs(tfce(z))', 2) + 1)' / (nS + 1);
end

function T = tfce(Z)
% TFCE along time of each row (E = 0.5, H = 2, dh = 0.1), signed
T = zeros(size(Z));
for sgn = [1 -1]
    Zs = max(sgn*Z, 0);
    for h = 0.1:0.1:max(Zs(:))
        M = Zs >= h;
        c = cumsum(M, 2); L = c - cummax(c .* ~M, 2);
        c = cumsum(M, 2, 'reverse'); L = (L + c - cummax(c .* ~M, 2, 'reverse') - 1) .* M;
        T = T + sgn * h^2 * 0.1 * sqrt(L) .* M;
    end
end
end

% Results (2026-11-05, 100 runs each; any false positive TFCE/tmax, then
% response detected TFCE/tmax):
%   R wave  response | cycle          | pre150         | pre300         | noQRS
%   0       0        | .02/.01        | .08/.07        | .08/.05        | .02/.02
%   0.5     0        | .08/.07        | .08/.07        | .04/.02        | .08/.05
%   1       0        | .29/.25        | .09/.06        | .04/.02        | .17/.13
%   0.5     0.15     | .19/.15 .29/.26| .09/.06 .38/.29| .04/.02 .38/.27| .16/.13 .37/.33
%   1       0.15     | .43/.34 .24/.17| .09/.06 .48/.35| .04/.02 .53/.35| .37/.31 .35/.26
%   0       0.5      | .47/.38 1/1    | .05/.07 1/1    | .05/.03 1/1    | .62/.51 1/1
% Against the cycle average, a heartbeat-locked waveform with an area (here
% the T wave and R; or a strong response) makes quiet latencies differ and
% weakens the response; a pre-R reference fixes both when the reference
% window is quiet (pre150 holds the end of the simulated P wave: 6-9%).
% On real recordings it is not reliably quiet: see CHECK_SURROGATE_REFERENCE.
