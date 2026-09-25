% COMPUTE_HEP_TF - Heartbeat-locked time-frequency measures (HRSP, HEPC) and
% the surrogate heartbeat control, from continuous (cleaned) EEG.
%
% HRSP (heartbeat-related spectral perturbation) and HEPC (heartbeat-evoked
% phase coupling) are computed from a zero-phase complex Morlet wavelet
% transform of the continuous signals, indexed at every heartbeat (Lee et
% al., 2024):
%   power - mean power across heartbeats (uV^2; wavelets normalized so that a
%           sinusoid of amplitude A gives A^2)
%   hrsp  - power in dB relative to its mean over the epoch window (the whole
%           cardiac cycle), for each channel and frequency
%   hepc  - pairwise phase consistency across heartbeats (Vinck et al., 2010):
%           unbiased, so its expected value is 0 for random phases whatever
%           the number of heartbeats (unlike ITC)
%
% Surrogate control (opts.nSurr > 0): the whole heartbeat train is shifted
% rigidly by a random delay (default -500 to 500 ms, excluding shifts shorter
% than a quarter of the median inter-beat interval), which keeps the number of
% beats and the inter-beat intervals but moves them to other cardiac phases
% (Park et al., 2018; Lee et al., 2024). The HEP, HRSP and HEPC are recomputed
% for each surrogate, and each point (channel x [frequency x] latency) is
% compared with its surrogate distribution at the same latency: z-score,
% p-value from the z-score (two-sided for HEP and HRSP, one-sided for HEPC),
% FDR-corrected across all points (Benjamini & Hochberg), and the empirical
% p-value (smallest possible: 1/(nSurr+1)). The question answered is whether
% the activity at a given latency after the heartbeat differs from the
% activity at other cardiac phases. A shifted train is still heartbeat-locked
% at another phase, so the maximum over latencies is not a valid null here.
% Note: the cardiac field artifact is heartbeat-locked too, so the control
% does not test a neural origin.
%
% Usage:
%   [tf, surr] = compute_hep_tf(X, fs, beats, win, opts)
%
% Inputs:
%   X     - continuous data (channels x samples)
%   fs    - sampling rate (Hz)
%   beats - heartbeat sample indices (the epochs' time-locking events)
%   win   - epoch window [start end] in ms, e.g. [-300 600]
%   opts  - optional struct:
%           .freqs      frequencies in Hz (default 4:30)
%           .cycles     wavelet cycles (default 5)
%           .tstep      output time step in ms (default ~10)
%           .tf         compute HRSP and HEPC (default true)
%           .nSurr      number of surrogates (default 0 = none)
%           .shift      surrogate shift range in ms (default [-500 500])
%           .boundaries sample positions of data discontinuities (EEGLAB
%                       'boundary' events): windows crossing them are excluded
%           .seed       random seed for the surrogate shifts (default 1)
%           .hep_times  time points of the HEP epochs in ms (default: every
%                       sample of win), e.g. HEP.times
%
% Outputs:
%   tf   - struct: times (ms), freqs (Hz), power, hrsp, hepc (channels x freqs
%          x times), hep (channels x times, the average of the same beats),
%          nBeats, cycles
%   surr - struct (empty without surrogates): nSurr, shifts (ms), and for hep,
%          hrsp, hepc: z, p, p_fdr, p_emp (see above), null_mean, null_sd,
%          null_lo/null_hi (2.5th and 97.5th percentiles of the surrogates)
%
% References:
%   Lee et al. (2024). Heartbeat-related spectral perturbation of
%   electroencephalogram reflects dynamic interoceptive attention states in
%   the trial-by-trial classification analysis. NeuroImage.
%   Park et al. (2018). Heartbeat-evoked cortical responses... (surrogate R-peaks)
%   Vinck et al. (2010). The pairwise phase consistency: a bias-free measure of
%   rhythmic neuronal synchronization. NeuroImage.
%
% Copyright (C) - Cedric Cannard, 2026

function [tf, surr] = compute_hep_tf(X, fs, beats, win, opts)

if nargin < 5, opts = struct; end
freqs  = getdef(opts, 'freqs', 4:30);
cycles = getdef(opts, 'cycles', 5);
tstep  = getdef(opts, 'tstep', 10);
doTF   = getdef(opts, 'tf', true);
nSurr  = getdef(opts, 'nSurr', 0);
shiftR = getdef(opts, 'shift', [-500 500]);
bnd    = getdef(opts, 'boundaries', []);
seed   = getdef(opts, 'seed', 1);

X = double(X);
[nChan, nPts] = size(X);
beats = round(beats(:));
nF = numel(freqs);

% Output time points (every k samples, ~tstep ms)
k = max(1, round(tstep/1000*fs));
tIdx = round(win(1)/1000*fs):k:round(win(2)/1000*fs)-1;
times = tIdx / fs * 1000;
hepIdx = round(win(1)/1000*fs):round(win(2)/1000*fs)-1;   % every sample for the HEP
if isfield(opts,'hep_times') && ~isempty(opts.hep_times)
    hepIdx = round(opts.hep_times(:)'/1000*fs);            % the epochs' own time points
end

% Wavelet support (samples): windows must stay in the data and not cross a
% discontinuity
sig = cycles ./ (2*pi*freqs);                  % temporal SD of each wavelet (s)
half = ceil(3*max(sig)*fs);
if ~doTF, half = 0; end
okBeat = @(b) valid_beats(b, [hepIdx(1) hepIdx(end)], half, nPts, bnd);

v = okBeat(beats);
if any(~v)
    fprintf('compute_hep_tf: %g/%g heartbeats too close to the edges or to a discontinuity were left out. \n', sum(~v), numel(v));
end
beats = beats(v);
nB = numel(beats);
if nB < 10
    error('compute_hep_tf: fewer than 10 heartbeats.')
end

% Surrogate shifts (same rigid shift of the whole train for all channels)
shifts = [];
if nSurr > 0
    rs = RandStream('twister','Seed',seed);
    lo = round(shiftR(1)/1000*fs); hi = round(shiftR(2)/1000*fs);
    minAbs = round(0.25*median(diff(beats)));
    minAbs = min(minAbs, floor(0.5*max(abs([lo hi]))));
    shifts = zeros(nSurr,1);
    for s = 1:nSurr
        d = 0;
        while abs(d) < minAbs
            d = lo + floor(rand(rs)*(hi-lo+1));
        end
        shifts(s) = d;
    end
    surrBeats = cell(nSurr,1);
    for s = 1:nSurr
        b = beats + shifts(s);
        surrBeats{s} = b(okBeat(b));
    end
end

% HEP: average of the same beats (and surrogates)
tf.hep = mean_epochs(X, beats, hepIdx);
tf.hep_times = hepIdx / fs * 1000;
surr = struct([]);
if nSurr > 0
    surr = struct('nSurr', nSurr, 'shifts', shifts/fs*1000);
    nullHep = zeros([size(tf.hep) nSurr]);
    for s = 1:nSurr
        nullHep(:,:,s) = mean_epochs(X, surrBeats{s}, hepIdx);
    end
    surr.hep = surrogate_stats(tf.hep, nullHep, 'both');
    clear nullHep
end

tf.times = times; tf.freqs = freqs; tf.cycles = cycles; tf.nBeats = nB;
if ~doTF
    return
end

% Zero-phase Morlet kernels (FFT), amplitude-normalized
nfft = 2^nextpow2(nPts + 2*half + 1);
W = complex(zeros(nF, nfft));
for iF = 1:nF
    tk = -ceil(3*sig(iF)*fs):ceil(3*sig(iF)*fs);            % samples, centered on 0
    w = exp(-(tk/fs).^2/(2*sig(iF)^2)) .* exp(1i*2*pi*freqs(iF)*tk/fs);
    w = 2 * w / sum(abs(w));                                 % sinusoid of amplitude A -> |coef| = A
    wp = zeros(1,nfft);
    wp(1:numel(w)) = w;
    wp = circshift(wp, -(numel(tk)-1)/2);                    % kernel center at index 1: no group delay
    W(iF,:) = fft(wp);
end

tf.power = nan(nChan, nF, numel(tIdx));
tf.hepc  = nan(nChan, nF, numel(tIdx));
if nSurr > 0
    nullP = zeros(nChan,nF,numel(tIdx),nSurr,'single'); nullC = nullP;
end
for iChan = 1:nChan
    FX = fft(X(iChan,:), nfft);
    for iF = 1:nF
        co = ifft(FX .* W(iF,:));
        co = co(1:nPts);
        [pw, pc] = tf_stats(co, beats, tIdx);
        tf.power(iChan,iF,:) = pw;
        tf.hepc(iChan,iF,:) = pc;
        if nSurr > 0
            for s = 1:nSurr
                [pwS, pcS] = tf_stats(co, surrBeats{s}, tIdx);
                nullP(iChan,iF,:,s) = pwS;
                nullC(iChan,iF,:,s) = pcS;
            end
        end
    end
end
tf.hrsp = 10*log10(tf.power ./ mean(tf.power,3));

if nSurr > 0
    % HRSP of each surrogate, relative to its own window mean
    nullH = 10*log10(double(nullP) ./ mean(double(nullP),3));
    clear nullP
    surr.hrsp = surrogate_stats(tf.hrsp, nullH, 'both');
    clear nullH
    surr.hepc = surrogate_stats(tf.hepc, double(nullC), 'right');
end
end

%% Subfunctions
function v = valid_beats(b, lims, half, nPts, bnd)
% Beats whose window (+ wavelet support) stays in the data without crossing
% a discontinuity
a = b + lims(1) - half;  z = b + lims(2) + half;
v = a >= 1 & z <= nPts;
for i = 1:numel(bnd)
    v = v & ~(bnd(i) > a & bnd(i) < z);
end
end

function m = mean_epochs(X, beats, idx)
% Average of the epochs X(:, beats + idx) (channels x times)
m = zeros(size(X,1), numel(idx));
for i = 1:numel(beats)
    m = m + X(:, beats(i) + idx);
end
m = m / numel(beats);
end

function [pw, pc] = tf_stats(co, beats, tIdx)
% Mean power and pairwise phase consistency across beats (1 x times)
C = co(beats + tIdx);                   % beats x times
if isvector(C) && numel(beats) > 1, C = C(:)'; end
N = numel(beats);
pw = mean(abs(C).^2, 1);
z = C ./ abs(C);
pc = (abs(sum(z,1)).^2 - N) / (N*(N-1));
end

function S = surrogate_stats(real, null, tail)
% Compare each point with its surrogate distribution (last dimension)
nd = ndims(null); nSurr = size(null, nd);
S.null_mean = mean(null, nd);
S.null_sd = std(null, 0, nd);
S.null_lo = prctile(null, 2.5, nd);
S.null_hi = prctile(null, 97.5, nd);
S.z = (real - S.null_mean) ./ S.null_sd;
dev = abs(null - S.null_mean);
if strcmp(tail, 'both')
    S.p = erfc(abs(S.z)/sqrt(2));
    S.p_emp = (sum(dev >= abs(real - S.null_mean), nd) + 1) / (nSurr + 1);
else
    S.p = 0.5*erfc(S.z/sqrt(2));
    S.p_emp = (sum(null >= real, nd) + 1) / (nSurr + 1);
end
S.p_fdr = fdr_bh(S.p);
end

function q = fdr_bh(p)
% Benjamini-Hochberg adjusted p-values (same size as p; NaN ignored)
q = nan(size(p));
v = find(~isnan(p));
[ps, ord] = sort(p(v));
m = numel(ps);
adj = min(1, flipud(cummin(flipud(ps(:) .* m ./ (1:m)'))));
q(v(ord)) = adj;
end

function v = getdef(s, f, d)
if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
