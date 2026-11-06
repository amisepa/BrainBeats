% COMPUTE_HEP_TF - Heartbeat-locked time-frequency measures (HRSP, HRPC) and
% the surrogate heartbeat control, from continuous (cleaned) EEG.
%
% HRSP (heartbeat-related spectral perturbation) and HRPC (heartbeat-evoked
% phase coupling) are computed from a zero-phase complex Morlet wavelet
% transform of the continuous signals, indexed at every heartbeat (Lee et
% al., 2024):
%   power - mean power across heartbeats (uV^2; wavelets normalized so that a
%           sinusoid of amplitude A gives A^2)
%   hrsp  - power in dB relative to its mean over the epoch window (the whole
%           cardiac cycle), for each channel and frequency. A heartbeat-locked
%           broadband transient (the QRS of the cardiac field artifact) raises
%           that mean, so latencies away from it show a power decrease, which
%           the 'random' surrogates (whose windows contain the QRS at random
%           latencies) find significant: normalizing by the mean of the
%           continuous data gives the same (-1.9 vs -1.7 dB for a transient
%           5 x the noise SD; none at 2 x)
%   power_induced, hrsp_induced - the same for the power that is not phase-
%           locked to the heartbeats: each heartbeat's wavelet coefficients
%           minus their mean across heartbeats, i.e. the HEP (cardiac field
%           artifact included) subtracted from every heartbeat before the
%           power is taken (Cohen & Donner, 2013). Computed as
%           N/(N-1) * (power - |mean coefficient|^2) over the whole wavelet
%           support, so the subtraction has no edge effect. A heartbeat-locked
%           transient of constant shape (the QRS) leaves no induced power, so
%           it no longer raises the window mean (see hrsp above).
%   hrpc  - pairwise phase consistency across heartbeats (Vinck et al., 2010):
%           unbiased, so its expected value is 0 for random phases whatever
%           the number of heartbeats (unlike ITC)
%
% Surrogate control (opts.nSurr > 0): one surrogate beat is drawn at a random
% position in each inter-beat interval of the whole heartbeat train
% (opts.allBeats), so the surrogate beats fall at independent random cardiac
% phases whatever the heart rate variability; as many surrogate beats as real
% ones are drawn. The null is the average of windows unrelated to the
% heartbeat, and the question answered is whether the activity at a given
% latency is locked to the heartbeats (one condition). Since the surrogate
% windows fall at all cardiac phases, the null mean at every latency is the
% average over the cardiac cycle (zero for high-passed data): a point differs
% from the surrogates when it differs from that average. A large
% heartbeat-locked deflection (the QRS, a late HEP) moves the rest of the
% cycle away from it with the opposite sign, so latencies without
% heartbeat-locked activity of their own can differ too; the family-wise
% control below holds when nothing is heartbeat-locked (see
% tests/sim_hep_fwer.m). (Surrogates made by shifting the whole train, as in
% Park et al., 2018, are still heartbeat-locked: they test contrasts between
% conditions, not one condition; removed in v1.7.)
% The HEP, HRSP and HRPC are recomputed for each surrogate (the HEP for at
% least 1000 surrogates, opts.nSurrHep: it is cheap, and its max statistics
% need many), and each point (channel x [frequency x] latency) is compared
% with its surrogate distribution at the same latency: z-score, p-value from
% the z-score, FDR-corrected across all points (Benjamini & Hochberg), and the
% empirical p-value (smallest possible: 1/(nSurr+1), 1/(nSurrHep+1) for the
% HEP). Family-wise error control over all points (channels x [frequencies
% x] latencies): the max statistic (p_tmax: share of surrogates whose
% largest |z| reaches the point's |z|), and with opts.stat 'tfce',
% threshold-free cluster enhancement along time (p_tfce; Smith & Nichols,
% 2009; E = 0.5, H = 2), which favors effects sustained over contiguous
% latencies over isolated peaks (detected a simulated response more often,
% tests/sim_hep_fwer.m, but took 28% of the runtime with 16 channels, 27
% frequencies and 100 surrogates: not computed by default). They do not
% depend on the shape of the null distribution: each
% surrogate is standardized by the mean and SD of the others, as the data
% are by all of them (no significant HRSP or HRPC point on white noise in
% tests/test_compute_hep_tf.m check 10). The FDR p-values come from
% assumed tails of the z-scores against the surrogates, a weaker basis
% with 100 surrogates and many points. HEP and HRSP: two-sided Gaussian tail. HRPC:
% one-sided exponential tail, exp(-(1+z)), since with random phases the
% pairwise phase consistency is a shifted exponential, not a Gaussian
% (Rayleigh; a Gaussian tail gave 2.5-4% FDR-significant points on pure
% noise).
% Note: this is a check of heartbeat locking. The cardiac field artifact is
% heartbeat-locked too, so it does not test a neural origin, and it does not
% compare conditions.
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
%           .tf         compute HRSP and HRPC (default true)
%           .nSurr      number of surrogates (default 0 = none)
%           .nSurrHep   number of surrogates for the HEP (default: nSurr, at
%                       least 1000; the first nSurr are those of HRSP/HRPC)
%           .stat       statistic the plots use: 'tmax' (default), 'tfce'
%                       or 'fdr' (tmax and fdr are always computed, tfce
%                       only when chosen)
%           .allBeats   all heartbeats of the recording (sample indices), from
%                       which the surrogates are built (default: beats)
%           .boundaries sample positions of data discontinuities (EEGLAB
%                       'boundary' events): windows crossing them are excluded
%           .exclude    [start end] sample spans of data rejected by the
%                       cleaning (e.g. the windows of the bad epochs): no
%                       surrogate window (with its wavelet support) overlaps
%                       them
%           .seed       random seed for the surrogate shifts (default 1)
%           .hep_times  time points of the HEP epochs in ms (default: every
%                       sample of win), e.g. HEP.times
%           .keepNull   also return the surrogate HEPs (surr.hep.null,
%                       channels x times x surrogates; default false)
%           .prestim    struct with .win (ms, e.g. [-300 -50]), .bands
%                       (bands x [low high) Hz, not overlapping) and .names:
%                       each heartbeat's induced power in each band, averaged
%                       over the window latencies and the band's frequencies
%                       (tf.prestim)
%
% Outputs:
%   tf   - struct: times (ms), freqs (Hz), power, hrsp, power_induced,
%          hrsp_induced, hrpc (channels x freqs x times), hep (channels x
%          times, the average of the same beats), nBeats, beats (the sample
%          indices used), cycles; with opts.prestim, prestim: win, names,
%          bands, power (channels x bands x beats, uV^2; NaN for a band with
%          no frequency computed)
%   surr - struct (empty without surrogates): nSurr, nSurrHep, and for hep,
%          hrsp, hrsp_induced, hrpc: z, p, p_fdr, p_emp (see above), null_mean, null_sd,
%          null_lo/null_hi (2.5th and 97.5th percentiles of the surrogates);
%          and family-wise: p_tmax and tmax_thresh (|z| to exceed for
%          p < .05), with stat 'tfce' also tfce (signed), p_tfce, tfce_thresh
%
% References:
%   Cohen & Donner (2013). Midfrontal conflict-related theta-band power
%   reflects neural oscillations that predict behavior. Journal of
%   Neurophysiology, 110(12), 2752-2763 (ERP subtracted from each trial
%   before the time-frequency decomposition).
%   Smith & Nichols (2009). Threshold-free cluster enhancement: addressing
%   problems of smoothing, threshold dependence and localisation in cluster
%   inference. NeuroImage, 44(1), 83-98.
%   Lee et al. (2024). Heartbeat-related spectral perturbation of
%   electroencephalogram reflects dynamic interoceptive attention states in
%   the trial-by-trial classification analysis. NeuroImage.
%   Park et al. (2018). Neural sources and underlying mechanisms of neural
%   responses to heartbeats, and their role in bodily self-consciousness: an
%   intracranial EEG study. Cerebral Cortex, 28(7), 2351-2364 (surrogate
%   R-peaks shifted by -500 to +500 ms, after Park et al., 2014).
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
nSurrHep = getdef(opts, 'nSurrHep', 1000*(nSurr > 0));   % HEP: cheap, and the max statistics need many
nSurrHep = max(nSurrHep, nSurr) * (nSurr > 0);
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

% Surrogate heartbeat trains (the same for all channels): nSurrHep trains,
% the first nSurr are used for HRSP and HRPC too (the same trains as with
% nSurr alone). One beat at a random position in each inter-beat interval
% of the whole train; nB valid surrogate beats drawn. When edges,
% discontinuities and excluded spans leave fewer than nB, more positions are
% drawn (another one per interval, up to 4 more rounds): with fewer beats
% than the HEP, the null would be wider.
if nSurr > 0
    rs = RandStream('twister','Seed',seed);
    surrBeats = cell(nSurrHep,1);
    allB = sort(round(getdef(opts, 'allBeats', beats)));
    allB = allB(:);
    ibi = diff(allB);
    excl = getdef(opts, 'exclude', zeros(0,2));
    lims = [hepIdx(1) hepIdx(end)] + [-half half];
    for s = 1:nSurrHep
        b = zeros(0,1);
        for r = 1:5
            c = allB(1:end-1) + floor(rand(rs, numel(ibi), 1) .* ibi);
            b = [b; c(okBeat(c) & ~any(c + lims(1) <= excl(:,2)' & c + lims(2) >= excl(:,1)', 2))]; %#ok<AGROW>
            if numel(b) >= nB, break, end
        end
        if numel(b) > nB
            b = b(randperm(rs, numel(b), nB));
        end
        b = sort(b);
        surrBeats{s} = b;
    end
end

% HEP: average of the same beats (and surrogates)
tf.hep = mean_epochs(X, beats, hepIdx);
tf.hep_times = hepIdx / fs * 1000;
surr = struct([]);
if nSurr > 0
    surr = struct('nSurr', nSurr, 'nSurrHep', nSurrHep, 'stat', getdef(opts, 'stat', 'tmax'));
    doTFCE = strcmpi(surr.stat, 'tfce');
    nullHep = zeros([size(tf.hep) nSurrHep]);
    for s = 1:nSurrHep
        nullHep(:,:,s) = mean_epochs(X, surrBeats{s}, hepIdx);
    end
    surr.hep = surrogate_stats(tf.hep, nullHep, 'both');
    surr.hep = fwer_stats(surr.hep, nullHep, doTFCE);
    surrBeats = surrBeats(1:nSurr);
    if getdef(opts, 'keepNull', false)
        surr.hep.null = nullHep;        % channels x times x surrogates
    end
    clear nullHep
end

tf.times = times; tf.freqs = freqs; tf.cycles = cycles; tf.nBeats = nB; tf.beats = beats;
if ~doTF
    return
end

% Pre-R-peak power of each heartbeat (induced, per band)
doPre = isfield(opts,'prestim') && ~isempty(opts.prestim);
if doPre
    pre = opts.prestim;
    tPre = tIdx(times >= pre.win(1) & times <= pre.win(2));
    nBand = size(pre.bands, 1);
    bandOf = zeros(1, nF);
    for iBand = 1:nBand
        in = freqs >= pre.bands(iBand,1) & freqs < pre.bands(iBand,2);
        if any(bandOf(in) > 0)
            error('compute_hep_tf: the pre-R-peak bands overlap.')
        end
        bandOf(in) = iBand;
    end
    nPerBand = accumarray(bandOf(bandOf > 0)', 1, [nBand 1])';
    prePow = zeros(nChan, nBand, nB);
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
tf.power_induced = tf.power;
tf.hrpc  = tf.power;
if nSurr > 0
    nullP = zeros(nChan,nF,numel(tIdx),nSurr,'single'); nullC = nullP; nullPi = nullP;
    % surrogate beats as columns (padded with their first beat), in blocks
    % of surrogates processed together
    nPer = cellfun(@numel, surrBeats);
    SB = zeros(max(nPer), nSurr);
    for s = 1:nSurr
        SB(:,s) = surrBeats{s}([1:nPer(s) ones(1, max(nPer)-nPer(s))]);
    end
    blk = max(1, min(nSurr, floor(2e6 / (size(SB,1)*numel(tIdx)))));
end
for iChan = 1:nChan
    FX = fft(X(iChan,:), nfft);
    for iF = 1:nF
        co = ifft(FX .* W(iF,:));
        co = co(1:nPts);
        [pw, pc, pwi] = tf_stats(co, beats, tIdx);
        tf.power(iChan,iF,:) = pw;
        tf.power_induced(iChan,iF,:) = pwi;
        tf.hrpc(iChan,iF,:) = pc;
        if doPre && bandOf(iF) > 0 && ~isempty(tPre)
            C = co(beats + tPre);                            % beats x latencies
            C = reshape(C, nB, []);
            prePow(iChan,bandOf(iF),:) = prePow(iChan,bandOf(iF),:) + ...
                reshape(mean(abs(C - mean(C,1)).^2, 2), 1, 1, []) / nPerBand(bandOf(iF));
        end
        if nSurr > 0
            for a = 1:blk:nSurr
                s = a:min(nSurr, a+blk-1);
                [pwS, pcS, pwiS] = tf_stats_many(co, SB(:,s), nPer(s), tIdx);
                nullP(iChan,iF,:,s) = reshape(pwS, 1, 1, [], numel(s));
                nullPi(iChan,iF,:,s) = reshape(pwiS, 1, 1, [], numel(s));
                nullC(iChan,iF,:,s) = reshape(pcS, 1, 1, [], numel(s));
            end
        end
    end
end
tf.hrsp = 10*log10(tf.power ./ mean(tf.power,3));
tf.hrsp_induced = 10*log10(tf.power_induced ./ mean(tf.power_induced,3));
if doPre
    prePow(:, nPerBand == 0 | isempty(tPre), :) = NaN;
    tf.prestim = struct('win', pre.win, 'names', {pre.names}, 'bands', pre.bands, 'power', prePow);
end

if nSurr > 0
    % HRSP of each surrogate, relative to its own window mean
    nullH = 10*log10(double(nullP) ./ mean(double(nullP),3));
    clear nullP
    surr.hrsp = fwer_tf(surrogate_stats(tf.hrsp, nullH, 'both'), nullH, doTFCE);
    nullH = 10*log10(double(nullPi) ./ mean(double(nullPi),3));
    clear nullPi
    surr.hrsp_induced = fwer_tf(surrogate_stats(tf.hrsp_induced, nullH, 'both'), nullH, doTFCE);
    clear nullH
    nullC = double(nullC);
    surr.hrpc = fwer_tf(surrogate_stats(tf.hrpc, nullC, 'exponential'), nullC, doTFCE);
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

function [pw, pc, pwi] = tf_stats(co, beats, tIdx)
% Mean power, pairwise phase consistency and induced power (the variance of
% the coefficients around their mean, unbiased) across beats (1 x times)
N = numel(beats);
C = reshape(co(beats + tIdx), N, numel(tIdx));   % beats x times
pw = mean(abs(C).^2, 1);
pwi = max(0, N/(N-1) * (pw - abs(mean(C,1)).^2));
z = C ./ abs(C);
pc = (abs(sum(z,1)).^2 - N) / (N*(N-1));
end

function [pw, pc, pwi] = tf_stats_many(co, B, n, tIdx)
% TF_STATS of several sets of beats at once: B (beats x sets, each column
% padded after its n(s) beats), outputs times x sets
[K, nS] = size(B);
w = (1:K)' <= n(:)';                                     % beats x sets: real beats
C = reshape(co(B(:) + tIdx), K, nS, numel(tIdx));       % beats x sets x times
z = (C ./ abs(C)) .* w;
C = C .* w;
n = n(:)';
pw = reshape(sum(abs(C).^2, 1), nS, []).' ./ n;
pwi = max(0, n./(n-1) .* (pw - abs(reshape(sum(C, 1), nS, []).' ./ n).^2));
pc = (abs(reshape(sum(z, 1), nS, []).').^2 - n) ./ (n.*(n-1));
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
    S.p = min(1, exp(-(1 + S.z)));      % z + 1 ~ Exp(1) under the null
    S.p_emp = (sum(null >= real, nd) + 1) / (nSurr + 1);
end
S.p_emp(isnan(S.z)) = NaN;              % e.g. a flat channel
S.p_fdr = fdr_bh(S.p);
end

function S = fwer_tf(S, null, doTFCE)
% Family-wise statistics of a time-frequency measure (FWER_STATS with the
% channels x frequencies as rows, TFCE along time)
sz = size(S.z);
nT = sz(end);
R = struct('z', reshape(S.z, [], nT), 'null_mean', reshape(S.null_mean, [], nT), ...
    'null_sd', reshape(S.null_sd, [], nT));
R = fwer_stats(R, reshape(null, [], nT, size(null, ndims(null))), doTFCE);
S.p_tmax = reshape(R.p_tmax, sz); S.tmax_thresh = R.tmax_thresh;
if doTFCE
    S.tfce = reshape(R.tfce, sz); S.p_tfce = reshape(R.p_tfce, sz); S.tfce_thresh = R.tfce_thresh;
end
end

function S = fwer_stats(S, null, doTFCE)
% Family-wise error control over all points (channels x latencies) of the
% HEP, from the surrogates: max statistic (tmax) and, if doTFCE,
% threshold-free cluster enhancement along time (TFCE; Smith & Nichols,
% 2009; E = 0.5, H = 2, dh = 0.1 in z units). Each surrogate's z-scores
% against the surrogate mean and SD give the null distribution of the
% largest |z| and of the largest |TFCE| over all points; a point's p-value
% is the share of surrogates whose maximum reaches its value. TFCE
% integrates the z-score over the extent of the run of contiguous latencies
% above each height, so sustained effects gain and isolated peaks lose.
nS = size(null, 3);
z = S.z; z(~isfinite(z)) = 0;
if doTFCE, S.tfce = tfce_time(z); end
maxZ = zeros(nS,1); maxT = zeros(nS,1);
chunk = max(1, floor(5e6 / numel(z)));          % surrogates per block (memory)
SS = (nS-1) * S.null_sd.^2;                     % sum of squared deviations
for a = 1:chunk:nS
    b = min(nS, a+chunk-1);
    % each surrogate standardized by the mean and SD of the others, as the
    % data are by all of them (with itself included, an extreme surrogate
    % value inflates its own SD: its z is shrunk and the max statistic is
    % liberal, 5 of 9 HRSP/HRPC noise tests significant with 100 surrogates)
    d = null(:,:,a:b) - S.null_mean;
    zn = (d * nS/(nS-1)) ./ sqrt(max(SS - d.^2 * nS/(nS-1), 0) / (nS-2));
    zn(~isfinite(zn)) = 0;
    maxZ(a:b) = squeeze(max(abs(zn), [], [1 2]));
    if doTFCE
        zn = permute(zn, [1 3 2]);               % (channels x surrogates) x times
        Tn = tfce_time(reshape(zn, [], size(zn,3)));
        maxT(a:b) = max(reshape(abs(Tn), size(zn,1), b-a+1, []), [], [1 3]);
    end
end
S.p_tmax = (sum(maxZ(:)' >= reshape(abs(z), [], 1), 2) + 1) / (nS + 1);
S.p_tmax = reshape(S.p_tmax, size(z));
% values to exceed for p < .05: the k-th largest surrogate maximum, with k
% the largest count of surrogates at or above a value that keeps p < .05
k = ceil(.05*(nS+1)) - 1;
mz = sort(maxZ, 'descend');
S.tmax_thresh = NaN;
if k >= 1, S.tmax_thresh = mz(k); end
if doTFCE
    S.p_tfce = (sum(maxT(:)' >= reshape(abs(S.tfce), [], 1), 2) + 1) / (nS + 1);
    S.p_tfce = reshape(S.p_tfce, size(z));
    mt = sort(maxT, 'descend');
    S.tfce_thresh = NaN;
    if k >= 1, S.tfce_thresh = mt(k); end
end
end

function T = tfce_time(Z)
% TFCE of each row of Z (series x times) along time, signed: positive and
% negative z-scores are enhanced separately. Each height is computed only
% for the series that reach it (sorted by their maximum).
E = 0.5; H = 2; dh = 0.1;
Zt = Z.';                                       % times x series
T = zeros(size(Zt));
for sgn = [1 -1]
    Zs = max(sgn * Zt, 0);
    [top, ord] = sort(max(Zs, [], 1), 'descend');
    if isempty(top), continue, end
    Zs = Zs(:, ord);
    Ts = zeros(size(Zs));
    for h = dh:dh:top(1)
        a = find(top >= h, 1, 'last');          % series 1:a reach h
        M = Zs(:, 1:a) >= h;
        Ts(:, 1:a) = Ts(:, 1:a) + (h^H * dh) .* (run_extent(M) .^ E) .* M;
    end
    T(:, ord) = T(:, ord) + sgn * Ts;
end
T = T.';
end

function L = run_extent(M)
% Length of the run of true values (down each column) that each element of
% M belongs to; 0 where M is false
c = cumsum(M, 1);
left = c - cummax(c .* ~M, 1);                  % true values up to here in the run
c = cumsum(M, 1, 'reverse');
right = c - cummax(c .* ~M, 1, 'reverse');      % from here to the run's end
L = (left + right - 1) .* M;
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
