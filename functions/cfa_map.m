% CFA_MAP - Scalp map of the cardiac field artifact at the QRS, and whether
% it stayed the same over the recording.
%
% The map is taken at the largest global field power (GFP) between -30 and
% 60 ms of the average heartbeat (20% trimmed mean, no baseline correction:
% the data must be high-pass filtered; a warning when the channels' offsets
% exceed the waveform, as in unfiltered data). Its shape is set by the orientation of the heart's
% electrical axis relative to the electrodes: it changes with posture (the
% heart's axis rotates between lying, sitting and standing), with the
% orientation of the head relative to the torso (turning or nodding; Dirlich
% et al., 1997), with the electrode positions and the reference, and between
% people (body build). A map that changes between conditions means that the
% cardiac field artifact left in the data differs between them.
%
% Stability over the recording (maps of the heartbeats averaged over 3
% samples around the QRS latency, each epoch's mean before -50 ms removed),
% two comparisons:
%   - before vs after: the map of the heartbeats before each time point
%     (20 to 80% of the heartbeats, in steps of 5%) correlated with the map
%     of the heartbeats after it: a lasting change (e.g. the head turned and
%     stayed turned) lowers the correlation at the time of the change;
%   - each block vs the rest: consecutive blocks (default 5) correlated
%     with the map of all other heartbeats: a temporary change.
% Noise alone lowers these correlations, so each is compared with the same
% correlation for the heartbeats in random order (default 1000
% permutations), as a z-score against them; p is the share of permutations
% whose lowest z-score is as low as the observed one (Bonferroni-corrected
% for the two comparisons). At that lowest point, rho = r / median r in
% random order estimates the correlation of the noise-free maps (1 = no
% change; NaN when the random-order correlation is below .5). The test is
% not testable when two random halves of the heartbeats correlate below .5
% (maps too noisy). A change is reported when p < .01 and rho < .95. It
% shows that the map at the QRS latency changed, which a moved head, body or
% electrodes do, and also a moved time-locking (a PPG pulse arrival time
% that drifts by ~20 ms or more, R-peaks detected on another wave).
% Sensitivity (tests/sim_cfa_drift.m; 12 resting recordings of 120-370
% heartbeats, raw and cleaned, the cardiac field rotated relative to the
% head for half of the heartbeats by spherical-spline interpolation): no
% false alarm on the 24 unchanged ones; a 30-deg head turn (rotated map
% correlating .51-.82 with the original) was detected in 3-4/12, 20 deg
% (.67-.91) in 1-2/12, 10 deg (.90-.97) or a nod of 30 deg in 0-1/12, a
% 30-deg turn lasting 30% of the recording in 1/12.
% A stable result does not exclude movements: only large lasting changes
% are seen. (Robust means and the whole QRS complex did not do better;
% p < .05 gave false alarms, 1/12.)
%
% Usage:
%   cfa = cfa_map(data, times, beatTime)
%   cfa = cfa_map(data, times, beatTime, opts)
%
% Inputs:
%   data     - EEG epochs (channels x times x heartbeats, in recording order),
%              uV, time-locked to the R-peaks
%   times    - epoch latencies (ms)
%   beatTime - time of each heartbeat in the recording (s)
%   opts     - optional struct: nBlocks (5), nPerm (1000), seed (1),
%              win ([-30 60] ms, where the QRS map is searched),
%              info (cell of text lines added to the notes, e.g. the
%              reference and the heart-artifact removal), unit ('uV')
%
% Output:
%   cfa - struct: latency (ms), map (channels x 1, uV), gfp (uV),
%         split (before vs after: minutes, r, r_random (median), r_random_lo
%         (2.5th percentile), rho, at each time point), blocks (each block vs
%         the rest: beats, minutes, r, r_random, r_random_lo, rho), p_change,
%         when (minutes of the largest change), rho (there), testable, changed,
%         summary (one line), notes (cell of text lines: values and how to
%         use them)
%
% Reference:
%   Dirlich, Vogl, Plaschke & Strian (1997). Cardiac field effects on the
%   EEG. Electroencephalography and Clinical Neurophysiology, 102(4), 307-315.
%
% Copyright (C) - Cedric Cannard, 2026

function cfa = cfa_map(data, times, beatTime, opts)

if nargin < 4, opts = struct; end
nBlocks = getdef(opts, 'nBlocks', 5);
nPerm   = getdef(opts, 'nPerm', 1000);
seed    = getdef(opts, 'seed', 1);
win     = getdef(opts, 'win', [-30 60]);
info    = getdef(opts, 'info', {});
unit    = getdef(opts, 'unit', 'uV');

[nChan, ~, nB] = size(data);
pre = times < -50;                                   % before the QRS onset
X = double(data);
if any(pre), X = X - mean(X(:,pre,:), 2); end        % each epoch's offset removed (change test)
% Reported map: from the epochs as they are (no baseline, like the HEP
% plots). Unfiltered data (the average heartbeat's mean exceeds its range
% on most channels) give the offsets' map: BRAINBEATS_PROCESS stops before.
erp = trimmean(double(data), 20, 3);
if median(abs(mean(erp,2))) > median(max(erp,[],2) - min(erp,[],2))
    warning('BrainBeats:unfiltered', ['cfa_map: the channels'' offsets exceed the heartbeat waveform: the data are ' ...
        'not high-pass filtered, and the QRS map shows the offsets. Filter the data first.'])
end
iQ = find(times >= win(1) & times <= win(2));
if isempty(iQ) || nChan < 3 || nB < 20
    error('cfa_map: no latency in [%g %g] ms, fewer than 3 channels, or fewer than 20 heartbeats.', win)
end
[~, k] = max(std(erp(:,iQ), 0, 1));
kq = iQ(k);
cfa.latency = times(kq);
cfa.map = erp(:,kq);
cfa.gfp = std(cfa.map);

% Map of each heartbeat around the QRS latency; split points and blocks
Y = squeeze(mean(X(:, max(1,kq-1):min(numel(times),kq+1), :), 2));   % channels x beats
ks = unique(round((.2:.05:.8) * nB));                % last beat before each split
nBlocks = max(2, min(nBlocks, floor(nB/10)));
edges = round(linspace(0, nB, nBlocks+1));
B = zeros(nB, nBlocks);
for b = 1:nBlocks, B(edges(b)+1:edges(b+1), b) = 1; end

[rS, rB] = map_r(Y, ks, B);
rs = RandStream('twister', 'Seed', seed);
nS = zeros(nPerm, numel(ks)); nBk = zeros(nPerm, nBlocks);
for p = 1:nPerm
    [nS(p,:), nBk(p,:)] = map_r(Y(:, randperm(rs, nB)), ks, B);
end
[zS, zSn] = zscore_null(rS, nS);
[zB, zBn] = zscore_null(rB, nBk);
pS = perm_p(zS, zSn); pB = perm_p(zB, zBn);
cfa.p_change = min(1, 2*min(pS, pB));

bt = beatTime(:)' / 60;
mS = median(nS, 1); mB = median(nBk, 1);
cfa.split = struct('minutes', bt(ks), 'r', rS, 'r_random', mS, 'r_random_lo', prctile(nS, 2.5, 1), 'rho', rS ./ mS);
for b = 1:nBlocks
    in = edges(b)+1:edges(b+1);
    cfa.blocks(b) = struct('beats', [in(1) in(end)], 'minutes', bt([in(1) in(end)]), 'r', rB(b), ...
        'r_random', mB(b), 'r_random_lo', prctile(nBk(:,b), 2.5), 'rho', rB(b) / mB(b));
end
% Largest change: the point with the lowest z-score of the comparison with
% the lower p; rho there (NaN if its random-order correlation is below .5)
cfa.testable = max(mS) >= .5;
if pS <= pB
    [~, i] = min(zS); cfa.when = bt(ks(i)); r0 = rS(i); m0 = mS(i);
else
    [~, i] = min(zB); cfa.when = cfa.blocks(i).minutes; r0 = rB(i); m0 = mB(i);
end
cfa.rho = NaN;
if m0 >= .5, cfa.rho = r0 / m0; end
cfa.changed = cfa.testable && cfa.p_change < .01 && cfa.rho < .95;

% Notes
when = strjoin(arrayfun(@(m) sprintf('%.1f', m), cfa.when, 'UniformOutput', false), '-');
if ~cfa.testable
    cfa.summary = sprintf(['Not testable: the QRS map of half the heartbeats is too noisy (correlation of two random ' ...
        'halves %.2f, < .5). A weak QRS map is common after cleaning or heart-artifact removal.'], max(mS));
elseif cfa.changed
    cfa.summary = sprintf(['CHANGED at %s min (p = %.3f; noise-corrected map correlation there %.2f): the head or ' ...
        'body position or the electrodes may have moved, or the time-locking (e.g. a PPG pulse arrival time ' ...
        'that drifts, or R-peaks detected on another wave).'], when, cfa.p_change, cfa.rho);
else
    cfa.summary = sprintf('Stable (p = %.3f).', cfa.p_change);
end
L = {sprintf('QRS map (cardiac field artifact): %g ms after the R-peak, GFP %.2f %s, %d channels, %d heartbeats.', ...
    cfa.latency, cfa.gfp, unit, nChan, nB)};
L = [L info(:)' {cfa.summary} {''} ...
    {'Map of the heartbeats before vs after each time point (r; random order: the same with the heartbeats shuffled, what noise alone gives):'}];
for i = 1:numel(ks)
    L{end+1} = sprintf('  %5.1f min: r = %.3f, random order %.3f (2.5th percentile %.3f)', bt(ks(i)), rS(i), mS(i), cfa.split.r_random_lo(i)); %#ok<AGROW>
end
L{end+1} = 'Map of each block vs all other heartbeats:';
for b = 1:nBlocks
    L{end+1} = sprintf('  %.1f-%.1f min: r = %.3f, random order %.3f (2.5th percentile %.3f)', ...
        cfa.blocks(b).minutes, rB(b), mB(b), cfa.blocks(b).r_random_lo); %#ok<AGROW>
end
L = [L {'' 'How to use it:' ...
    '- Compare conditions or sessions on this map (spatial correlation of their maps). Different maps mean that the cardiac field artifact left in the data differs, which alone can produce HEP differences, mostly at latencies where the HEP maps resemble the QRS map (r_CFA in the HEP figure).' ...
    '- The map follows the orientation of the heart''s electrical axis relative to the electrodes: posture (lying, sitting, standing, lying on one side), the orientation of the head relative to the torso (turning, nodding: Dirlich et al., 1997), electrode positions and the reference, and body build.' ...
    '- Record all conditions in the same posture, with the head supported in the same position. Posture also changes the HEP itself (larger lying down than upright; Plantard et al., 2025).' ...
    '- The stability test only sees large lasting changes: in a simulation on 12 recordings of 120-370 heartbeats, a 30-degree head turn for half of the recording was detected in 3-4 of 12, a 20-degree turn in 1-2, a 10-degree turn in 0-1, with no false alarm. A stable result does not exclude movements. A change can also come from the time-locking: a PPG pulse arrival time that drifts by ~20 ms or more, or R-peaks detected on another wave.' ...
    '- References: Dirlich, Vogl, Plaschke & Strian (1997), Electroencephalography and Clinical Neurophysiology 102:307-315. Plantard et al. (2025), Journal of Neurophysiology 134:830-842.'}];
cfa.notes = L;
end

function [rS, rB] = map_r(Y, ks, B)
% Correlation of the mean map before vs after each split point, and of each
% block's mean map vs the mean map of the other beats
n = size(Y, 2);
C = cumsum(Y, 2);
before = C(:,ks) ./ ks;
after = (C(:,end) - C(:,ks)) ./ (n - ks);
rS = colcorr(before, after);
nb = sum(B, 1);
Bs = Y * B;
rB = colcorr(Bs ./ nb, (C(:,end) - Bs) ./ (n - nb));
end

function r = colcorr(A, B)
A = A - mean(A,1); B = B - mean(B,1);
r = sum(A.*B,1) ./ sqrt(sum(A.^2,1) .* sum(B.^2,1));
end

function p = perm_p(z, zn)
% Share of permutations whose lowest z-score is as low as the observed one
% (1 if nothing could be computed)
if all(isnan(z)), p = 1; return, end
p = (sum(min(zn, [], 2) <= min(z)) + 1) / (size(zn,1) + 1);
end

function [z, zn] = zscore_null(r, rn)
% z-scores of the observed and permuted correlations against the permuted ones
m = median(rn, 1); s = std(rn, 0, 1);
s(s == 0) = eps;
z = (r - m) ./ s;
zn = (rn - m) ./ s;
end

function v = getdef(s, f, d)
if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
