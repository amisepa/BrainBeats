function [nn_intervals, nn_t, nPeaks, idx_bad, idx_interp] = clean_rr(rr_t, rr_intervals, peak_amp, rPeaks, varargin)
% CLEAN_RR - Clean an RR interval series into NN intervals and fill missed beats.
%
% Usage:
%   [nn, nn_t, nPeaks, idx_bad, idx_interp] = clean_rr(rr_t, rr, peak_amp, rPeaks)
%   [...] = clean_rr(..., 'ecg_signal', sig, 'fs', fs)
%
% Inputs (one value per RR interval, as returned by get_RR):
%   rr_t         - time of the peak closing each interval (s)
%   rr_intervals - RR intervals (s)
%   peak_amp     - signal amplitude at each peak ([] to skip the amplitude check)
%   rPeaks       - peak sample indices
% Optional name-value pairs:
%   'interpolate_missing' - fill gaps with synthetic beats (default true, logical)
%   'max_missing'         - max beats inserted per gap (default Inf)
%   'mad_threshold'       - MAD multiplier for the lower RR limit (default 10)
%   'amp_mad_threshold'   - MAD multiplier for amplitude outliers (default 15;
%                           Inf disables)
%   'win_size'            - local-median window (s, default 15)
%   'gap_ratio'           - gap if interval > gap_ratio x local median (default 1.5)
%   'ectopic_rule'        - 'adaptive' (default since v1.7): ectopic beats of
%                           the beat classification of Lipponen & Tarvainen
%                           (2019), thresholds scaled to the beat-to-beat
%                           variability of the recording (see Notes);
%                           'ratio': fixed rule ('ectopic_ratio', the v1.6
%                           default); 'none'
%   'ectopic_ratio'       - 'ratio' rule: ectopic beat if its interval is
%                           shorter than (1 - ratio) x the local median and
%                           the next one longer than (1 + ratio) x (default
%                           0.2); Inf disables the ectopic step (both rules)
%   'ectopic_min_shift'   - 'adaptive' rule: a beat is replaced only if that
%                           moves it by at least this fraction of the mean
%                           of its two intervals (default 0.1; 0 replaces
%                           every beat classified ectopic, as published)
%   'ecg_signal'          - ECG/PPG signal, to place synthetic beats on samples
%   'fs'                  - its sampling rate (Hz)
%   'sig_t'               - its time axis in s (default (0:N-1)/fs)
%
% Outputs (column vectors):
%   nn_intervals - cleaned NN intervals (s)
%   nn_t         - their timestamps (s)
%   nPeaks       - peak sample indices; synthetic beats get the nearest
%                  sample of ecg_signal, or NaN if no signal is given
%   idx_bad      - true for gaps that remain unfilled
%   idx_interp   - true for synthetic (inserted) intervals
%
% Steps:
%   1. Limits from the data: min RR = median - mad_threshold x MAD,
%      bounded to [0.30 2.00] s (200-30 bpm).
%   2. Beats with outlying amplitude are removed.
%   3. Intervals shorter than min RR are removed; a removed interval is
%      added to the next kept one so NN still spans two kept beats.
%   3b. Ectopic beats (premature beat then compensatory pause) are removed
%      the same way. 'ratio': a short-long pair against the local median;
%      the merged interval (~2 beats) is filled by step 5 when it is a gap.
%      'adaptive': beats classified ectopic from their interval difference
%      (dRR) and the neighboring ones, divided by 5.2 x the quartile
%      deviation of |dRR| over 91 beats (Lipponen & Tarvainen, 2019; as in
%      NeuroKit2's 'Kubios' method), one at a time, the series classified
%      again after each; the merged interval always gets one beat in its
%      middle (step 5). Without this step, a pair at 0.75/1.25 x passed
%      steps 3-4.
%   4. Gaps: intervals > gap_ratio x the median of neighbors within win_size,
%      and the merged intervals of 'adaptive'.
%   5. Each gap is split into equal intervals by inserting
%      round(gap / local median) - 1 synthetic beats (at least one for the
%      merged intervals of 'adaptive').
%
% Notes:
%   amp_mad_threshold is deliberately lenient (15 MAD, not isoutlier's 3):
%   ECG amplitude varies with respiration, posture and electrode contact;
%   only disconnections or rail artifacts should be caught here. Long
%   intervals are never removed, only flagged or filled (steps 4-5).
%
%   Ectopic rules on 40 ds003838 and sample recordings with premature beats
%   injected (full compensatory pause; tests/compare_ectopic_rules.m,
%   2026-11-05): 'ratio' detected 83-89% of them at 0.70 x the interval,
%   60-65% at 0.75 x and 13-15% at 0.80 x, 'adaptive' 98-99%, 96-97% and
%   92-96% (RMSSD error 0.6-6.9 vs 5.7-22.7 ms). On the recordings as they
%   are (12,800 beats), 'adaptive' replaced 5 beats 'ratio' kept: 2 false
%   detections (an electrode pop, a spike), 1 sinus beat and 2 unclear on
%   the ECG; 'ratio' replaced 1 (a false detection, replaced by 'adaptive'
%   too). With 'ectopic_min_shift' 0 (the published rule), 'adaptive' also
%   replaced 8 sinus beats, each by less than 10% of its interval.
%
% Reference:
%   Lipponen, J. A., & Tarvainen, M. P. (2019). A robust algorithm for heart
%   rate variability time series artefact correction using novel beat
%   classification. Journal of Medical Engineering & Technology, 43(3),
%   173-181.
%
% Copyright (C) - Cedric Cannard, 2023

% Parse inputs
p = inputParser;
addParameter(p, 'interpolate_missing', true,  @islogical);
addParameter(p, 'max_missing',         Inf,   @isnumeric);
addParameter(p, 'mad_threshold',       10,     @isnumeric);
addParameter(p, 'amp_mad_threshold',   15,    @isnumeric);
addParameter(p, 'win_size',            15,    @isnumeric);
addParameter(p, 'gap_ratio',           1.5,   @isnumeric);
addParameter(p, 'ectopic_rule',        'adaptive', @(x) any(strcmpi(x, {'ratio','adaptive','none'})));
addParameter(p, 'ectopic_ratio',       0.2,   @isnumeric);
addParameter(p, 'ectopic_min_shift',   0.1,   @isnumeric);
addParameter(p, 'ecg_signal',          [],    @isnumeric);
addParameter(p, 'fs',                  [],    @isnumeric);
addParameter(p, 'sig_t',               [],    @isnumeric);
parse(p, varargin{:});
interpolate_missing = p.Results.interpolate_missing;
max_missing         = p.Results.max_missing;
mad_threshold       = p.Results.mad_threshold;
amp_mad_threshold   = p.Results.amp_mad_threshold;
win_size            = p.Results.win_size;
gap_ratio           = p.Results.gap_ratio;
ectopic_rule        = lower(p.Results.ectopic_rule);
ectopic_ratio       = p.Results.ectopic_ratio;
ectopic_min_shift   = p.Results.ectopic_min_shift;
if ~isfinite(ectopic_ratio), ectopic_rule = 'none'; end
ecg_signal          = p.Results.ecg_signal(:);
fs                  = p.Results.fs;
sig_t               = p.Results.sig_t(:);

% Build signal time axis if not provided
use_ecg = ~isempty(ecg_signal) && ~isempty(fs);
if use_ecg && isempty(sig_t)
    sig_t = (0:length(ecg_signal)-1)' / fs;
end
if ~isempty(ecg_signal) && isempty(fs)
    warning('clean_rr: ecg_signal provided but fs missing. Synthetic nPeaks will be NaN.');
    use_ecg = false;
end

% Initialize
n_orig       = length(rr_intervals);
nn_intervals = rr_intervals(:);
nn_t         = rr_t(:);
nPeaks       = double(rPeaks(:));
idx_bad      = false(n_orig, 1);
idx_interp   = false(n_orig, 1);

% Step 1: Data-driven RR limits (max_rr is only reported; long intervals
% are handled as gaps in Step 4)
med_rr = median(nn_intervals, 'omitnan');
mad_rr = mad(nn_intervals, 1);
min_rr = med_rr - mad_threshold * mad_rr;
max_rr = med_rr + mad_threshold * mad_rr;

% Bound the limits to human extremes (30-200 bpm)
min_rr = max(min_rr, 0.30);   % 200 bpm
max_rr = min(max_rr, 2.00);   % 30 bpm
fprintf('\n--- clean_rr ---\n');
fprintf('  Median RR:    %.3f s (%.0f bpm)\n', med_rr, 60/med_rr);
fprintf('  MAD:          %.3f s\n', mad_rr);
fprintf('  Min RR (data-derived): %.3f s (%.0f bpm)\n', min_rr, 60/min_rr);
fprintf('  Max RR (data-derived): %.3f s (%.0f bpm)\n', max_rr, 60/max_rr);

% Step 2: Remove beats with abnormal amplitude
% Lenient MAD threshold (default 15x) to catch only disconnection artifacts
% (flat line, signal rail). Amplitude varies normally with respiration and
% posture; the interval-based Steps 3-5 handle ectopic and missed beats.
if ~isempty(peak_amp) && length(peak_amp) == length(nn_intervals) && isfinite(amp_mad_threshold)
    amp_med      = median(peak_amp(:), 'omitnan');
    amp_mad      = mad(peak_amp(:), 1);
    amp_outliers = abs(peak_amp(:) - amp_med) > amp_mad_threshold * amp_mad;
    if any(amp_outliers)
        fprintf('  Removing %d beats with abnormal amplitude (>%.0f MAD from median)\n', ...
            sum(amp_outliers), amp_mad_threshold);
    end
else
    amp_outliers = false(size(nn_intervals));
end

% Step 3: Remove intervals shorter than min_rr (and amplitude outliers)
too_short = nn_intervals < min_rr;
if any(too_short)
    fprintf('  Removing %d intervals below min RR (%.3f s)\n', sum(too_short), min_rr);
end
to_remove = amp_outliers | too_short;
n_removed = sum(to_remove);
% A removed beat ends its interval: fold that interval into the next kept one,
% so each NN interval still spans the time between two kept beats (a false
% beat splitting 0.8 s into 0.35 + 0.45 s gives back 0.8 s, not 0.45 s).
% Step 4 then flags the merged interval as a gap if a real beat was removed.
carry = 0;
for i = 1:length(nn_intervals)
    if to_remove(i)
        carry = carry + nn_intervals(i);
    elseif carry > 0
        nn_intervals(i) = nn_intervals(i) + carry;
        carry = 0;
    end
end
nn_intervals(to_remove) = [];
nn_t(to_remove)         = [];
nPeaks(to_remove)       = [];
idx_bad(to_remove)      = [];
idx_interp(to_remove)   = [];
if isempty(nn_intervals)
    warning('All intervals removed - check your data quality');
    return;
end

% Step 3b: Ectopic beats. A removed beat's interval is merged with the next
% one, as in Step 3; to_fill marks the intervals Step 5 must fill.
to_fill = false(size(nn_intervals));
switch ectopic_rule
    case 'ratio'
        % a premature beat (its interval shorter than (1 - ectopic_ratio) x
        % the local median) followed by a compensatory pause (longer than
        % (1 + ectopic_ratio) x): the merged interval is filled if it is a gap
        local_med = local_median(nn_intervals, nn_t, win_size, med_rr);
        ect = find(nn_intervals(1:end-1) < (1-ectopic_ratio)*local_med(1:end-1) & ...
            nn_intervals(2:end) > (1+ectopic_ratio)*local_med(1:end-1));
        if ~isempty(ect)
            fprintf('  Removing %d ectopic beats (short-long pair, %.0f%% from the local median)\n', ...
                numel(ect), 100*ectopic_ratio);
        end
    case 'adaptive'
        % ectopic beats of the classification of Lipponen & Tarvainen
        % (2019); each is replaced in the middle of its two intervals
        ect = find(adaptive_ectopic(nn_intervals, ectopic_min_shift));
        to_fill(ect+1) = true;
        if ~isempty(ect)
            fprintf('  Removing %d ectopic beats (Lipponen & Tarvainen, 2019)\n', numel(ect));
        end
    otherwise
        ect = [];
end
if ~isempty(ect)
    % flagged beats are never adjacent, so each merge is independent
    nn_intervals(ect+1) = nn_intervals(ect+1) + nn_intervals(ect);
    nn_intervals(ect) = [];
    nn_t(ect)         = [];
    nPeaks(ect)       = [];
    idx_bad(ect)      = [];
    idx_interp(ect)   = [];
    to_fill(ect)      = [];
    n_removed = n_removed + numel(ect);
end

% Step 4: Detect gaps against the local median of neighboring intervals
% (+/- win_size/2 s, excluding the interval itself)
local_med = local_median(nn_intervals, nn_t, win_size, med_rr);

idx_bad = nn_intervals > gap_ratio * local_med | to_fill;

if any(idx_bad)
    fprintf('  Flagging %d gaps (interval > %.1fx local median)\n', sum(idx_bad), gap_ratio);
else
    fprintf('  No gaps detected\n');
end

% Step 5 (optional): Fill gaps with synthetic beats, equally spaced
if interpolate_missing && any(idx_bad)

    gap_indices    = find(idx_bad);
    n_interp_total = 0;

    % last gap first, so earlier indices stay valid after insertion; the
    % first interval is skipped (no preceding beat time)
    for k = length(gap_indices):-1:1
        i = gap_indices(k);
        if i < 2, continue; end
        gap       = nn_intervals(i);
        ref_rr_i  = local_med(i);
        n_missing = round(gap / ref_rr_i) - 1;
        if to_fill(i), n_missing = max(n_missing, 1); end
        if n_missing < 1, continue; end
        if n_missing > max_missing
            fprintf('  Gap at t=%.2f s: %d beats missing, exceeds max_missing=%d - skipping\n', ...
                nn_t(i), n_missing, max_missing);
            continue;
        end

        j_vec = (1:n_missing)';
        new_t = nn_t(i-1) + j_vec * (gap / (n_missing + 1));
        sub_iv  = gap / (n_missing + 1);
        sub_ivs = repmat(sub_iv, n_missing + 1, 1);
        if use_ecg
            new_peaks = arrayfun(@(t) find_nearest_sample(sig_t, t), new_t);
        else
            new_peaks = NaN(n_missing, 1);
        end

        % the gap becomes the last sub-interval; the others are inserted before it
        nn_intervals(i) = sub_ivs(end);
        idx_bad(i)      = false;
        nn_intervals = [nn_intervals(1:i-1); sub_ivs(1:end-1);              nn_intervals(i:end)];
        nn_t         = [nn_t(1:i-1);         new_t;                         nn_t(i:end)];
        nPeaks       = [nPeaks(1:i-1);        new_peaks;                    nPeaks(i:end)];
        idx_bad      = [idx_bad(1:i-1);       false(n_missing,1);           idx_bad(i:end)];
        idx_interp   = [idx_interp(1:i-1);    true(n_missing,1);            idx_interp(i:end)];
        local_med    = [local_med(1:i-1);     repmat(ref_rr_i,n_missing,1); local_med(i:end)];

        n_interp_total = n_interp_total + n_missing;
        fprintf('  Gap at t=%.2f s (%.3f s): inserted %d beat(s) at %.3f s intervals\n', ...
            nn_t(i + n_missing), gap, n_missing, sub_iv);
    end

    if n_interp_total > 0
        fprintf('  Total synthetic beats inserted: %d\n', n_interp_total);
        if ~use_ecg
            fprintf('  Pass ''ecg_signal'' and ''fs'' for accurate peak placement.\n');
        end
    else
        fprintf('  No gaps were suitable for filling\n');
    end
end

% Summary
n_final = length(nn_intervals);
fprintf('  Removed:   %d/%d (%.1f%%)\n', n_removed,                   n_orig,  100*n_removed/n_orig);
fprintf('  Flagged:   %d/%d (%.1f%%)\n', sum(idx_bad),                n_final, 100*sum(idx_bad)/n_final);
fprintf('  Synthetic: %d/%d (%.1f%%)\n', sum(idx_interp),             n_final, 100*sum(idx_interp)/n_final);
fprintf('  Clean:     %d/%d (%.1f%%)\n', sum(~idx_bad & ~idx_interp), n_final, 100*sum(~idx_bad & ~idx_interp)/n_final);
fprintf('----------------\n');

end

%% Helpers
function idx = find_nearest_sample(sig_t, t)
[~, idx] = min(abs(sig_t - t));
end

function moved = adaptive_ectopic(rr, min_shift)
% Beats to replace (true at the interval each one closes), one at a time:
% classify, move the first ectopic beat to the middle of its neighbors,
% classify again. A flag at beat j marks interval j; the beat out of place
% is j or the one before (a premature beat often flags the beat closing its
% compensatory pause): the one whose move leaves the smaller successive
% differences over intervals j-2 to j+2. A beat is moved only if that
% shifts it by at least min_shift x the mean of its two intervals. Moved
% beats and their neighbors are not moved again; the first and last beats
% are left alone.
w = rr(:);  n = numel(w);
moved = false(n, 1);  skip = false(n, 1);
for it = 1:n
    c = classify_beats(w);
    j = find(c.ectopic & ~skip, 1);
    if isempty(j), break; end
    skip(j) = true;
    cand = [j j-1];  cost = inf(1, 2);
    span = max(1, j-2):min(n, j+2);
    for q = 1:2
        bq = cand(q);
        if bq < 2 || bq > n-1 || any(moved(bq-1:bq+1)), continue; end
        v = w(span);  i0 = bq - span(1) + 1;
        v(i0:i0+1) = (v(i0) + v(i0+1)) / 2;
        cost(q) = sum(abs(diff(v)));
    end
    [cmin, q] = min(cost);
    if isinf(cmin), continue; end
    b = cand(q);
    if abs(w(b) - w(b+1)) / (w(b) + w(b+1)) < min_shift, continue; end
    moved(b) = true;
    w(b:b+1) = (w(b) + w(b+1)) / 2;
end
end

function c = classify_beats(rr)
% Beat classification of Lipponen & Tarvainen (2019), following NeuroKit2's
% implementation (signal_fixpeaks, method 'Kubios'). rr: intervals (s); each
% output is true where the beat closing that interval is ectopic, long or
% short, extra (a false beat: its interval and the next sum to one beat), or
% where a beat is missing before it (interval ~ 2 beats); mrr: deviation
% from the 11-beat median (medrr) in units of Th2.
alpha = 5.2;  c1 = 0.13;  c2 = 0.17;  win = 91;  medwin = 11;
rr = rr(:);  n = numel(rr);
c = struct('ectopic', false(n,1), 'longshort', false(n,1), 'extra', false(n,1), 'missed', false(n,1), ...
    'medrr', movmedian(rr, medwin), 'mrr', zeros(n,1));
if n < 5, return; end

% dRR (difference with the previous interval), in units of Th1
drr = [0; diff(rr)];  drr(1) = mean(drr(2:end));
th1 = alpha * quartile_dev(abs(drr), win);
d   = drr ./ max(th1, eps);
% S12: the neighbor dRR of opposite direction (max if dRR > 0, min if < 0);
% S22: the next two dRR (min if dRR >= 0, max if < 0). Ends reflected.
dp   = [d(3); d(2); d; d(end-1); d(end-2)];
prv  = dp(2:n+1);  nxt = dp(4:n+3);  nxt2 = dp(5:n+4);
s12  = zeros(n,1);
up   = d > 0;  dn = d < 0;
s12(up) = max(prv(up), nxt(up));
s12(dn) = min(prv(dn), nxt(dn));
s22     = min(nxt, nxt2);
s22(dn) = max(nxt(dn), nxt2(dn));
% mRR: deviation from the 11-beat median, shorter intervals counted double,
% in units of Th2
medrr = c.medrr;
mrr   = rr - medrr;
mrr(mrr < 0) = 2 * mrr(mrr < 0);
th2   = alpha * quartile_dev(abs(mrr), win);
m     = mrr ./ max(th2, eps);
c.mrr = m;

% Decision flow (Lipponen & Tarvainen, 2019, Fig. 1)
i = 1;
while i <= n-2
    if abs(d(i)) <= 1, i = i + 1; continue; end
    if (d(i) > 1 && s12(i) < -c1*d(i) - c2) || (d(i) < -1 && s12(i) > -c1*d(i) + c2)
        c.ectopic(i) = true;
        i = i + 1;
        continue
    end
    cand = i;
    if abs(d(i+1)) < abs(d(i+2)), cand = [i i+1]; end
    for j = cand
        long  = d(j) > 1  && s22(j) < -1;
        short = d(j) < -1 && s22(j) > 1;
        if ~(long || short || abs(m(j)) > 3), continue; end
        if short && abs(rr(j) + rr(j+1) - medrr(j)) < th2(j)
            c.extra(j) = true;
        elseif long && abs(rr(j)/2 - medrr(j)) < th2(j)
            c.missed(j) = true;
        else
            c.longshort(j) = true;
        end
    end
    i = i + numel(cand);
end
end

function qd = quartile_dev(x, win)
% (Q3 - Q1) / 2 over a centered window of win values (shorter at the ends);
% quartiles interpolated linearly between order statistics
n  = numel(x);  h = (win - 1) / 2;
qd = zeros(n, 1);
for k = 1:n
    s = sort(x(max(1, k-h):min(n, k+h)));
    qd(k) = (lin_quantile(s, 0.75) - lin_quantile(s, 0.25)) / 2;
end
end

function q = lin_quantile(s, p)
pos = p * (numel(s) - 1) + 1;
lo  = floor(pos);
hi  = min(lo + 1, numel(s));
q   = s(lo) + (pos - lo) * (s(hi) - s(lo));
end

function local_med = local_median(nn, nn_t, win_size, med_rr)
% Median of the neighboring intervals within +/- win_size/2 s (excluding the
% interval itself), or med_rr without neighbors
n = length(nn);
local_med = zeros(n, 1);
for i = 1:n
    in_win     = abs(nn_t - nn_t(i)) <= win_size/2;
    in_win(i)  = false;
    if any(in_win)
        local_med(i) = median(nn(in_win), 'omitnan');
    else
        local_med(i) = med_rr;
    end
end
end