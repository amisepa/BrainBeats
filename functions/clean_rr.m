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
%   4. Gaps: intervals > gap_ratio x the median of neighbors within win_size.
%   5. Each gap is split into equal intervals by inserting
%      round(gap / local median) - 1 synthetic beats.
%
% Notes:
%   amp_mad_threshold is deliberately lenient (15 MAD, not isoutlier's 3):
%   ECG amplitude varies with respiration, posture and electrode contact;
%   only disconnections or rail artifacts should be caught here. Long
%   intervals are never removed, only flagged or filled (steps 4-5).
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

% Step 4: Detect gaps against the local median of neighboring intervals
% (+/- win_size/2 s, excluding the interval itself)
n  = length(nn_intervals);
local_med = zeros(n, 1);
for i = 1:n
    in_win     = abs(nn_t - nn_t(i)) <= win_size/2;
    in_win(i)  = false;
    if any(in_win)
        local_med(i) = median(nn_intervals(in_win), 'omitnan');
    else
        local_med(i) = med_rr;
    end
end

idx_bad = nn_intervals > gap_ratio * local_med;

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

%% Helper
function idx = find_nearest_sample(sig_t, t)
[~, idx] = min(abs(sig_t - t));
end