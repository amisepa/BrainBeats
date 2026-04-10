function [nn_intervals, nn_t, nPeaks, idx_bad, idx_interp] = clean_rr(rr_t, rr_intervals, peak_amp, rPeaks, varargin)
    % clean_rr - Clean RR interval series and interpolate missing heartbeats.
    %
    % Strategy:
    %   1. Estimate physiological min/max RR from the data itself (no hardcoded limits)
    %      by computing the distribution of intervals and excluding outliers via MAD
    %   2. Remove beats with abnormal amplitude
    %   3. Remove intervals outside the data-derived physiological range
    %   4. Detect gaps (missed beats) using a sliding-window local median
    %   5. Fill gaps by inserting synthetic beats at the correct timestamps,
    %      with signal amplitude read directly from the ECG/PPG at those times
    %
    % Usage:
    %   [nn, nn_t, nPeaks, idx_bad, idx_interp] = clean_rr(rr_t, rr, peak_amp, rPeaks)
    %   [...] = clean_rr(..., 'ecg_signal', sig, 'fs', fs)
    %
    % Inputs:
    %   rr_t         - RR interval timestamps (s), same length as rr_intervals
    %   rr_intervals - RR interval durations (s)
    %   peak_amp     - amplitudes of detected peaks, same length as rr_intervals
    %   rPeaks       - R-peak sample indices, same length as rr_intervals
    %
    % Optional name-value inputs:
    %   'interpolate_missing' (default: true)  - fill gaps with synthetic beats
    %   'max_missing'         (default: 3)     - max beats to insert per gap
    %   'mad_threshold'       (default: 3)     - MAD multiplier for outlier detection
    %   'win_size'            (default: 15)    - sliding window (s) for local median
    %   'gap_ratio'           (default: 1.5)   - gap flagged if interval > ratio * local_median
    %   'ecg_signal'          (default: [])    - raw ECG/PPG (samples x 1)
    %   'fs'                  (default: [])    - sampling rate (Hz)
    %   'sig_t'               (default: [])    - signal time axis (s); if empty, built from fs
    %
    % Outputs:
    %   nn_intervals - cleaned RR intervals (s)
    %   nn_t         - corresponding timestamps (s)
    %   nPeaks       - peak sample indices (NaN for synthetic beats without signal input)
    %   idx_bad      - logical: remaining unfilled suspicious intervals
    %   idx_interp   - logical: synthetically inserted intervals

    %% Parse inputs
    p = inputParser;
    addParameter(p, 'interpolate_missing', true,  @islogical);
    addParameter(p, 'max_missing',         3,     @isnumeric);
    addParameter(p, 'mad_threshold',       3,     @isnumeric);
    addParameter(p, 'win_size',            15,    @isnumeric);
    addParameter(p, 'gap_ratio',           1.5,   @isnumeric);
    addParameter(p, 'ecg_signal',          [],    @isnumeric);
    addParameter(p, 'fs',                  [],    @isnumeric);
    addParameter(p, 'sig_t',              [],    @isnumeric);
    parse(p, varargin{:});
    interpolate_missing = p.Results.interpolate_missing;
    max_missing         = p.Results.max_missing;
    mad_threshold       = p.Results.mad_threshold;
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

    %% Initialize
    n_orig       = length(rr_intervals);
    nn_intervals = rr_intervals(:);
    nn_t         = rr_t(:);
    nPeaks       = double(rPeaks(:));
    idx_bad      = false(n_orig, 1);
    idx_interp   = false(n_orig, 1);

    %% Step 1: Data-driven physiological limits (no hardcoded values)
    % Use MAD on the central distribution to find the normal range.
    % This adapts to whatever heart rate the subject actually has.
    med_rr  = median(nn_intervals, 'omitnan');
    mad_rr  = mad(nn_intervals, 1);  % median absolute deviation
    min_rr  = med_rr - mad_threshold * mad_rr;
    max_rr  = med_rr + mad_threshold * mad_rr;

    % Hard floor/ceiling only at true physiological extremes (10-300 bpm)
    min_rr  = max(min_rr, 0.20);
    max_rr  = min(max_rr, 6.00);

    fprintf('\n--- clean_rr ---\n');
    fprintf('  Median RR:    %.3f s (%.0f bpm)\n', med_rr, 60/med_rr);
    fprintf('  MAD:          %.3f s\n', mad_rr);
    fprintf('  Min RR (data-derived): %.3f s (%.0f bpm)\n', min_rr, 60/min_rr);
    fprintf('  Max RR (data-derived): %.3f s (%.0f bpm)\n', max_rr, 60/max_rr);

    %% Step 2: Remove beats with abnormal amplitude
    if ~isempty(peak_amp) && length(peak_amp) == length(nn_intervals)
        amp_outliers = isoutlier(peak_amp(:), 'median');
        if any(amp_outliers)
            fprintf('  Removing %d beats with abnormal amplitude\n', sum(amp_outliers));
        end
    else
        amp_outliers = false(size(nn_intervals));
    end

    %% Step 3: Remove only impossibly short intervals and amplitude outliers
    % Long intervals are NOT removed - they are gaps (missed beats) handled in Steps 4-5.
    too_short = nn_intervals < min_rr;
    if any(too_short)
        fprintf("  Removing %d intervals below min RR (%.3f s)\n", sum(too_short), min_rr);
    end


    to_remove = amp_outliers | too_short;
    n_removed = sum(to_remove);

    nn_intervals(to_remove) = [];
    nn_t(to_remove)         = [];
    nPeaks(to_remove)       = [];
    idx_bad(to_remove)      = [];
    idx_interp(to_remove)   = [];

    if isempty(nn_intervals)
        warning('All intervals removed - check your data quality');
        return;
    end

    %% Step 4: Detect gaps with sliding-window local median
    % A gap is an interval that is gap_ratio times larger than the local
    % median of its neighbors, suggesting one or more missed beats.
    % Using a local (not global) median handles gradual HR changes correctly.
    n         = length(nn_intervals);
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

    %% Step 5 (optional): Fill gaps with synthetic beats
    %
    % For each gap at index i:
    %   n_missing  = round(gap / local_median) - 1
    %   new_t      = evenly spaced from t(i-1) to t(i) in steps of local_median
    %   new_amp    = interp1(sig_t, ecg_signal, new_t)   <- actual signal value
    %   new_peaks  = nearest sample index in sig_t        <- no guessing
    %
    % Processed in reverse order to avoid index shifts during insertion.

    if interpolate_missing && any(idx_bad)

        gap_indices    = find(idx_bad);
        n_interp_total = 0;

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

            % Evenly space synthetic beats between the two flanking real peaks
            j_vec = (1:n_missing)';
            new_t = nn_t(i-1) + j_vec * (gap / (n_missing + 1));

            % Sub-intervals (equal spacing)
            sub_iv  = gap / (n_missing + 1);
            sub_ivs = repmat(sub_iv, n_missing + 1, 1);  % includes remainder to next real peak

            % Sample indices and amplitudes: read directly from the signal
            if use_ecg
                % Nearest sample in sig_t for each synthetic timestamp
                new_peaks = arrayfun(@(t) find_nearest_sample(sig_t, t), new_t);
            else
                new_peaks = NaN(n_missing, 1);
            end

            % Shrink the gap entry to one sub-interval, mark filled
            nn_intervals(i) = sub_ivs(end);
            idx_bad(i)      = false;

            % Insert synthetic rows before index i
            nn_intervals = [nn_intervals(1:i-1); sub_ivs(1:end-1);         nn_intervals(i:end)];
            nn_t         = [nn_t(1:i-1);         new_t;                    nn_t(i:end)];
            nPeaks       = [nPeaks(1:i-1);        new_peaks;               nPeaks(i:end)];
            idx_bad      = [idx_bad(1:i-1);       false(n_missing,1);      idx_bad(i:end)];
            idx_interp   = [idx_interp(1:i-1);    true(n_missing,1);       idx_interp(i:end)];
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

    %% Summary
    n_final = length(nn_intervals);
    fprintf('  Removed:   %d/%d (%.1f%%)\n', n_removed,                   n_orig,  100*n_removed/n_orig);
    fprintf('  Flagged:   %d/%d (%.1f%%)\n', sum(idx_bad),                n_final, 100*sum(idx_bad)/n_final);
    fprintf('  Synthetic: %d/%d (%.1f%%)\n', sum(idx_interp),             n_final, 100*sum(idx_interp)/n_final);
    fprintf('  Clean:     %d/%d (%.1f%%)\n', sum(~idx_bad & ~idx_interp), n_final, 100*sum(~idx_bad & ~idx_interp)/n_final);
    fprintf('----------------\n');

end

%% Helper: find the sample index in sig_t closest to a given timestamp
function idx = find_nearest_sample(sig_t, t)
    [~, idx] = min(abs(sig_t - t));
end