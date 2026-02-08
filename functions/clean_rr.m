function [nn_intervals, nn_t, idx_rem, idx_corr] = clean_rr(rr_t, rr_intervals, peak_amp)
    % clean_rr - Detects and corrects abnormal heartbeats in RR interval time series.
    %
    % Strategy:
    %   1. Remove beats with abnormal amplitude (likely noise, not true R-peaks)
    %   2. Remove physiologically impossible short intervals (ectopic/false positive)
    %   3. Interpolate long gaps (missed beats)
    %   4. Flag remaining ectopics using adaptive thresholding based on local RR
    %
    % INPUT:
    %   rr_t          - Vector of timestamps (seconds) for each RR interval
    %   rr_intervals  - Vector of RR intervals in seconds
    %   peak_amp      - Vector of signal amplitude at detected R-peaks
    %
    % OUTPUT:
    %   nn_intervals  - Corrected NN intervals (seconds)
    %   nn_t          - Updated timestamps
    %   idx_rem       - Indices of removed intervals (from original input)
    %   idx_corr      - Indices of interpolated intervals (in output)

    %% Parameters (physiologically motivated)
    min_rr = 0.33;          % 180 bpm - absolute physiological minimum
    max_rr = 2.0;           % 30 bpm - reasonable maximum for resting
    amp_threshold = 5;      % MAD multipliers for amplitude outliers
    ectopic_threshold = 0.25; % 25% deviation from local median = likely ectopic
    local_window = 11;      % Window for local median (odd number)
    
    %% Initialize
    n_orig = length(rr_intervals);
    nn_intervals = rr_intervals(:);
    nn_t = rr_t(:);
    peak_amp = peak_amp(:);
    
    idx_rem = false(n_orig, 1);
    
    %% Step 1: Remove beats with abnormal amplitude (false R-peak detections)
    % Use MAD-based outlier detection - robust to non-normal distributions
    med_amp = median(peak_amp, 'omitnan');
    mad_amp = mad(peak_amp, 1);  % MAD with scaling factor
    
    if mad_amp > 0
        amp_outliers = abs(peak_amp - med_amp) > amp_threshold * mad_amp * 1.4826;
    else
        amp_outliers = false(size(peak_amp));
    end
    
    if any(amp_outliers)
        fprintf('Removing %d beats with abnormal amplitude\n', sum(amp_outliers));
    end
    
    %% Step 2: Remove physiologically impossible short intervals
    % These are almost certainly false positives or PVCs
    short_rr = nn_intervals < min_rr;
    
    if any(short_rr)
        fprintf('Removing %d intervals < %.0f ms (>%.0f bpm)\n', ...
            sum(short_rr), min_rr*1000, 60/min_rr);
    end
    
    %% Combine removal criteria
    to_remove = amp_outliers | short_rr;
    idx_rem(to_remove) = true;
    
    nn_intervals(to_remove) = [];
    nn_t(to_remove) = [];
    
    if isempty(nn_intervals)
        warning('All intervals removed - check your data quality');
        idx_corr = [];
        return;
    end
    
    %% Step 3: Detect ectopic beats using adaptive local threshold
    % An ectopic beat creates a short-long or long-short pattern
    % We look for intervals that deviate significantly from local median
    
    idx_corr = false(length(nn_intervals), 1);
    
    % Compute local median RR (robust estimate of local heart rate)
    local_med = movmedian(nn_intervals, local_window, 'omitnan', 'Endpoints', 'shrink');
    
    % Relative deviation from local median
    rel_deviation = abs(nn_intervals - local_med) ./ local_med;
    
    % Flag as ectopic if deviation exceeds threshold
    ectopics = rel_deviation > ectopic_threshold;
    
    % Additional check: ectopics typically come in pairs (short-long pattern)
    % A true ectopic followed by compensatory pause: RR1 + RR2 ≈ 2 * normal RR
    % This helps distinguish from normal sinus arrhythmia
    for i = 2:length(ectopics)-1
        if ectopics(i)
            % Check if this is part of a compensatory pattern
            sum_rr = nn_intervals(i-1) + nn_intervals(i);
            expected_sum = 2 * local_med(i);
            
            % If the sum is close to expected (compensatory), likely true ectopic
            % If not, might just be normal variability - be conservative
            if abs(sum_rr - expected_sum) / expected_sum > 0.15
                % Doesn't fit ectopic pattern - might be normal variability
                % Only keep flagged if deviation is very large
                if rel_deviation(i) < ectopic_threshold * 1.5
                    ectopics(i) = false;
                end
            end
        end
    end
    
    if any(ectopics)
        fprintf('Interpolating %d suspected ectopic beats\n', sum(ectopics));
        idx_corr(ectopics) = true;
        
        % Interpolate ectopics using PCHIP (shape-preserving, no overshoot)
        valid_idx = ~ectopics;
        if sum(valid_idx) >= 2
            nn_intervals = interp1(nn_t(valid_idx), nn_intervals(valid_idx), ...
                nn_t, 'pchip', 'extrap');
        end
    end
    
    %% Step 4: Handle long gaps (missed beats)
    % If RR > max_rr, beats were likely missed
    long_gaps = nn_intervals > max_rr;
    
    if any(long_gaps)
        fprintf('Detected %d long gaps (>%.1f s) - likely missed beats\n', ...
            sum(long_gaps), max_rr);
        
        % For long gaps, we can either:
        % (a) Interpolate (simple but loses timing info)
        % (b) Insert estimated beats (better for HRV)
        % Here we interpolate, but flag them
        
        idx_corr(long_gaps) = true;
        
        valid_idx = ~long_gaps;
        if sum(valid_idx) >= 2
            nn_intervals = interp1(nn_t(valid_idx), nn_intervals(valid_idx), ...
                nn_t, 'pchip', 'extrap');
        end
    end
    
    %% Final bounds check
    % Clamp any remaining extreme values (from extrapolation artifacts)
    nn_intervals = max(nn_intervals, min_rr);
    nn_intervals = min(nn_intervals, max_rr);
    
    %% Summary
    if isempty(idx_corr) || ~any(idx_corr)
        idx_corr = [];
        fprintf('No intervals required interpolation\n');
    else
        fprintf('Total corrected: %d/%d intervals (%.1f%%)\n', ...
            sum(idx_corr), length(nn_intervals), 100*sum(idx_corr)/length(nn_intervals));
    end
    
    fprintf('Total removed: %d/%d intervals (%.1f%%)\n', ...
        sum(idx_rem), n_orig, 100*sum(idx_rem)/n_orig);
end