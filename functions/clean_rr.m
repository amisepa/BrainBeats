function [nn_intervals, nn_t, nPeaks, idx_bad] = clean_rr(rr_t, rr_intervals, peak_amp, rPeaks)
    % clean_rr - Detects and removes abnormal heartbeats in RR interval time series.
    %
    % Strategy:
    %   1. Remove beats with abnormal amplitude (likely noise, not true R-peaks)
    %   2. Remove physiologically impossible short intervals (ectopic/false positive)
    %   3. FLAG (don't interpolate) remaining suspicious intervals
    %
    % NOTE: We mark bad intervals rather than interpolate, because interpolating
    % interval durations creates artificial peak timings that don't match the ECG.
    
    %% Parameters (more conservative)
    min_rr = 0.4;              % 150 bpm - physiological minimum
    max_rr = 1.5;              % 40 bpm - reasonable maximum
    % amp_threshold = 3;         % MAD multipliers 
    % ectopic_threshold = 0.35;  % 35% deviation (was too strict at 25%)
    % local_window = 61;         % Larger window for more stable local median
    
    %% Initialize
    n_orig = length(rr_intervals);
    nn_intervals = rr_intervals(:);
    nn_t = rr_t(:);
    peak_amp = peak_amp(:);
    nPeaks = rPeaks(:);
    
    idx_rem = false(n_orig, 1);
    idx_bad = false(n_orig, 1);  % Flagged as suspicious but not removed
    
    %% Step 1: Remove beats with EXTREMELY abnormal amplitude
    % Be very conservative - only remove clear false detections
    % med_amp = median(peak_amp, 'omitnan');
    % mad_amp = mad(peak_amp, 1);
    % 
    % if mad_amp > 0
    %     % Only remove if amplitude is VERY different (extreme outliers)
    %     amp_outliers = abs(peak_amp - med_amp) > amp_threshold * mad_amp * 1.4826;
    % else
    %     amp_outliers = false(size(peak_amp));
    % end
    peak_amp(1) = [];
    amp_outliers = isoutlier(peak_amp, 'median'); 
    
    if any(amp_outliers)
        fprintf('Removing %d beats with extreme amplitude outliers\n', sum(amp_outliers));
    end
    
    %% Step 2: Remove physiologically impossible short intervals
    short_rr = nn_intervals < min_rr;
    
    if any(short_rr)
        fprintf('Removing %d intervals < %.0f ms (>%.0f bpm)\n', ...
            sum(short_rr), min_rr*1000, 60/min_rr);
    end
    
    %% Combine removal criteria and remove

    to_remove = amp_outliers | short_rr;
    idx_rem(to_remove) = true;
    
    nn_intervals(to_remove) = [];
    nn_t(to_remove) = [];
    nPeaks(to_remove) = [];
    idx_bad(to_remove) = [];  % Update idx_bad to match
    
    if isempty(nn_intervals)
        warning('All intervals removed - check your data quality');
        return;
    end
    
    %% Step 3: FLAG (don't interpolate) suspicious ectopic-like intervals
    % Compute local median using larger window
    % local_med = movmedian(nn_intervals, local_window, 'omitnan', 'Endpoints', 'shrink');
    
    % % Relative deviation from local median
    % rel_deviation = abs(nn_intervals - local_med) ./ local_med;
    
    % % Flag as suspicious if deviation exceeds threshold
    % suspicious = rel_deviation > ectopic_threshold;
    % 
    % % Additional constraint: must also be part of short-long or long-short pattern
    % confirmed_ectopic = false(size(suspicious));
    % for i = 2:length(suspicious)-1
    %     if suspicious(i)
    %         % Check for compensatory pattern
    %         % Either: short followed by long, or long followed by short
    %         prev_ratio = nn_intervals(i) / nn_intervals(i-1);
    %         next_ratio = nn_intervals(i) / nn_intervals(i+1);
    % 
    %         % Classic ectopic pattern: big change followed by compensation
    %         if (prev_ratio < 0.7 && next_ratio > 1.3) || ...
    %            (prev_ratio > 1.3 && next_ratio < 0.7)
    %             confirmed_ectopic(i) = true;
    %         end
    %     end
    % end
    
    % Also flag very long gaps (missed beats)
    long_gaps = nn_intervals > max_rr;
    
    % Combine all "bad" intervals
    % idx_bad = confirmed_ectopic | long_gaps;
    idx_bad = long_gaps;
    
    if any(idx_bad)
        fprintf('Flagging %d suspicious intervals (%.1f%%)\n', ...
            sum(idx_bad), 100*sum(idx_bad)/length(idx_bad));
        % fprintf('  - %d ectopic-like patterns\n', sum(confirmed_ectopic));
        fprintf('  - %d long gaps (missed beats)\n', sum(long_gaps));
    else
        fprintf('No suspicious intervals detected\n');
    end
    
    %% Summary
    fprintf('Total removed: %d/%d intervals (%.1f%%)\n', ...
        sum(idx_rem), n_orig, 100*sum(idx_rem)/n_orig);
    fprintf('Total flagged as suspicious: %d/%d remaining intervals (%.1f%%)\n', ...
        sum(idx_bad), length(idx_bad), 100*sum(idx_bad)/length(idx_bad));
    
    % NOTE: idx_bad can be used downstream to exclude intervals from HRV analysis
    % or to guide manual inspection, but peaks remain at their detected locations
end