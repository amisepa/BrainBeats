% Detect R peaks from raw ECG with automatic boundary detection and offset correction
%
% Handles concatenated multi-session recordings by:
% 1. Detecting boundaries (large baseline/offset changes)
% 2. De-meaning each segment to remove DC offset
% 3. Processing each segment independently with P&T algorithm
% 4. Combining results
%
% Based on get_RR_v2.m with EEGLAB-style segment handling
% February 8, 2026 - Cedric Cannard

function [RR, RR_t, peaks, sig, tm, sign, HR] = get_RR_adaptive(signal, tm, params)

% Parameters
fs = params.fs;
sig_type = params.heart_signal;

if size(signal,1) < size(signal,2)
    signal = signal';
end

sign = [];
nSamp = size(signal,1);
tm = tm(:) / 1000;  % convert to seconds

if numel(tm) ~= nSamp
    error('tm must have same number of samples as signal')
end

%% ECG with Boundary Detection
if strcmpi(sig_type, 'ecg')

    % Parameters
    if isfield(params,'ecg_peakthresh')
        peakThresh = params.ecg_peakthresh;
    else
        peakThresh = .6;
    end
    if isfield(params,'ecg_searchback')
        search_back = params.ecg_searchback;
    else
        search_back = true;
    end
    if isfield(params,'ecg_refperiod')
        ref_period = params.ecg_refperiod;
    else
        ref_period = 0.25;
    end

    sig = signal(:);

    % Flatline check
    if prctile(abs(sig), 95) < 0.05
        error('ECG time series amplitude too small (likely flat line)')
    end

    fprintf('\n=== BOUNDARY DETECTION & OFFSET CORRECTION ===\n');
    
    % Detect boundaries (large baseline/offset shifts)
    if isfield(params, 'boundary_indices') && ~isempty(params.boundary_indices)
        boundary_indices = params.boundary_indices(:);
        fprintf(' - Using manual boundaries at samples: %s\n', mat2str(boundary_indices'));
    else
        boundary_indices = detect_baseline_shifts(sig, fs);
        if isempty(boundary_indices)
            fprintf(' - No boundaries detected\n');
        else
            fprintf(' - Detected %d boundaries at: %s\n', ...
                length(boundary_indices), mat2str(boundary_indices'));
        end
    end
    
    % Create segment definitions
    if isempty(boundary_indices)
        segments = struct('start', 1, 'end', nSamp);
    else
        boundaries = [1; boundary_indices(:); nSamp+1];
        segments = struct('start', {}, 'end', {});
        for i = 1:length(boundaries)-1
            segments(i).start = boundaries(i);
            segments(i).end = boundaries(i+1) - 1;
        end
    end
    
    fprintf('\n=== PROCESSING %d SEGMENTS ===\n', length(segments));
    
    % Process each segment
    all_peaks = [];
    all_polarity = [];
    sig_corrected = zeros(size(sig));
    
    for seg_idx = 1:length(segments)
        seg_start = segments(seg_idx).start;
        seg_end = segments(seg_idx).end;
        seg_sig = sig(seg_start:seg_end);
        seg_tm = tm(seg_start:seg_end);
        
        fprintf('\n--- Segment %d: samples %d-%d (%.1f s) ---\n', ...
            seg_idx, seg_start, seg_end, (seg_end-seg_start+1)/fs);
        
        % Skip if too short
        if length(seg_sig) < 2*fs
            fprintf('  Skipping (too short)\n');
            sig_corrected(seg_start:seg_end) = seg_sig;
            continue
        end
        
        % Remove DC offset (de-mean)
        seg_median = median(seg_sig);
        seg_sig_corrected = seg_sig - seg_median;
        sig_corrected(seg_start:seg_end) = seg_sig_corrected;
        
        fprintf('  DC offset removed: %.2f µV\n', seg_median);
        
        % Process with P&T algorithm
        [seg_peaks, seg_pol] = process_ecg_segment_pt(seg_sig_corrected, seg_tm, fs, ...
            peakThresh, search_back, ref_period);
        
        % Adjust peaks to global indices
        if ~isempty(seg_peaks)
            global_peaks = seg_peaks + seg_start - 1;
            all_peaks = [all_peaks; global_peaks(:)];
            all_polarity = [all_polarity; seg_pol];
        end
    end
    
    % Sort peaks chronologically
    [peaks, sort_idx] = sort(all_peaks);
    
    % Overall polarity (mode)
    if ~isempty(all_polarity)
        pol = mode(all_polarity);
        sign = pol * median(abs(sig_corrected(peaks)));
    else
        sign = [];
        pol = 1;
    end
    
    % Use corrected signal
    sig = sig_corrected;
    
    fprintf('\n=== COMBINED RESULTS ===\n');
    fprintf(' - Total R-peaks: %d\n', length(peaks));
    if ~isempty(all_polarity) && length(unique(all_polarity)) > 1
        fprintf(' - NOTE: Polarity varied across segments\n');
    end
    if pol < 0
        fprintf(' - Overall polarity: negative\n');
    else
        fprintf(' - Overall polarity: positive\n');
    end

    % RR intervals and HR
    RR = diff(peaks) ./ fs;
    RR_t = tm(peaks);
    HR = 60 ./ diff(RR_t);

    if ~isempty(RR)
        fprintf(' - Median RR: %.3f s (%.1f bpm)\n', median(RR), 60/median(RR));
    end

else
    error('Only ECG supported. Use get_RR.m for PPG.')
end

end

%% Detect baseline shifts (boundaries)
function boundary_indices = detect_baseline_shifts(sig, fs)
    % Detect MAJOR baseline/offset changes (like session concatenations)
    % Only triggers on very large shifts (hundreds of µV)
    
    window_size = round(5 * fs);  % 5-second windows (was 2s)
    step_size = round(2 * fs);    % 2-second steps (was 0.5s)
    
    n_windows = floor((length(sig) - window_size) / step_size) + 1;
    medians = zeros(n_windows, 1);
    mads = zeros(n_windows, 1);  % Median absolute deviation
    times = zeros(n_windows, 1);
    
    % Compute rolling median and MAD (baseline and variability)
    for i = 1:n_windows
        idx_start = (i-1) * step_size + 1;
        idx_end = min(idx_start + window_size - 1, length(sig));
        medians(i) = median(sig(idx_start:idx_end));
        mads(i) = mad(sig(idx_start:idx_end), 1);  % Median absolute deviation
        times(i) = idx_start + floor(window_size/2);
    end
    
    % Find large changes in median (baseline shifts)
    median_changes = abs(diff(medians));
    
    % Much more conservative threshold
    % Typical ECG amplitude is ~100-2000 µV, we want to catch shifts > 200 µV
    baseline_noise = median(median_changes);
    
    % Use both statistical and absolute thresholds
    stat_threshold = 15 * std(median_changes);  % 15 sigma (was 3)
    absolute_threshold = max(200, 2*median(mads));  % At least 200 µV or 2x typical variability
    threshold = max(stat_threshold, absolute_threshold);
    
    fprintf('  Boundary detection threshold: %.1f µV\n', threshold);
    fprintf('  Median baseline change: %.1f µV\n', baseline_noise);
    
    % Find peaks in median changes
    [~, shift_locs] = findpeaks(median_changes, ...
        'MinPeakHeight', threshold, ...
        'MinPeakDistance', round(30*fs/step_size));  % At least 30s apart (was 10s)
    
    if isempty(shift_locs)
        boundary_indices = [];
    else
        boundary_indices = round(times(shift_locs + 1));
        
        % Validate each boundary: check amplitude difference is substantial
        keep = true(size(boundary_indices));
        for i = 1:length(boundary_indices)
            b = boundary_indices(i);
            
            % Compare 10s before and after boundary
            before_start = max(1, b - 10*fs);
            before_end = b - 1;
            after_start = b;
            after_end = min(length(sig), b + 10*fs);
            
            if (before_end - before_start) < fs || (after_end - after_start) < fs
                keep(i) = false;
                continue;
            end
            
            med_before = median(sig(before_start:before_end));
            med_after = median(sig(after_start:after_end));
            shift_size = abs(med_before - med_after);
            
            % Require at least 150 µV shift (very conservative)
            if shift_size < 150
                keep(i) = false;
                fprintf('  Rejected boundary at %d: shift only %.1f µV\n', b, shift_size);
            else
                fprintf('  Validated boundary at %d: shift %.1f µV\n', b, shift_size);
            end
        end
        
        boundary_indices = boundary_indices(keep);
        
        % Filter boundaries too close to start/end
        boundary_indices(boundary_indices < 10*fs) = [];
        boundary_indices(boundary_indices > length(sig) - 10*fs) = [];
    end
end

%% Process single ECG segment with P&T algorithm
function [peaks, pol] = process_ecg_segment_pt(sig, tm, fs, peakThresh, search_back, ref_period)
    % Full P&T algorithm on a single segment
    
    nSamp = length(sig);
    
    % P&T constants
    med_smooth_nb_coef = round(fs/100);
    int_nb_coef = round(7*fs/256);
    if mod(med_smooth_nb_coef, 2) == 0
        med_smooth_nb_coef = med_smooth_nb_coef + 1;
    end
    
    % P&T operations
    dffecg = [0; diff(sig)];
    sqrecg = dffecg.^2;
    b_int = ones(int_nb_coef,1) / int_nb_coef;
    intecg = filtfilt(b_int, 1, sqrecg);
    mdfint = medfilt1(intecg, med_smooth_nb_coef);
    
    % Match length
    if numel(mdfint) < numel(sig)
        mdfint = [mdfint; zeros(numel(sig) - numel(mdfint), 1)];
    elseif numel(mdfint) > numel(sig)
        mdfint = mdfint(1:numel(sig));
    end
    
    % P&T threshold
    if nSamp/fs > 90
        xs = sort(mdfint(fs:min(fs*90, end)));
    else
        xs = sort(mdfint(max(1,fs):end));
    end
    
    if nSamp/fs > 10
        ind_xs = ceil(98/100*length(xs));
    else
        ind_xs = ceil(99/100*length(xs));
    end
    en_thres = xs(ind_xs);
    
    % Candidate regions
    poss_reg = mdfint > (peakThresh * en_thres);
    
    if ~any(poss_reg)
        peaks = [];
        pol = 1;
        return
    end
    
    % Search-back for missed beats
    if search_back
        indAboveThreshold = find(poss_reg);
        if numel(indAboveThreshold) > 2
            RRv = diff(tm(indAboveThreshold));
            RRv = RRv(RRv > 0.01);
            if ~isempty(RRv)
                medRRv = median(RRv);
                indMissedBeat = find(diff(tm(indAboveThreshold)) > 1.5*medRRv);
                indStart = indAboveThreshold(indMissedBeat);
                indEnd = indAboveThreshold(indMissedBeat+1);
                for i = 1:numel(indStart)
                    poss_reg(indStart(i):indEnd(i)) = ...
                        mdfint(indStart(i):indEnd(i)) > (0.5 * peakThresh * en_thres);
                end
            end
        end
    end
    
    % Segment boundaries
    left = find(diff([0; poss_reg])==1);
    right = find(diff([poss_reg; 0])==-1);
    nb_peaks = numel(left);
    
    if nb_peaks == 0
        peaks = [];
        pol = 1;
        return
    end
    
    % Find peaks and determine polarity
    peaks = zeros(nb_peaks, 1);
    for i = 1:nb_peaks
        a = left(i);
        b = right(i);
        seg = sig(a:b);
        
        [vmax, imax] = max(seg);
        [vmin, imin] = min(seg);
        
        % Simple polarity: which is larger in absolute value?
        if abs(vmax) > abs(vmin)
            peaks(i) = a + imax - 1;
        else
            peaks(i) = a + imin - 1;
        end
    end
    
    % Determine overall polarity
    peak_values = sig(peaks);
    if median(peak_values) > 0
        pol = 1;
    else
        pol = -1;
    end
    
    % Keep only peaks of dominant polarity
    if pol > 0
        keep = sig(peaks) > 0;
    else
        keep = sig(peaks) < 0;
    end
    peaks = peaks(keep);
    
    % Refine peaks
    for i = 1:length(peaks)
        a = left(find(left <= peaks(i), 1, 'last'));
        b = right(find(right >= peaks(i), 1, 'first'));
        seg = sig(a:b);
        if pol > 0
            [~, idx] = max(seg);
        else
            [~, idx] = min(seg);
        end
        peaks(i) = a + idx - 1;
    end
    
    % Enforce refractory period
    if length(peaks) > 1
        refSamples = round(ref_period * fs);
        keep = true(size(peaks));
        for k = 2:length(peaks)
            if (peaks(k) - peaks(k-1)) < refSamples
                if abs(sig(peaks(k))) > abs(sig(peaks(k-1)))
                    keep(k-1) = false;
                else
                    keep(k) = false;
                end
            end
        end
        peaks = peaks(keep);
    end
    
    % Micro-refinement
    micro = max(1, round(0.015 * fs));
    for k = 1:length(peaks)
        p = peaks(k);
        w1 = max(1, p - micro);
        w2 = min(nSamp, p + micro);
        if pol > 0
            [~, ii] = max(sig(w1:w2));
        else
            [~, ii] = min(sig(w1:w2));
        end
        peaks(k) = w1 + ii - 1;
    end
    peaks = unique(peaks, 'stable');
    
    fprintf('  Detected %d peaks (polarity: %s)\n', ...
        length(peaks), iif(pol > 0, 'pos', 'neg'));
end

% Inline if
function out = iif(condition, true_val, false_val)
    if condition
        out = true_val;
    else
        out = false_val;
    end
end