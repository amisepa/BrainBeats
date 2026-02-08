% Detect R peaks from raw ECG signals using ADAPTIVE THRESHOLDS
% 
% Adapted from PPG detection methods: continuously updates peak height 
% and distance thresholds based on recent signal characteristics
%
% February 8, 2026 - Cedric Cannard
%   - Implemented adaptive threshold approach for multi-session recordings
%   - Updates thresholds every 5s based on past 20s of detected peaks
%   - Handles amplitude and baseline changes automatically without segmentation
%
% Copyright (C), BrainBeats, Cedric Cannard, 2023-2026

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

%% ECG with Adaptive Thresholds
if strcmpi(sig_type, 'ecg')

    % Parameters
    if isfield(params,'ecg_refperiod')
        ref_period = params.ecg_refperiod;
    else
        ref_period = 0.25; % refractory period
    end

    % Use column vector
    sig = signal(:);

    % flatline check 
    if prctile(abs(sig), 95) < 0.05
        error('ECG time series amplitude too small (likely flat line)')
    end

    fprintf('\n=== ADAPTIVE PEAK DETECTION ===\n');
    
    % High-pass filter entire signal to remove baseline drift
    if fs >= 100
        [b_hp, a_hp] = butter(2, 0.5/(fs/2), 'high');
        sig = filtfilt(b_hp, a_hp, sig);
        fprintf(' - Applied 0.5 Hz high-pass filter\n');
    end
    
    % Calibration period: use 10-30 seconds (skip first 5s for stability)
    calib_start = min(round(5*fs), round(nSamp*0.1));  % Skip first 5s or 10%
    calib_duration = min(25, floor(nSamp/fs) - calib_start/fs - 1);
    calib_end = min(calib_start + round(calib_duration * fs), nSamp);
    calib_data = sig(calib_start:calib_end);
    
    fprintf(' - Calibration: %.1f-%.1f s (%d samples)\n', ...
        calib_start/fs, calib_end/fs, length(calib_data));
    
    % Estimate initial thresholds
    [init_height, pol] = estimate_peak_height(calib_data);
    
    % Conservative initial distance: assume HR between 40-150 bpm
    % Min distance = 60/150 bpm = 0.4s, Max = 60/40 bpm = 1.5s
    init_distance = round(0.4 * fs);  % ~150 bpm max (was 0.33s/180bpm)
    
    fprintf(' - Initial height threshold: %.2f\n', init_height);
    fprintf(' - Initial distance: %d samples (%.2f s, max ~%.0f bpm)\n', ...
        init_distance, init_distance/fs, 60/(init_distance/fs));
    fprintf(' - Polarity: %s\n', iif(pol > 0, 'positive', 'negative'));
    
    % Adaptive parameters
    update_interval = 3 * fs;  % Update every 5 seconds
    lookback_window = 10 * fs;  % Use past 20 seconds
    
    % P&T constants for QRS detection
    med_smooth_nb_coef = round(fs/100);
    int_nb_coef = round(7*fs/256);
    if mod(med_smooth_nb_coef, 2) == 0
        med_smooth_nb_coef = med_smooth_nb_coef + 1;
    end
    
    % Initialize
    peaks = [];
    current_height = init_height;
    current_dist = init_distance;
    
    % Process entire signal with P&T to get QRS energy
    fprintf(' - Computing P&T QRS energy envelope...\n');
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
    
    % Estimate initial energy threshold
    if nSamp/fs > 90
        xs = sort(mdfint(fs:fs*90));
    else
        xs = sort(mdfint(fs:end));
    end
    if nSamp/fs > 10
        ind_xs = ceil(98/100*length(xs));
    else
        ind_xs = ceil(99/100*length(xs));
    end
    en_thres = xs(ind_xs);
    current_energy_thresh = 0.3 * en_thres;  % Start conservative
    
    fprintf(' - Initial energy threshold: %.2f\n', current_energy_thresh);
    
    % Process in chunks
    chunk_size = update_interval;
    n_chunks = ceil(nSamp / chunk_size);
    
    fprintf(' - Processing %d chunks with adaptive P&T thresholds...\n', n_chunks);
    
    for i_chunk = 1:n_chunks
        chunk_start = (i_chunk - 1) * chunk_size + 1;
        chunk_end = min(i_chunk * chunk_size, nSamp);
        
        % Find candidate QRS regions based on energy
        chunk_energy = mdfint(chunk_start:chunk_end);
        poss_reg = chunk_energy > current_energy_thresh;
        
        if ~any(poss_reg)
            continue
        end
        
        % Get segment boundaries
        left = find(diff([0; poss_reg])==1);
        right = find(diff([poss_reg; 0])==-1);
        
        % Find peak in each QRS segment
        chunk_data = sig(chunk_start:chunk_end);
        chunk_peaks = [];
        
        for i_seg = 1:length(left)
            a = left(i_seg);
            b = right(i_seg);
            seg = chunk_data(a:b);
            
            if pol > 0
                [~, pk_idx] = max(seg);
            else
                [~, pk_idx] = min(seg);
            end
            
            chunk_peaks = [chunk_peaks; a + pk_idx - 1];
        end
        
        % Enforce minimum distance within chunk
        if length(chunk_peaks) > 1
            keep = true(size(chunk_peaks));
            for k = 2:length(chunk_peaks)
                if chunk_peaks(k) - chunk_peaks(k-1) < current_dist
                    % Keep peak with higher energy
                    if mdfint(chunk_start + chunk_peaks(k) - 1) > mdfint(chunk_start + chunk_peaks(k-1) - 1)
                        keep(k-1) = false;
                    else
                        keep(k) = false;
                    end
                end
            end
            chunk_peaks = chunk_peaks(keep);
        end
        
        % Adjust to global indices
        chunk_peaks = chunk_peaks + chunk_start - 1;
        peaks = [peaks; chunk_peaks(:)];
        
        % Update thresholds based on past 20s
        if chunk_end > lookback_window
            lookback_start = chunk_end - lookback_window;
            recent_peaks = peaks(peaks >= lookback_start & peaks <= chunk_end);
            
            if length(recent_peaks) >= 5
                % Update height threshold
                if pol > 0
                    recent_heights = sig(recent_peaks);
                else
                    recent_heights = -sig(recent_peaks);
                end
                current_height = 0.4 * median(recent_heights);  % 40% of median
                
                % Update energy threshold based on recent QRS energy
                recent_energy = mdfint(recent_peaks);
                current_energy_thresh = 0.4 * median(recent_energy);
                
                % Update distance threshold
                if length(recent_peaks) >= 2
                    recent_dist = diff(recent_peaks);
                    median_dist = median(recent_dist);
                    current_dist = round(0.6 * median_dist);
                    
                    % Enforce physiological limits (40-150 bpm)
                    min_dist = round(0.4*fs);   % Max 150 bpm
                    max_dist = round(1.5*fs);   % Max 40 bpm
                    current_dist = max(current_dist, min_dist);
                    current_dist = min(current_dist, max_dist);
                end
            end
        end
    end
    
    % Remove duplicates
    peaks = unique(peaks);
    
    % Enforce refractory period
    if length(peaks) > 1
        refSamples = round(ref_period * fs);
        keep = true(size(peaks));
        for k = 2:length(peaks)
            if (peaks(k) - peaks(k-1)) < refSamples
                if pol > 0
                    if sig(peaks(k)) > sig(peaks(k-1))
                        keep(k-1) = false;
                    else
                        keep(k) = false;
                    end
                else
                    if sig(peaks(k)) < sig(peaks(k-1))
                        keep(k-1) = false;
                    else
                        keep(k) = false;
                    end
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
    
    sign = pol * median(abs(sig(peaks)));
    
    fprintf('\n=== RESULTS ===\n');
    fprintf(' - Detected %d R-peaks\n', length(peaks));
    if pol < 0
        fprintf(" - Polarity: negative\n");
    else
        fprintf(" - Polarity: positive\n");
    end

    % RR intervals and HR
    RR = diff(peaks) ./ fs;
    RR_t = tm(peaks);
    HR = 60 ./ diff(RR_t);

    if ~isempty(RR)
        fprintf(' - Median RR: %.3f s (%.1f bpm)\n', median(RR), 60/median(RR));
    end

else
    error('Only ECG supported in adaptive version. Use get_RR.m for PPG.')
end

end

%% Helper: Estimate peak height threshold
function [peak_height, polarity] = estimate_peak_height(data)
    % Robust estimation of peak height and polarity
    
    % Use findpeaks with minimum distance to avoid noise
    fs_approx = 250;  % Assume ~250 Hz for min distance
    min_pk_dist = round(0.4 * fs_approx);  % At least 0.4s apart
    
    % Check both positive and negative peaks
    [pos_peaks, ~] = findpeaks(data, 'MinPeakDistance', min_pk_dist);
    [neg_peaks, ~] = findpeaks(-data, 'MinPeakDistance', min_pk_dist);
    
    if isempty(pos_peaks)
        pos_height = 0;
    else
        pos_height = median(pos_peaks);
    end
    
    if isempty(neg_peaks)
        neg_height = 0;
    else
        neg_height = median(neg_peaks);
    end
    
    % Determine polarity
    if pos_height > neg_height
        polarity = 1;
        peak_height = 0.4 * pos_height;  % 40% of median peak (more conservative)
    else
        polarity = -1;
        peak_height = 0.4 * neg_height;
    end
    
    % Fallback if no peaks found
    if peak_height == 0
        peak_height = 0.25 * prctile(abs(data), 95);
        polarity = 1;
    end
end

% Simple inline if
function out = iif(condition, true_val, false_val)
    if condition
        out = true_val;
    else
        out = false_val;
    end
end