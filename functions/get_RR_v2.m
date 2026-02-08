% Detect R peaks from raw ECG signals and heartbeat onsets (pulse
% waveforms) from PPG signals.
%
% ECG:
%   ECG signal is bandpassed filtered using a custom filter that provides
%   great performance. QRS detector is based on the P&T method. Energy
%   threshold is estimated at 99% of amplitude distribution to avoid crash
%   due to large bumps. 1 s removed for choosing thresh because of filter
%   lag and often contains artifacts. A search-back algorithm detects missed
%   peaks by lowering the threshold in periods where the RR interval
%   variability (RRv) is > 1.5*medianRRv. The polarity of the R-peaks is
%   detected over the first 30 s using max of abs values. The rest of the
%   peaks are compared to this sign, preventing from alternating between
%   positive/negative detections.
%
% PPG:
%   Algorithms designed for adult human signals sampled at 125 Hz, but works
%   with any sampling rate (using on-the-fly resampling). Signals shorter
%   than 5 min long are rescaled.
%   Original code from qppg from the Physionet Cardiovascular Signal Processing
%   Toolbox. Original authors: W. Zong (1998), revised by G. Moody (2010), Qiao Li
%   (2010-2011), Adriana Vest (2011), and Giulia Da Poian (2018).
%
% INPUTS:
%   signal      - raw ECG or PPG signal
%   tm          - time vector (in milliseconds)
%   params      - params structure containing the following fields:
%                   fs (sample rate) and heart_signal ('ecg' or 'ppg')
%
% OUTPUTS:
%   RR          - RR intervals (in s)
%   RR_t        - time vector
%   peaks       - R peaks (samples)
%   sig         - ECG signal after processing
%   pol         - ECG signal polarity
%   HR          - Heart rate (in beats/min; bpm)
%
% Example:
%
% When using this code, please cite:
%   Vest et al. (2018). An Open Source Benchmarked Toolbox for Cardiovascular
%   Waveform and Interval Analysis. Physiological measurement.
%
% ECG detection improvements: January 12, 2026 by Cedric Cannard
%   - Enforced consistent signal orientation by converting the ECG signal
%     to a column vector once at the start of the ECG branch, eliminating
%     downstream dimension and indexing errors.
%   - Reworked the Pan–Tompkins processing to be fully zero-phase:
%       - Padded the differentiated signal to preserve length.
%       - Replaced the integration filter with filtfilt to remove phase delay.
%       - Ensured the median filter window length is odd.
%       - Removed all manual delay compensation and circshift logic.
%   - Corrected polarity detection and peak finding to operate on the
%     original ECG signal rather than intermediate or transposed variables,
%     ensuring correct peak selection regardless of signal inversion.
%   - Added a post-detection peak refinement step that recenters each coarse
%     peak by searching for the true local extremum within a short QRS-scale window.
%   - Fixed search-back logic inconsistencies:
%       - Used logical masks instead of mixed index spaces.
%       - Properly flipped backward-search results back into forward time
%         before combining.
%       - Clarified forward vs backward variable naming.
%   - Clamped all search windows to valid signal bounds and aligned
%     intermediate vectors to identical lengths to prevent off-by-one and
%     assignment errors.
%   - Fixed minor typos and incorrect variable references in diagnostic
%     and warning messages.
%
% February 7, 2026 - Cedric Cannard
%   - Added automatic boundary detection for concatenated session recordings
%   - Implemented segment-wise filtering following EEGLAB approach (Widmann et al. 2015)
%   - "Filters must not be applied across signal discontinuities"
%   - Each segment processed independently with adaptive thresholds
%   - Handles amplitude and baseline changes between sessions automatically
%
% Copyright (C), BrainBeats, Cedric Cannard, 2023-2026

function [RR, RR_t, peaks, sig, tm, sign, HR] = get_RR_v2(signal, tm, params)

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


%% ECG
if strcmpi(sig_type, 'ecg')

    % Parameters
    if isfield(params,'ecg_peakthresh')
        peakThresh = params.ecg_peakthresh;
    else
        peakThresh = .3;
    end
    if isfield(params,'ecg_searchback')
        search_back = params.ecg_searchback;
    else
        search_back = true;
    end
    if isfield(params,'ecg_refperiod')
        ref_period = params.ecg_refperiod;
    else
        ref_period = 0.25; % refractory period
    end

    % Use a consistent orientation: COLUMN vectors everywhere in ECG branch
    sig = signal(:);

    % flatline check 
    if prctile(abs(sig), 95) < 0.05
        error('ECG time series amplitude too small (likely flat line)')
    end

    % ============================================================
    % AUTOMATIC BOUNDARY DETECTION (like EEGLAB)
    % ============================================================
    fprintf('\n=== AUTOMATIC BOUNDARY DETECTION ===\n');
    
    % Check if user manually specified boundaries
    if isfield(params, 'boundary_indices') && ~isempty(params.boundary_indices)
        boundary_indices = params.boundary_indices(:);
        fprintf(' - Using user-specified boundaries at samples: %s\n', mat2str(boundary_indices'));
    elseif isfield(params, 'disable_boundary_detection') && params.disable_boundary_detection
        boundary_indices = [];
        fprintf(' - Boundary detection disabled by user\n');
    else
        % Automatic detection
        boundary_indices = detect_signal_boundaries(sig, fs);
    end
    
    if isempty(boundary_indices)
        fprintf(' - No boundaries detected: processing as single continuous segment\n');
        segments = {struct('start', 1, 'end', length(sig))};
    else
        fprintf(' - Detected %d boundary/boundaries at sample(s): %s\n', ...
            length(boundary_indices), mat2str(boundary_indices));
        % Create segment definitions
        boundaries = [1; boundary_indices(:); length(sig)+1];
        segments = cell(1, length(boundaries)-1);
        for i = 1:length(boundaries)-1
            segments{i} = struct('start', boundaries(i), 'end', boundaries(i+1)-1);
        end
    end
    
    % ============================================================
    % SEGMENT-WISE FILTERING (like EEGLAB eeg_filtnew)
    % "Filters must not be applied across signal discontinuities"
    % - Widmann et al. 2015
    % ============================================================
    fprintf('\n=== SEGMENT-WISE FILTERING ===\n');
    sig_filtered = zeros(size(sig));
    
    for seg_idx = 1:length(segments)
        seg_start = segments{seg_idx}.start;
        seg_end = segments{seg_idx}.end;
        seg_data = sig(seg_start:seg_end);
        
        fprintf(' - Segment %d: samples %d-%d (%.1f s)\n', ...
            seg_idx, seg_start, seg_end, (seg_end-seg_start+1)/fs);
        
        % High-pass filter this segment only (0.5 Hz to remove baseline drift)
        if fs >= 100 && length(seg_data) > 3*fs  % Need sufficient data for filter
            [b_hp, a_hp] = butter(2, 0.5/(fs/2), 'high');
            seg_filtered = filtfilt(b_hp, a_hp, seg_data);
            fprintf('   Applied 0.5 Hz high-pass filter\n');
        else
            seg_filtered = seg_data;
            fprintf('   No filtering (segment too short or low fs)\n');
        end
        
        sig_filtered(seg_start:seg_end) = seg_filtered;
    end
    
    % Replace original signal with filtered version
    sig = sig_filtered;
    
    % ============================================================
    % SEGMENT-WISE R-PEAK DETECTION
    % ============================================================
    fprintf('\n=== SEGMENT-WISE R-PEAK DETECTION ===\n');
    
    all_peaks = zeros(0, 1);  % Empty column vector
    all_polarity = [];
    
    for seg_idx = 1:length(segments)
        seg_start = segments{seg_idx}.start;
        seg_end = segments{seg_idx}.end;
        seg_sig = sig(seg_start:seg_end);
        seg_tm = tm(seg_start:seg_end);
        
        fprintf('\n--- Processing segment %d/%d ---\n', seg_idx, length(segments));
        
        % Skip if segment too short
        if length(seg_sig) < 2*fs
            fprintf(' - Skipping (too short: %.1f s)\n', length(seg_sig)/fs);
            continue
        end
        
        % Process this segment with P&T algorithm
        [seg_peaks, seg_pol] = process_ecg_segment(seg_sig, seg_tm, fs, peakThresh, search_back, ref_period);
        
        % Adjust peak indices back to original signal
        if ~isempty(seg_peaks)
            seg_peaks = seg_peaks(:);  % Ensure column vector
            all_peaks = [all_peaks; seg_peaks + seg_start - 1];
            all_polarity = [all_polarity; seg_pol];
        end
    end
    
    % Sort combined peaks chronologically
    [peaks, ~] = sort(all_peaks);
    
    % Overall polarity based on most common across segments
    if ~isempty(all_polarity)
        pol_mode = mode(all_polarity);
        sign = pol_mode * median(abs(sig(peaks)));
    else
        sign = [];
        pol_mode = 1;
    end
    
    fprintf('\n=== COMBINED RESULTS ===\n');
    fprintf(' - Total R-peaks detected: %d\n', length(peaks));
    if ~isempty(all_polarity)
        unique_pols = unique(all_polarity);
        if length(unique_pols) > 1
            fprintf(' - NOTE: Polarity varied across segments\n');
        end
        if pol_mode < 0
            fprintf(" - Overall peaks' polarity: negative\n");
        else
            fprintf(" - Overall peaks' polarity: positive\n");
        end
    end

    % RR intervals and HR
    RR = diff(peaks) ./ fs;
    RR_t = tm(peaks);
    HR = 60 ./ diff(RR_t);

    if ~isempty(RR)
        fprintf(' - Median RR: %.3f s (%.1f bpm)\n', median(RR), 60/median(RR));
    end
    

%%  PPG
elseif strcmpi(sig_type, 'ppg')

    if ~exist('fs','var')
        error("You must provide your signal' sampling rate as 2nd input")
    end

    % PARAMETERS
    if isfield(params,'ppg_buffer')
        BUFLN = params.ppg_buffer;
    else
        BUFLN = 4096;
    end

    if isfield(params,'ppg_learnperiod')
        LPERIOD = fs*params.ppg_learnperiod;
    else
        LPERIOD  = fs*5;
    end

    if isfield(params,'ppg_learnthresh')
        minthresh = params.ppg_learnthresh;
    else
        minthresh = 5;
    end

    if isfield(params,'ppg_eyeclosing')
        EyeClosing =  round(fs*params.ppg_eyeclosing);
    else
        EyeClosing = round(fs*0.65);
    end

    if isfield(params,'ppg_expctperiod')
        ExpectPeriod = round(fs*params.ppg_expctperiod);
    else
        ExpectPeriod = round(fs*5);
    end

    if isfield(params,'ppg_slopewindow')
        SLPwindow = round(fs*params.ppg_slopewindow);
    else
        SLPwindow = round(fs*0.1);
    end

    INVALID_signal = -32758;

    timer = 0;
    peaks = [];
    beat_n = 1;
    from = 1;
    to = length(signal);

    if signal(1) <= INVALID_signal
        signal(1) = mean(signal);
    end
    inv = find(signal<=INVALID_signal);
    for i = 1:length(inv)
        signal(inv(i)) = signal(inv(i)-1);
    end

    if length(signal) < 5*60*fs
        signal = (signal-min(signal))./(max(signal)-min(signal)).*4000-2000;
    else
        n = 1;
        for i=1:5*60*fs:length(signal)
            max_signal(n)=max(signal(i:min(i+5*60*fs-1,length(signal))));
            min_signal(n)=min(signal(i:min(i+5*60*fs-1,length(signal))));
            n=n+1;
        end
        signal = (signal-median(min_signal))./(median(max_signal)-median(min_signal)).*4000-2000;
    end

    ebuf(1:BUFLN) = 0;
    lbuf = ebuf;
    if from > BUFLN
        tt_2 = from-BUFLN;
    else
        tt_2 = 0;
    end

    t1 = 8*fs;
    t1 = t1+from;
    T0 = 0;
    n = 0;
    for t = from:t1
        [temp,ebuf,lbuf,tt_2] = slpsamp(t,signal,BUFLN,ebuf,lbuf,tt_2,SLPwindow);
        if temp > INVALID_signal
            T0 = T0+temp;
            n=n+1;
        end
    end
    T0 = T0/n;
    Ta = 3*T0;

    learning = 1;

    t = from;
    while t <= to

        if learning
            if t > from + LPERIOD
                learning = 0;
                T1 = T0;
                t = from;
            else
                T1 = 2*T0;
            end
        end

        [temp,ebuf,lbuf,tt_2] = slpsamp(t,signal,BUFLN,ebuf,lbuf,tt_2,SLPwindow);

        if temp > T1
            timer = 0;
            maxd = temp;
            mind = maxd;
            tmax = t;
            for tt = t + 1: t + EyeClosing-1
                [temp2 ,ebuf,lbuf,tt_2] = slpsamp(tt,signal,BUFLN,ebuf,lbuf,tt_2,SLPwindow);
                if temp2 > maxd
                    maxd=temp2;
                    tmax=tt;
                end
            end
            if maxd == temp
                t = t+1;
                continue
            end

            for tt = tmax :-1: t-EyeClosing/2+1
                [temp2 ,ebuf,lbuf,tt_2] = slpsamp(tt,signal,BUFLN,ebuf,lbuf,tt_2,SLPwindow);
                if temp2< mind
                    mind=temp2;
                end
            end
            if maxd > mind+10
                onset = (maxd-mind)/100+2;
                tpq = t-round(0.04*fs);
                maxmin_2_3_threshold=(maxd-mind)*2.0/3;
                for tt = tmax:-1:t-EyeClosing/2+1
                    [temp2, ebuf,lbuf,tt_2] = slpsamp(tt,signal,BUFLN,ebuf,lbuf,tt_2,SLPwindow);
                    if temp2 < maxmin_2_3_threshold
                        break
                    end
                end
                for tt = tt:-1:t - EyeClosing / 2 + round(0.024*fs)
                    [temp2 ,ebuf,lbuf,tt_2] = slpsamp(tt,signal,BUFLN,ebuf,lbuf,tt_2,SLPwindow);
                    [temp3 ,ebuf,lbuf,tt_2] = slpsamp(tt-round(0.024*fs),signal,BUFLN,ebuf,lbuf,tt_2,SLPwindow);
                    if temp2 - temp3<onset
                        tpq = tt-round(0.016*fs);
                        break
                    end
                end

                valley_v = round(tpq);
                for valley_i = round(max(2,tpq-round(0.20*fs))):round(min(tpq+round(0.05*fs),length(signal)-1))

                    if valley_v <= 0
                        t = t + 1;
                        continue;
                    end

                    if signal(valley_v)>signal(valley_i) && signal(valley_i)<=signal(valley_i-1) && signal(valley_i)<=signal(valley_i+1)
                        valley_v = valley_i;
                    end
                end


                if ~learning

                    if beat_n == 1
                        if round(valley_v) > 0
                            peaks(beat_n) = round(valley_v);
                            beat_n = beat_n + 1;
                        end
                    else
                        if round(valley_v) > peaks(beat_n-1)
                            peaks(beat_n) = round(valley_v);
                            beat_n = beat_n + 1;
                        end
                    end
                end

                Ta = Ta + (maxd - Ta)/10;
                T1 = Ta / 3;

                t = tpq+EyeClosing;
            end
        else
            if ~learning
                timer = timer+1;
                if timer > ExpectPeriod && Ta > minthresh
                    Ta = Ta-1;
                    T1 = Ta / 3;
                end
            end
        end
        t=t+1;
    end

    sig = signal;
    RR = diff(peaks) ./ fs;
    RR_t = peaks ./ fs;
    HR = 60 ./ diff(tm(peaks));

else
    error('Signal type must be ECG or PPG')
end

end

%% Helper function: Automatic boundary detection
function boundary_indices = detect_signal_boundaries(sig, fs)
    % Automatically detect session boundaries based on amplitude and baseline changes
    % Similar to how EEGLAB detects discontinuities for filtering
    % Returns indices where boundaries should be placed
    
    % Parameters
    window_size = round(3 * fs);  % 3-second window (increased from 2)
    min_segment_duration = 30 * fs;  % Minimum 30 seconds between boundaries (increased from 10)
    threshold_factor = 5;  % Sensitivity (increased from 3.5 - higher = less sensitive)
    
    n = length(sig);
    boundary_indices = [];
    
    % Need sufficient data
    if n < 2 * min_segment_duration
        return
    end
    
    % Compute rolling statistics on absolute signal
    sig_abs = abs(sig);
    half_win = floor(window_size/2);
    
    % Preallocate
    rolling_mean = zeros(n, 1);
    rolling_std = zeros(n, 1);
    
    for i = 1:n
        win_start = max(1, i - half_win);
        win_end = min(n, i + half_win);
        window_data = sig_abs(win_start:win_end);
        rolling_mean(i) = mean(window_data);
        rolling_std(i) = std(window_data);
    end
    
    % Smooth to reduce noise
    smooth_window = round(fs);  % 1 second
    rolling_mean_smooth = movmean(rolling_mean, smooth_window);
    rolling_std_smooth = movmean(rolling_std, smooth_window);
    
    % Detect large changes
    mean_changes = abs(diff(rolling_mean_smooth));
    std_changes = abs(diff(rolling_std_smooth));
    
    % Normalize by median (robust to outliers)
    mean_change_norm = mean_changes / (median(mean_changes) + eps);
    std_change_norm = std_changes / (median(std_changes) + eps);
    
    % Combined change metric
    combined_change = max(mean_change_norm, std_change_norm);
    
    % Find peaks in change metric
    [~, locs] = findpeaks(combined_change, ...
        'MinPeakHeight', threshold_factor, ...
        'MinPeakDistance', min_segment_duration);
    
    % Adjust for diff operation (add 1)
    candidate_boundaries = locs + 1;
    
    % Additional validation: check if there's actually a significant change
    for i = 1:length(candidate_boundaries)
        idx = candidate_boundaries(i);
        
        % Compare 10-second windows before and after (increased from 5 for stability)
        before_win = max(1, idx - 10*fs):idx-1;
        after_win = idx:min(n, idx + 10*fs - 1);
        
        if ~isempty(before_win) && ~isempty(after_win)
            before_amp = mean(abs(sig(before_win)));
            after_amp = mean(abs(sig(after_win)));
            
            % Significant amplitude change (>100% difference, i.e., 2x ratio)
            amp_ratio = max(before_amp, after_amp) / (min(before_amp, after_amp) + eps);
            
            if amp_ratio > 2.0  % Increased from 1.5 to be less sensitive
                boundary_indices = [boundary_indices; idx];
            end
        end
    end
    
    % Sort and ensure unique
    boundary_indices = unique(boundary_indices);
end

%% Helper function: Process single ECG segment
function [peaks, pol] = process_ecg_segment(sig, tm, fs, peakThresh, search_back, ref_period)
    % Process a single continuous ECG segment using P&T algorithm
    % This contains the core detection logic
    
    % Ensure column vectors
    sig = sig(:);
    tm = tm(:);
    nSamp = length(sig);
    
    % Constants
    med_smooth_nb_coef = round(fs/100);
    int_nb_coef = round(7*fs/256);
    
    % P&T operations (zero-phase)
    dffecg = [0; diff(sig)];
    sqrecg = dffecg.^2;
    
    % Integrate using filtfilt
    b_int = ones(int_nb_coef,1) / int_nb_coef;
    intecg = filtfilt(b_int, 1, sqrecg);
    
    % Median filter (odd length)
    if mod(med_smooth_nb_coef, 2) == 0
        med_smooth_nb_coef = med_smooth_nb_coef + 1;
    end
    mdfint = medfilt1(intecg, med_smooth_nb_coef);
    
    % Match length to sig
    if numel(mdfint) < numel(sig)
        mdfint = [mdfint; zeros(numel(sig) - numel(mdfint), 1)];
    elseif numel(mdfint) > numel(sig)
        mdfint = mdfint(1:numel(sig));
    end
    
    % P&T threshold
    if nSamp/fs > 90
        xs = sort(mdfint(fs:fs*90));
    else
        xs = sort(mdfint(fs:end));
    end
    
    if nSamp/fs > 10
        ind_xs = ceil(98/100*length(xs));
        en_thres = xs(ind_xs);
    else
        ind_xs = ceil(99/100*length(xs));
        en_thres = xs(ind_xs);
    end
    
    % Candidate regions
    poss_reg = mdfint > (peakThresh * en_thres);
    if ~any(poss_reg)
        peaks = [];
        pol = 1;
        fprintf('   No candidate regions found\n');
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
        fprintf('   No peaks after search-back\n');
        return
    end
    
    % Extract features for polarity detection (temporal ordering method)
    idxMax = zeros(1, nb_peaks);
    idxMin = zeros(1, nb_peaks);
    valMax = zeros(1, nb_peaks);
    valMin = zeros(1, nb_peaks);
    
    for i = 1:nb_peaks
        a = left(i);
        b = right(i);
        seg = sig(a:b);
        
        [valMax(i), imx] = max(seg);
        [valMin(i), imn] = min(seg);
        
        idxMax(i) = a + imx - 1;
        idxMin(i) = a + imn - 1;
    end
    
    % Polarity detection: which extremum comes first?
    nb_pol = min(nb_peaks, max(3, round(30*fs/median(diff(left)))));
    polVotes = ones(1, nb_pol);
    for i = 1:nb_pol
        if idxMin(i) < idxMax(i)
            polVotes(i) = -1;
        else
            polVotes(i) = +1;
        end
    end
    pol = median(polVotes);
    if pol >= 0
        pol = +1;
    else
        pol = -1;
    end
    
    % Peak selection with polarity constraint
    peaks = zeros(nb_peaks, 1);  % Column vector
    pkval = zeros(nb_peaks, 1);  % Column vector
    
    for i = 1:nb_peaks
        a = left(i);
        b = right(i);
        
        if pol > 0
            % R expected positive: search before min (S)
            stop = min(idxMin(i), b);
            if stop <= a + 1
                stop = b;
            end
            [pkval(i), ii] = max(sig(a:stop));
            peaks(i) = a + ii - 1;
            
            if peaks(i) > idxMin(i) && idxMax(i) >= a && idxMax(i) <= b
                peaks(i) = idxMax(i);
                pkval(i) = valMax(i);
            end
        else
            % R expected negative: search before max
            stop = min(idxMax(i), b);
            if stop <= a + 1
                stop = b;
            end
            [pkval(i), ii] = min(sig(a:stop));
            peaks(i) = a + ii - 1;
            
            if peaks(i) > idxMax(i) && idxMin(i) >= a && idxMin(i) <= b
                peaks(i) = idxMin(i);
                pkval(i) = valMin(i);
            end
        end
    end
    
    % Enforce refractory period
    [peaks, ord] = sort(peaks, 'ascend');
    pkval = pkval(ord);
    
    keep = true(size(peaks));
    for k = 2:numel(peaks)
        if (peaks(k) - peaks(k-1)) < round(ref_period * fs)
            if abs(pkval(k)) > abs(pkval(k-1))
                keep(k-1) = false;
            else
                keep(k) = false;
            end
        end
    end
    peaks = peaks(keep);
    
    % Micro-refinement
    Lsig = numel(sig);
    micro = max(1, round(0.015 * fs));
    
    for k = 1:numel(peaks)
        p = peaks(k);
        w1 = max(1, p - micro);
        w2 = min(Lsig, p + micro);
        
        if pol > 0
            [~, ii] = max(sig(w1:w2));
        else
            [~, ii] = min(sig(w1:w2));
        end
        peaks(k) = w1 + ii - 1;
    end
    
    peaks = unique(peaks, 'stable');
    peaks = peaks(:);  % Ensure column vector
    
    fprintf('   Detected %d peaks (polarity: %s, threshold: %.2f)\n', ...
        length(peaks), iif(pol==1, 'pos', 'neg'), en_thres);
end

% Simple inline if helper
function out = iif(condition, true_val, false_val)
    if condition
        out = true_val;
    else
        out = false_val;
    end
end

%% PPG Subfunction
function [beat1,ebuf,lbuf,tt_2] = slpsamp(t,signal,BUFLN,ebuf,lbuf,tt_2, SLPwindow)

while t > tt_2
    prevVal = 0;

    if tt_2>0 && tt_2-1>0 && tt_2<length(signal) && tt_2-1<length(signal)
        val2 = signal(tt_2 - 1);
        val1 = signal(tt_2);
    else
        val2 = prevVal;
        val1 = val2;
    end

    dy =  val1-val2;
    if dy < 0
        dy = 0;
    end
    tt_2 = tt_2+1;
    M = round(mod(tt_2,(BUFLN-1))+1);
    et = dy;
    ebuf(M) = et;
    aet = 0;
    for i = 0:SLPwindow-1
        p = M-i;
        if p <= 0
            p = p+BUFLN;
        end
        aet = aet+ebuf(p);
    end
    lbuf(M) = aet;

end
M3 = round(mod(t,(BUFLN-1))+1);
beat1 = lbuf(M3);

end