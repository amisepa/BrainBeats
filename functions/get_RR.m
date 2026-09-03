% Detect R-peaks from ECG signals or pulse wave onsets from PPG signals.
%
% ECG method (Pan-Tompkins):
%   Zero-phase FIR bandpass filter (default 3-35 Hz) for QRS isolation.
%   QRS detection via differentiation, squaring, and moving-average
%   integration (Hann window, 150 ms, zero-phase). Energy threshold at
%   98th percentile. Search-back recovers missed beats when an RR interval
%   exceeds 1.5x the local median. Polarity is determined by a single
%   global vote across all detected QRS segments: for each segment the
%   dominant deflection amplitude (positive vs negative) casts a vote,
%   and the majority sets the global polarity (medical standard, consistent
%   with Pan-Tompkins). An adaptive chunked mode is available for
%   concatenated recordings with DC shifts between sessions (see
%   params.ecg_adaptive_pol). Each coarse peak is refined to the nearest
%   local extremum within a 15 ms window.
%
% PPG method:
%   Zero-phase FIR bandpass filter (0.5-3 Hz). Detects valleys (default)
%   or peaks using MATLAB's findpeaks with adaptive MinPeakHeight (MAD-
%   based) and MinPeakDistance (from median of initial RR estimates).
%
% Usage:
%   [RR, RR_t, peaks, signal, times, polarity, HR] = get_RR(signal, times, params)
%
% Inputs:
%   signal  - raw ECG or PPG signal (1 x N or N x 1)
%   times   - time vector in milliseconds (1 x N or N x 1)
%   params  - struct with required and optional fields (see below)
%
% Required params fields:
%   .fs           - sampling rate (Hz)
%   .heart_signal - signal type: 'ecg' or 'ppg'
%
% Optional ECG params fields:
%   .ecg_bandpass      - [hp lp] bandpass cutoffs in Hz (default: [3 35])
%                        set to false to skip filtering (if pre-filtered externally)
%   .ecg_peakthresh    - P&T energy threshold multiplier (default: 0.35)
%                        lower = more sensitive, higher = more conservative
%   .ecg_searchback    - enable search-back for missed beats (default: true)
%   .ecg_refperiod     - refractory period in s; beats closer than this are
%                        deduplicated by amplitude (default: 0.25 s)
%
%   /!\ THESE TWO DEFAULTS WERE TRANSPOSED IN THIS HEADER UNTIL 2026-08-28.
%       It read peakthresh 0.25 / refperiod 0.35, while the code below has
%       always assigned peakthresh 0.35 (line ~182) and refperiod 0.25
%       (line ~183). The CODE is what ran, so 0.35 / 0.25 are the values
%       behind every peak train this function has ever produced. Corrected
%       in the header, deliberately -- NOT in the code.
%
%   /!\ DO NOT "RESTORE" refperiod TO 0.35 s WITHOUT RE-VALIDATING. On a slow
%       heart the T-wave can fall BETWEEN the two values -- e.g. a recording
%       with the T at ~348 ms after R, inside 0.35 s and outside 0.25 s -- so
%       0.35 s would dedupe a real beat against its own T-wave and delete one
%       of the pair. It would also shift peaks on ordinary recordings and
%       invalidate any validation against existing trains. Note this parameter
%       is NOT what keeps T-waves out: the Pan-Tompkins energy front end
%       (differentiate, square, 150 ms Hann integration, 98th-percentile
%       threshold) leaves the low-slope T below threshold, which is why get_RR
%       is robust where gradient-based detectors are not. Lowering peakthresh
%       to 0.25 would likewise admit more candidates, not fewer.
%   .ecg_polarity      - force polarity: 1 (upright R-wave) or -1 (inverted lead)
%                        overrides all automatic detection (default: [], auto)
%   .ecg_adaptive_pol  - use adaptive chunked polarity instead of global QRS vote
%                        (default: false). Useful for concatenated recordings
%                        with DC shifts between sessions. When true, polarity is
%                        estimated in escalating chunk sizes (nSamp/4 -> nSamp/2
%                        -> full recording) until chunks are consistent (<30%
%                        minority); each P&T segment inherits its chunk polarity.
%
% Optional PPG params fields:
%   .ppg_bandpass      - apply bandpass filter (default: true)
%                        set to false if pre-filtered externally
%   .ppg_detect_mode   - 'valleys' (default) or 'peaks'
%   .ppg_height_method - MinPeakHeight estimation method: 'mad' (default),
%                        'std', or 'trimmean'
%
% Outputs:
%   RR       - RR intervals (s), length N-1
%   RR_t     - timestamps of RR interval onsets (s), length N-1
%   peaks    - R-peak or pulse-onset sample indices, length N-1
%   signal   - signal after bandpass filtering (or raw if filtering skipped)
%   times    - time vector converted to seconds
%   polarity - median peak amplitude; sign indicates ECG lead polarity
%              (positive = upright R-wave, negative = inverted lead); [] for PPG
%   HR       - instantaneous heart rate (bpm), length N-1
%
% Notes:
%   - The first element of RR, RR_t, peaks, and HR is removed on output
%     so that each RR interval is paired with the timestamp of its
%     trailing peak (consistent with HRV convention).
%   - For HEP analysis: epoch the output signal using output peaks as
%     triggers. The signal returned is bandpass-filtered (3-35 Hz by
%     default), suitable for direct use as HEP trigger input.
%   - To skip bandpass filtering (e.g. signal already preprocessed):
%       params.ecg_bandpass = false;
%
% Examples:
%   % Basic ECG usage (medical-standard polarity, auto-detected):
%   params.fs           = 256;
%   params.heart_signal = 'ecg';
%   [RR, RR_t, peaks, ecg_filt, times_s, polarity, HR] = get_RR(ecg, EEG.times, params);
%
%   % Force polarity manually (e.g. known inverted lead):
%   params.ecg_polarity = -1;
%   [RR, RR_t, peaks, ecg_filt, times_s, polarity, HR] = get_RR(ecg, EEG.times, params);
%
%   % Use adaptive chunked polarity for concatenated sessions:
%   params.ecg_adaptive_pol = true;
%   [RR, RR_t, peaks, ecg_filt, times_s, polarity, HR] = get_RR(ecg, EEG.times, params);
%
%   % PPG with default valley detection:
%   params.fs           = 256;
%   params.heart_signal = 'ppg';
%   [RR, RR_t, peaks, ppg_filt, times_s, ~, HR] = get_RR(ppg, EEG.times, params);
%
%   % PPG with peak detection and custom bandpass:
%   params.ppg_detect_mode = 'peaks';
%   params.ppg_bandpass    = false;   % already filtered externally
%   [RR, RR_t, peaks, ppg_filt, times_s, ~, HR] = get_RR(ppg, EEG.times, params);
%
% References:
%   Pan & Tompkins (1985). IEEE Trans Biomed Eng, 32(3), 230-236.
%   Vest et al. (2018). Physiological Measurement, 39(10).
%
% Changelog:
%   v2.3 - June 2026 
%     ECG:
%       - Replaced adaptive chunked polarity (default) with medical-standard
%         global QRS-segment vote, consistent with Pan-Tompkins (1985).
%         For each detected QRS window, positive vs negative peak amplitude
%         determines a local vote; the majority across all beats sets the
%         global polarity. Operates on bandpass-filtered QRS windows, making
%         it robust to baseline wander and T/P waves.
%       - Adaptive chunked polarity preserved as opt-in via
%         params.ecg_adaptive_pol = true (recommended only for concatenated
%         recordings with DC shifts between sessions).
%       - Added params.ecg_polarity for manual polarity override (+1/-1),
%         taking priority over all automatic detection.
%       - Updated header documentation and added usage examples.
%
%   v2.2 - April 2026
%     ECG:
%       - Replaced single global polarity vote with adaptive chunked polarity:
%         starts at nSamp/4s chunks, escalates to nSamp/2 then full recording
%         if polarity is inconsistent across chunks (minority fraction >30%).
%       - prctile(99/1) used within each chunk for spike robustness.
%       - Each P&T segment inherits the polarity of its chunk; polarity
%         can flip at chunk boundaries within one recording.
%       - Micro-refinement uses per-beat chunk polarity throughout.
%
%   v2.1 - April 2026 
%     ECG:
%       - Fixed polarity detection: replaced fragile QRS-ordinal vote
%         (idxMax < idxMin) with global amplitude comparison (|max| vs
%         |min| over full bandpass-filtered signal).
%       - Fixed peak localization for negative polarity: search now spans
%         the full segment [left..right].
%       - Simplified peak localization loop.
%
%   v2.0 - April 2026 (Cedric Cannard)
%       - Replaced ecg_highpass with ecg_bandpass ([hp lp] Hz, default [3 35]).
%       - Replaced causal integration filter with zero-phase filtfilt (Hann window).
%       - Removed medfilt1 smoothing step.
%       - Fixed HR output length.
%       - Output variable renamed from 'sign' to 'polarity'.
%     PPG:
%       - Replaced legacy Physionet qppg with findpeaks-based detector.
%       - Zero-phase FIR bandpass filter (0.5-3 Hz).
%       - Adaptive MinPeakHeight (MAD-based by default).
%
% Copyright (C), BrainBeats, Cedric Cannard, 2023

function [RR, RR_t, peaks, signal, times, polarity, HR] = get_RR(signal, times, params)

fs       = params.fs;
sig_type = params.heart_signal;

% Sanity check: times should be in milliseconds
if max(times) < 1000 && max(times) > 0
    error('get_RR: times appears to be in seconds (max=%.3f). Pass times in milliseconds (e.g. EEG.times).', max(times));
end

% Enforce column vector and convert time to seconds
signal = signal(:);
times  = times(:) ./ 1000;
nSamp  = numel(signal);

if numel(times) ~= nSamp
    error('get_RR: times must have the same number of samples as signal.');
end

polarity = [];
%=========================================================================
%  ECG - Pan-Tompkins QRS detector
%=========================================================================
if strcmpi(sig_type, 'ecg')

    % --- Default parameters -------------------------------------------
    bp_cutoff        = [3 35];   % bandpass [hp lp] Hz
    peakThresh       = 0.35;      % P&T energy threshold multiplier
    search_back      = true;     % search-back for missed beats
    ref_period       = 0.25;     % refractory period (s)
    ecg_polarity     = [];       % [] = auto; 1 = force positive; -1 = force negative
    adaptive_pol     = false;    % false = medical-standard global vote (default)
                                 % true  = adaptive chunked polarity (legacy)

    if isfield(params, 'ecg_bandpass')     && numel(params.ecg_bandpass) == 2, bp_cutoff   = params.ecg_bandpass;  end
    if isfield(params, 'ecg_peakthresh'),    peakThresh  = params.ecg_peakthresh; end
    if isfield(params, 'ecg_searchback'),    search_back = params.ecg_searchback; end
    if isfield(params, 'ecg_refperiod'),     ref_period  = params.ecg_refperiod;  end
    if isfield(params, 'ecg_adaptive_pol'),  adaptive_pol = params.ecg_adaptive_pol; end
    if isfield(params, 'ecg_polarity') && ~isempty(params.ecg_polarity)
        ecg_polarity = sign(params.ecg_polarity);
    end
    % ------------------------------------------------------------------

    % Bandpass filter
    if isfield(params, 'ecg_bandpass') && isequal(params.ecg_bandpass, false)
        fprintf('  Bandpass filter: skipped (pre-filtered externally)\n');
    else
        % Interpolate non-finite samples before filtering
        non_finite = ~isfinite(signal);
        if any(non_finite)
            n_bad = sum(non_finite);
            warning('get_RR: %d non-finite samples (%.2f%%) - interpolating before filtering.', ...
                n_bad, 100 * n_bad / numel(signal));
            t      = (1:numel(signal))';
            signal = interp1(t(~non_finite), signal(~non_finite), t, 'pchip', 'extrap');
        end

        hp_order = 3 * round(fs / bp_cutoff(1));
        b_hp     = fir1(hp_order, bp_cutoff(1) / (fs/2), 'high');
        signal   = filtfilt(b_hp, 1, signal);

        lp_order = 5 * round(fs / bp_cutoff(2));
        b_lp     = fir1(lp_order, bp_cutoff(2) / (fs/2), 'low');
        signal   = filtfilt(b_lp, 1, signal);
        fprintf('  Bandpass filter: %.1f-%.1f Hz (order-%d/%d FIR, zero-phase)\n', ...
            bp_cutoff(1), bp_cutoff(2), hp_order, lp_order);
    end

    % Flatline check
    if prctile(abs(signal), 95) < 0.05
        error('get_RR: ECG amplitude too small - likely a flat line.');
    end

    % P&T pipeline: differentiate, square, integrate
    dffecg      = [0; diff(signal)];
    sqrecg      = dffecg .^ 2;
    int_nb_coef = round(0.150 * fs);  % 150 ms standard window
    b_int       = hann(int_nb_coef) / sum(hann(int_nb_coef));
    signal_filt = filtfilt(b_int, 1, sqrecg);
    fprintf('  P&T integration filter: lowpass ~%.1f Hz (Hann window, zero-phase)\n', fs / int_nb_coef);

    % Energy threshold (98th percentile, skip first second)
    xs        = sort(signal_filt(fs : min(fs*90, nSamp)));
    en_thres  = xs(ceil(0.98 * numel(xs)));
    poss_reg  = signal_filt > (peakThresh * en_thres);

    if ~any(poss_reg)
        warning('get_RR: no QRS candidates found - check signal quality.');
        [RR, RR_t, peaks, HR] = deal([]);
        return
    end

    % Search-back for missed beats
    if search_back
        indAT = find(poss_reg);
        if numel(indAT) > 2
            RRv = diff(times(indAT));
            RRv = RRv(RRv > 0.01);
            if ~isempty(RRv)
                medRRv    = median(RRv);
                missedIdx = find(diff(times(indAT)) > 1.5 * medRRv);
                indStart  = indAT(missedIdx);
                indEnd    = indAT(missedIdx + 1);
                for i = 1:numel(indStart)
                    poss_reg(indStart(i):indEnd(i)) = ...
                        signal_filt(indStart(i):indEnd(i)) > (0.5 * peakThresh * en_thres);
                end
            end
        end
    end

    % Segment boundaries
    left     = find(diff([0; poss_reg]) ==  1);
    right    = find(diff([poss_reg; 0]) == -1);
    nb_peaks = numel(left);

    if nb_peaks == 0
        warning('get_RR: no QRS segments found after thresholding.');
        [RR, RR_t, peaks, HR] = deal([]);
        return
    end

    % -----------------------------------------------------------------
    % Polarity detection
    %
    % Priority order:
    %   1. params.ecg_polarity: manual override (+1 or -1), skips all detection.
    %   2. params.ecg_adaptive_pol = true: adaptive chunked polarity (legacy
    %      mode, useful for concatenated sessions with DC shifts between them).
    %   3. Default (medical standard, P&T): single global vote across all
    %      detected QRS segments. For each segment, the dominant deflection
    %      (positive vs negative peak amplitude) determines the local vote;
    %      the majority across all beats sets the global polarity. This is
    %      robust to baseline wander and T/P waves because it operates on
    %      the bandpass-filtered QRS windows themselves, not the raw signal.
    % -----------------------------------------------------------------

    if ~isempty(ecg_polarity)
        % 1. Manual override
        pol_seg = repmat(ecg_polarity, 1, nb_peaks);
        fprintf('  Polarity: forced %s (params.ecg_polarity = %d)\n', ...
            ternary(ecg_polarity > 0, 'positive', 'negative'), ecg_polarity);

    elseif adaptive_pol
        % 2. Adaptive chunked polarity (legacy)
        chunk_sizes   = [floor(nSamp/fs/4), floor(nSamp/fs/2), Inf];
        inconsist_thr = 0.3;
        chunk_pol     = [];
        chunk_n       = 0;
        used_chunk_s  = NaN;

        for cs = chunk_sizes
            if isinf(cs)
                chunk_n  = nSamp;
                n_chunks = 1;
            else
                chunk_n  = round(cs * fs);
                n_chunks = max(1, floor(nSamp / chunk_n));
            end

            chunk_pol = zeros(1, n_chunks);
            for c = 1:n_chunks
                i1  = (c-1) * chunk_n + 1;
                i2  = min(c * chunk_n, nSamp);
                seg = signal(i1:i2);
                seg_hi = prctile(seg, 99);
                seg_lo = prctile(seg,  1);
                chunk_pol(c) = sign(abs(seg_hi) - abs(seg_lo));
                if chunk_pol(c) == 0, chunk_pol(c) = 1; end
            end

            n_minority    = sum(chunk_pol ~= mode(chunk_pol));
            minority_frac = n_minority / n_chunks;

            if isinf(cs)
                used_chunk_s = Inf;
                break;
            elseif minority_frac <= inconsist_thr
                used_chunk_s = cs;
                break;
            else
                fprintf('  Polarity inconsistent at %ds chunks (%.0f%% minority) - escalating\n', ...
                    cs, 100 * minority_frac);
            end
        end

        n_pos_chunks = sum(chunk_pol > 0);
        n_neg_chunks = sum(chunk_pol < 0);
        if isinf(used_chunk_s)
            fprintf('  Chunk polarity (full recording): %d/%d positive, %d/%d negative\n', ...
                n_pos_chunks, n_chunks, n_neg_chunks, n_chunks);
        else
            fprintf('  Chunk polarity (%ds): %d/%d positive, %d/%d negative\n', ...
                used_chunk_s, n_pos_chunks, n_chunks, n_neg_chunks, n_chunks);
        end

        seg_mid = round((left + right) / 2);
        pol_seg = ones(1, nb_peaks);
        for c = 1:n_chunks
            i1 = (c-1) * chunk_n + 1;
            i2 = min(c * chunk_n, nSamp);
            in_chunk = seg_mid >= i1 & seg_mid <= i2;
            pol_seg(in_chunk) = chunk_pol(c);
        end

    else
        % 3. Medical-standard global vote (default): one vote per QRS segment
        % % Based on amplitude:
        votes = zeros(1, nb_peaks);
        for i = 1:nb_peaks
            seg       = signal(left(i):right(i));
            peak_pos  = max(seg);
            peak_neg  = abs(min(seg));
            votes(i)  = sign(peak_pos - peak_neg);
            if votes(i) == 0, votes(i) = 1; end
        end
        pol_global = sign(sum(votes));
        if pol_global == 0, pol_global = 1; end
        pol_seg    = repmat(pol_global, 1, nb_peaks);
        fprintf('  Polarity: %s (global QRS vote, %d/%d segments)\n', ...
            ternary(pol_global > 0, 'positive', 'negative'), ...
            sum(votes == pol_global), nb_peaks);

        % Sometimes S is bigger than R, leading to incorrect R-peak detections.
        % R always comes before S regardless of their relative amplitudes, 
        % so timing can be more robust than amplitude comparison:
        % votes = zeros(1, nb_peaks);
        % for i = 1:nb_peaks
        %     seg          = signal(left(i):right(i));
        %     [~, idx_pos] = max(seg);
        %     [~, idx_neg] = min(seg);
        %     % Upright QRS: positive R-peak precedes negative S-trough
        %     votes(i) = sign(idx_neg - idx_pos);
        %     if votes(i) == 0, votes(i) = 1; end
        % end
        % pol_global = sign(sum(votes));
        % if pol_global == 0, pol_global = 1; end
        % pol_seg    = repmat(pol_global, 1, nb_peaks);
        % fprintf('  Polarity: %s (timing-based QRS vote, %d/%d segments)\n', ...
        %     ternary(pol_global > 0, 'positive', 'negative'), ...
        %     sum(votes == pol_global), nb_peaks);
    end

    % Auto-flip inverted signal so R-peaks are always positive.
    % Applies to all polarity modes. Disable with params.ecg_flip_signal = false.
    flip_signal = true;
    if isfield(params, 'ecg_flip_signal'), flip_signal = params.ecg_flip_signal; end

    pol_majority = sign(sum(pol_seg));
    if pol_majority == 0, pol_majority = 1; end

    if pol_majority < 0 && flip_signal
        signal   = -signal;
        pol_seg(:) = 1;
        fprintf('  Signal flipped: inverted lead corrected to positive polarity\n');
    end
    
    % Peak localization: max or min of each P&T segment per its polarity
    peaks = zeros(1, nb_peaks);
    pkval = zeros(1, nb_peaks);
    for i = 1:nb_peaks
        a = left(i);  b = right(i);
        if pol_seg(i) > 0
            [pkval(i), ii] = max(signal(a:b));
        else
            [pkval(i), ii] = min(signal(a:b));
        end
        peaks(i) = a + ii - 1;
    end

    polarity = median(pkval);
    pol      = sign(median(pol_seg));  % majority - used as fallback only
    if pol == 0, pol = 1; end

    % Refractory period: keep larger |amplitude| when peaks too close
    [peaks, ord] = sort(peaks, 'ascend');
    pkval = pkval(ord);
    keep  = true(size(peaks));
    for k = 2:numel(peaks)
        if peaks(k) - peaks(k-1) < round(ref_period * fs)
            if abs(pkval(k)) > abs(pkval(k-1)), keep(k-1) = false;
            else,                                keep(k)   = false;
            end
        end
    end
    peaks = peaks(keep);

    % Micro-refinement: snap to nearest local extremum within 15 ms.
    % Uses pol_seg (per-segment polarity) aligned to the kept peaks.
    % After refractory filtering peaks is a subset - reindex pol_seg to match.
    pol_seg_kept = pol_seg(ord);
    pol_seg_kept = pol_seg_kept(keep);
    micro = max(1, round(0.015 * fs));
    for k = 1:numel(peaks)
        w1 = max(1, peaks(k) - micro);
        w2 = min(nSamp, peaks(k) + micro);
        if pol_seg_kept(k) > 0, [~, ii] = max(signal(w1:w2));
        else,                    [~, ii] = min(signal(w1:w2));
        end
        peaks(k) = w1 + ii - 1;
    end
    peaks = unique(peaks, 'stable');

    fprintf('  P&T energy threshold: %.2f\n', en_thres);

    % RR intervals and HR
    RR      = diff(peaks) / fs;
    RR_t    = times(peaks);
    HR      = 60 ./ RR;
    RR_t(1)  = [];
    peaks(1) = [];



% =========================================================================
%  PPG - findpeaks with adaptive threshold and distance
% =========================================================================
elseif strcmpi(sig_type, 'ppg')

    % --- Default parameters -------------------------------------------
    detect_mode   = 'valleys';  % 'valleys' or 'peaks'
    height_method = 'mad';      % MinPeakHeight method: 'mad', 'std', 'trimmean'

    if isfield(params, 'ppg_detect_mode'),   detect_mode   = params.ppg_detect_mode;   end
    if isfield(params, 'ppg_height_method'), height_method = params.ppg_height_method; end
    % ------------------------------------------------------------------

    % Zero-phase FIR bandpass (skipped if pre-filtered externally)
    if ~isfield(params, 'ppg_bandpass') || params.ppg_bandpass
        hp_order = 2 * round(fs / 0.5);
        b_hp     = fir1(hp_order, 0.5 / (fs/2), 'high');
        signal   = filtfilt(b_hp, 1, signal);
        lp_order = 2 * round(fs / 3);
        b_lp     = fir1(lp_order, 3 / (fs/2), 'low');
        signal   = filtfilt(b_lp, 1, signal);
        fprintf('  Highpass: 0.5 Hz + Lowpass: 3 Hz (FIR, zero-phase)\n');
    else
        fprintf('  Filtering: skipped (pre-filtered externally)\n');
    end

    % Flip signal for valley detection
    if strcmpi(detect_mode, 'valleys')
        det_sig = -signal;
    else
        det_sig = signal;
    end

    % Adaptive MinPeakHeight
    minPeakHeight = estimateMinPeakHeight(det_sig, height_method);

    % Estimate MinPeakDistance from initial detection
    [~, tmp_peaks] = findpeaks(det_sig, 'MinPeakHeight', minPeakHeight);
    if numel(tmp_peaks) < 2
        warning('get_RR: too few PPG peaks detected - check signal quality or detection mode.');
        [RR, RR_t, peaks, HR] = deal([]);
        return
    end
    tmp_rr      = diff(tmp_peaks) / fs;
    minPeakDist = round(median(tmp_rr) * fs / 2);

    % Final detection
    [~, peaks] = findpeaks(det_sig, 'MinPeakHeight', minPeakHeight, 'MinPeakDistance', minPeakDist);
    fprintf('  PPG detection mode: %s\n', detect_mode);
    fprintf('  MinPeakHeight: %.2f | MinPeakDistance: %d samples (%.2f s)\n', ...
        minPeakHeight, minPeakDist, minPeakDist/fs);

    % RR intervals and HR
    RR      = diff(peaks) / fs;
    RR_t    = times(peaks);
    HR      = 60 ./ RR;
    RR_t(1)  = [];
    peaks(1) = [];

else
    error('get_RR: params.heart_signal must be ''ecg'' or ''ppg''.');
end

end  % main function

%% =========================================================================
%  Helpers
% =========================================================================
function minPeakHeight = estimateMinPeakHeight(sig, method)
switch lower(method)
    case 'mad'
        minPeakHeight = median(sig, 'omitnan') + 0.5 * mad(sig);
    case 'std'
        minPeakHeight = median(sig, 'omitnan') + 0.5 * std(sig, 'omitnan');
    case 'trimmean'
        pos = sig(sig > 0);
        if isempty(pos), pos = sig; end
        minPeakHeight = trimmean(pos, 20);
    otherwise
        error('estimateMinPeakHeight: unknown method ''%s''.', method);
end
minPeakHeight = round(minPeakHeight, 2);
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end