% Detect R-peaks from ECG signals or pulse wave onsets from PPG signals.
%
% ECG method (Pan-Tompkins):
%   Optional zero-phase FIR highpass filter for baseline wander removal.
%   QRS detection via differentiation, squaring, and moving-average
%   integration (Hann window, zero-phase). Energy threshold at 98th
%   percentile. Search-back algorithm recovers missed beats when RR
%   variability exceeds 1.5x the median. Polarity is determined from
%   QRS ordering over the first 40 beats. Each coarse peak is refined
%   to the nearest local extremum within a 15 ms window.
%
% PPG method:
%   Zero-phase FIR bandpass filter (0.5-3 Hz). Detects valleys (default)
%   or peaks using MATLAB's findpeaks with adaptive MinPeakHeight (MAD-
%   based) and MinPeakDistance (trimmed mean of initial RR estimates).
%
% Usage:
%   [RR, RR_t, peaks, signal, times, sign, HR] = get_RR(signal, times, params)
%
% Inputs:
%   signal  - raw ECG or PPG signal (1 x N or N x 1)
%   times   - time vector in milliseconds (1 x N or N x 1)
%   params  - struct with fields:
%               .fs             - sampling rate (Hz)
%               .heart_signal   - 'ecg' or 'ppg'
%
%             ECG optional fields:
%               .ecg_highpass   - highpass cutoff (Hz); skipped if absent or 0
%               .ecg_peakthresh - P&T energy threshold multiplier (default 0.6)
%               .ecg_searchback - enable search-back (default true)
%               .ecg_refperiod  - refractory period in s (default 0.25)
%
%             PPG optional fields:
%               .ppg_detect_mode - 'valleys' (default) or 'peaks'
%               .ppg_height_method - MinPeakHeight method: 'mad' (default),
%                                    'std', or 'trimmean'
%
% Outputs:
%   RR      - RR intervals (s), length N-1
%   RR_t    - timestamps of RR intervals (s), length N-1
%   peaks   - R-peak or pulse-onset sample indices, length N-1
%   signal  - signal after any filtering applied inside this function
%   times   - time vector in seconds
%   polarity- median peak amplitude (ECG polarity indicator; [] for PPG)
%   HR      - heart rate (bpm), length N-2
%
% Notes:
%   - All outputs have the first element removed to match RR interval length.
%   - For HEP analysis, epoch the output signal (already filtered if
%     params.ecg_highpass is set) using the output peaks as triggers.
%
% References:
%   Pan & Tompkins (1985). IEEE Trans Biomed Eng.
%   Vest et al. (2018). Physiological Measurement.
%
% % Changelog:
%   v2.0 - April 2026 (Cedric Cannard)
%     ECG:
%       - Highpass filter now applied by default (params.ecg_highpass = 0.5 Hz)
%         for baseline wander removal. Set params.ecg_highpass = 0 to skip
%         if signal is already pre-filtered externally.
%       - Replaced causal integration filter with zero-phase filtfilt
%         (Hann window) to eliminate phase delay in QRS detection.
%       - Removed medfilt1 smoothing step which was distorting signal morphology.
%       - Fixed polarity detection: voting now over first 40 beats (was
%         incorrectly using 30*fs segments instead of beats).
%       - Fixed HR output length: now matches RR (length N-1) instead of N-2.
%       - Peak refinement now correctly uses the pre-filtered signal throughout.
%       - Output variable renamed from 'sign' to 'polarity' to avoid
%         shadowing MATLAB's built-in sign() function.
%     PPG:
%       - Completely replaced the legacy Physionet sliding-buffer algorithm
%         (qppg) with a simpler, more transparent findpeaks-based detector.
%       - Zero-phase FIR bandpass filter applied (0.5-3 Hz) before detection.
%       - Adaptive MinPeakHeight estimated via robust statistics (MAD-based
%         by default; configurable via params.ppg_height_method).
%       - Adaptive MinPeakDistance estimated from median of initial RR.
%       - Detection mode (valleys or peaks) configurable via params.ppg_detect_mode.
%     General:
%       - Removed dead code (max_force, stale length-match block, slpsamp).
%       - HR now always computed before first-element removal for consistency.
%       - Improved input validation and error/warning messages throughout.
% 
% Copyright (C), BrainBeats, Cedric Cannard, 2023

function [RR, RR_t, peaks, signal, times, polarity, HR] = get_RR(signal, times, params)

fs       = params.fs;
sig_type = params.heart_signal;

% Sanity check: times should be in milliseconds
% If max value < 1000, it's likely already in seconds
if max(times) < 1000 && max(times) > 0
    error('get_RR: times appears to be in seconds (max=%.3f). Pass times in milliseconds (e.g. EEG.times).', max(times));
end

% Enforce column vector and convert time to seconds
signal = signal(:);
times  = times(:) / 1000;
nSamp  = numel(signal);

if numel(times) ~= nSamp
    error('get_RR: times must have the same number of samples as signal.');
end


polarity = [];

%=========================================================================
%  ECG - Pan-Tompkins QRS detector
%=========================================================================
if strcmpi(sig_type, 'ecg')

    % Parameters
    peakThresh  = 0.6;   if isfield(params,'ecg_peakthresh'), peakThresh  = params.ecg_peakthresh; end
    search_back = true;  if isfield(params,'ecg_searchback'),  search_back = params.ecg_searchback; end
    ref_period  = 0.25;  if isfield(params,'ecg_refperiod'),   ref_period  = params.ecg_refperiod;  end

    % Bandpass filter for ECG (default: [3 25] Hz)
    % Set params.ecg_bandpass = false to skip if pre-filtered externally
    if isfield(params, 'ecg_bandpass') && ~params.ecg_bandpass
        fprintf('  Bandpass filter: skipped (pre-filtered externally)\n');
    else
        % Check for non-finite values before filtering
        non_finite = ~isfinite(signal);
        if any(non_finite)
            n_bad = sum(non_finite);
            warning('  Warning: %d non-finite samples detected (%.2f%%) — interpolating before filtering!!!\n', ...
                n_bad, 100 * n_bad / numel(signal));
            t = (1:numel(signal))';   % column vector
            signal = interp1(t(~non_finite), signal(~non_finite), t, 'pchip', 'extrap');
        end

        bp_cutoff = [3 25];  % default [highpass lowpass] in Hz
        if isfield(params, 'ecg_bandpass') && numel(params.ecg_bandpass) == 2
            bp_cutoff = params.ecg_bandpass;
        end
    
        % Highpass
        hp_order = 3 * round(fs / bp_cutoff(1));
        b_hp     = fir1(hp_order, bp_cutoff(1) / (fs/2), 'high');
        signal   = filtfilt(b_hp, 1, signal);
    
        % Lowpass
        lp_order = 5 * round(fs / bp_cutoff(2));
        b_lp     = fir1(lp_order, bp_cutoff(2) / (fs/2), 'low');
        signal   = filtfilt(b_lp, 1, signal);
        fprintf('  Bandpass filter: %.1f–%.1f Hz (order-%d/%d FIR, zero-phase)\n', ...
            bp_cutoff(1), bp_cutoff(2), hp_order, lp_order);
    end

    % Flatline check
    if prctile(abs(signal), 95) < 0.05
        error('get_RR: ECG amplitude too small - likely a flat line.');
    end

    % P&T pipeline: differentiate, square, integrate
    dffecg      = [0; diff(signal)];
    sqrecg      = dffecg .^ 2;
    % int_nb_coef = round(7 * fs / 256); % 7 samples at 256 Hz is only ~27 ms, which is too narrow. The validated window is 150 ms
    int_nb_coef = round(0.150 * fs);  % 150 ms, more standard
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
                medRRv        = median(RRv);
                missedIdx     = find(diff(times(indAT)) > 1.5 * medRRv);
                indStart      = indAT(missedIdx);
                indEnd        = indAT(missedIdx + 1);
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

    % Per-segment extrema
    idxMax = zeros(1, nb_peaks);  idxMin = zeros(1, nb_peaks);
    valMax = zeros(1, nb_peaks);  valMin = zeros(1, nb_peaks);
    for i = 1:nb_peaks
        seg = signal(left(i):right(i));
        [valMax(i), imx] = max(seg);  idxMax(i) = left(i) + imx - 1;
        [valMin(i), imn] = min(seg);  idxMin(i) = left(i) + imn - 1;
    end

    % Polarity from first 40 beats (max before min = positive R)
    nb_pol   = min(nb_peaks, 40);
    polVotes = arrayfun(@(i) 2*(idxMax(i) < idxMin(i)) - 1, 1:nb_pol);
    pol      = sign(median(polVotes));
    if pol == 0, pol = 1; end
    fprintf('  Peaks polarity: %s\n', ternary(pol > 0, 'positive', 'negative'));

    % Peak localization with polarity constraint
    peaks = zeros(1, nb_peaks);
    pkval = zeros(1, nb_peaks);
    for i = 1:nb_peaks
        a = left(i);  b = right(i);
        if pol > 0
            stop = max(a+1, min(idxMin(i), b));
            [pkval(i), ii] = max(signal(a:stop));
            peaks(i) = a + ii - 1;
            if peaks(i) > idxMin(i) && idxMax(i) >= a && idxMax(i) <= b
                peaks(i) = idxMax(i);  pkval(i) = valMax(i);
            end
        else
            stop = max(a+1, min(idxMax(i), b));
            [pkval(i), ii] = min(signal(a:stop));
            peaks(i) = a + ii - 1;
            if peaks(i) > idxMax(i) && idxMin(i) >= a && idxMin(i) <= b
                peaks(i) = idxMin(i);  pkval(i) = valMin(i);
            end
        end
    end
    polarity = median(pkval);

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

    % Micro-refinement: snap to nearest local extremum within 15 ms
    micro = max(1, round(0.015 * fs));
    for k = 1:numel(peaks)
        w1 = max(1, peaks(k) - micro);
        w2 = min(nSamp, peaks(k) + micro);
        if pol > 0, [~, ii] = max(signal(w1:w2));
        else,       [~, ii] = min(signal(w1:w2));
        end
        peaks(k) = w1 + ii - 1;
    end
    peaks = unique(peaks, 'stable');

    fprintf('  P&T energy threshold: %.2f\n', en_thres);

    % RR intervals and HR
    RR     = diff(peaks) / fs;
    RR_t   = times(peaks);
    HR     = 60 ./ RR;
    RR_t(1) = [];
    peaks(1) = [];

% =========================================================================
%  PPG - findpeaks with adaptive threshold and distance
% =========================================================================
elseif strcmpi(sig_type, 'ppg')

    % Parameters
    detect_mode   = 'valleys';  if isfield(params,'ppg_detect_mode'),    detect_mode   = params.ppg_detect_mode;    end
    height_method = 'mad';      if isfield(params,'ppg_height_method'),  height_method = params.ppg_height_method;  end

    % Zero-phase FIR highpass then lowpass (skipped if pre-filtered externally)
    if ~isfield(params, 'ppg_bandpass') || params.ppg_bandpass
        % Highpass: 0.5 Hz
        hp_order = 2 * round(fs / 0.5);
        b_hp     = fir1(hp_order, 0.5 / (fs/2), 'high');
        signal   = filtfilt(b_hp, 1, signal);
        % Lowpass: 3 Hz
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
    tmp_rr       = diff(tmp_peaks) / fs;
    % minPeakDist  = round(trimmean(tmp_rr, 20) * fs / 2);
    minPeakDist  = round(median(tmp_rr) * fs / 2);

    % Final detection
    [~, peaks] = findpeaks(det_sig, 'MinPeakHeight', minPeakHeight, 'MinPeakDistance', minPeakDist);
    fprintf('  PPG detection mode: %s\n', detect_mode);
    fprintf('  MinPeakHeight: %.2f | MinPeakDistance: %d samples (%.2f s)\n', ...
        minPeakHeight, minPeakDist, minPeakDist/fs);

    % RR intervals and HR
    RR     = diff(peaks) / fs;
    RR_t   = times(peaks);
    HR     = 60 ./ RR;
    RR_t(1)  = [];
    peaks(1) = [];

else
    error('get_RR: params.heart_signal must be ''ecg'' or ''ppg''.');
end

end  % main function

%% =========================================================================
%  Helpers: adaptive MinPeakHeight estimation
% =========================================================================
function minPeakHeight = estimateMinPeakHeight(sig, method)
% Estimate minimum peak height for findpeaks using robust statistics.
%
% Inputs:
%   sig    - signal vector (already flipped for valley detection if needed)
%   method - 'mad' (default), 'std', 'percentile', or 'trimmean'
%
% Output:
%   minPeakHeight - scalar threshold

switch lower(method)
    case 'mad'
        % minPeakHeight = trimmean(sig, 20) + 0.5 * mad(sig);
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