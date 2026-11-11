% GET_RR - Detect R-peaks (ECG) or pulse-wave onsets (PPG) and return RR intervals.
%
% Usage:
%   [RR, RR_t, peaks, signal, times, polarity, HR] = get_RR(signal, times, params)
%
% Inputs:
%   signal - raw ECG or PPG signal, 1 x N or N x 1
%   times  - time vector in ms, same length as signal (e.g. EEG.times);
%            an error is thrown if it looks like seconds
%   params - struct. Required fields:
%     .fs           - sampling rate (Hz)
%     .heart_signal - 'ecg' or 'ppg'
%   Optional ECG fields:
%     .ecg_bandpass     - [hp lp] FIR bandpass in Hz (default [3 35]); false to skip
%     .ecg_peakthresh   - energy threshold multiplier; lower = more sensitive (default 0.35)
%     .ecg_searchback   - re-search long gaps at half threshold (default true)
%     .ecg_refperiod    - refractory period (s); of two closer peaks the larger
%                         is kept (default 0.25, see Notes)
%     .ecg_polarity     - force polarity, 1 or -1 (default [] = automatic)
%     .ecg_adaptive_pol - per-chunk polarity, for concatenated recordings with
%                         DC shifts (default false = one global QRS vote)
%     .ecg_flip_signal  - flip an inverted lead so R-peaks are positive (default true)
%   Optional PPG fields:
%     .ppg_polarity      - 1 (blood volume: systolic peak up) or -1 (light
%                          intensity, inverted); default [] = automatic (see
%                          Method)
%     .ppg_bandpass      - 0.5-3 Hz FIR bandpass (default true); false to skip
%     .ppg_detect_mode   - 'valleys' (pulse onsets, default) or 'peaks'
%                          (systolic peaks)
%     .ppg_height_method - 'prominence' (default): each pulse must stand out by
%                          30% of the local (2-s) pulse amplitude, so beats are
%                          kept where the pulse amplitude drops; or a global
%                          MinPeakHeight: 'mad', 'std', 'trimmean' (older behavior)
%
% Outputs (the first detected beat is dropped, so each interval is paired
% with the peak that closes it):
%   RR       - RR intervals (s)
%   RR_t     - time of the peak closing each interval (s), to a fraction of
%              a sample (see Method)
%   peaks    - sample index of the peak closing each interval
%   signal   - filtered signal, column (ECG or PPG flipped if inverted); use
%              it with peaks for HEP epoching
%   times    - input time vector converted to s, column
%   polarity - ECG: median amplitude of the detected peaks in the returned
%              signal (positive once an inverted lead is flipped); PPG: 1 or
%              -1 (inverted input, flipped)
%   HR       - instantaneous heart rate (bpm), 60 ./ RR
%
% Method:
%   ECG: Pan-Tompkins. Bandpass, derivative, squaring and 150 ms Hann
%   integration (all zero-phase). Candidates exceed ecg_peakthresh x the
%   98th percentile of the integrated energy (1-90 s); gaps > 1.5x the
%   median spacing are searched again. Each peak is the extremum of its
%   segment, refined within +/-15 ms. RR and RR_t use the peak times to a
%   fraction of a sample (vertex of the parabola through the peak and its two
%   neighbors); peaks stay sample indices. On the sample grid, RR intervals
%   were off by 1.65 ms RMS at 250 Hz and 3.2 ms at 128 Hz (64 ds003838
%   ECGs resampled from 1000 Hz), which biased RMSSD (+0.35 ms at 128 Hz) and
%   pNN50 (up to 3 points per recording); with the parabola 0.04 and 0.19 ms
%   (MIT-BIH, from 360 Hz to 180 and 128 Hz: 2.29 -> 0.11, 3.22 -> 0.24 ms;
%   tests/validate_subsample_rr.m).
%   PPG: polarity from the quantile skewness of the slope (5th and 95th
%   percentiles, 0.5-8 Hz; the systolic upstroke is the fast edge; a warning
%   below 0.2), inverted signals flipped; then findpeaks on
%   the filtered signal (negated for valleys), with
%   MinPeakDistance = half the median interval of a first pass, and a
%   prominence of at least 30% of the local peak-to-peak amplitude. A global
%   height threshold missed low-amplitude pulses (first 3 s of the sample
%   PPG: 301 beats vs 306 R-peaks; 306 with the local prominence). Pulse
%   times are refined with the same parabola (61 ds003838 PPGs resampled
%   from 1000 Hz: interval error 1.87 -> 0.90 ms RMS at 250 Hz, 3.30 -> 0.64
%   ms at 128 Hz; RMSSD bias +0.13 and +0.40 ms -> 0).
%
% Notes:
%   Shortest recordings: ECG 3 s, PPG 12 s (filtfilt needs more than three
%   lengths of the high-pass filters; shorter ones stop with an error).
%   ecg_refperiod is 0.25 s, not 0.35 s: on a slow heart a T-wave detected
%   ~350 ms after R would be deduplicated against its own beat and one of
%   the pair deleted; 0.25 s is also the value all validation used.
%   T-waves are kept out by the energy threshold, not by this parameter.
%
% References:
%   Pan & Tompkins (1985). IEEE Trans Biomed Eng, 32(3), 230-236.
%   Vest et al. (2018). Physiological Measurement, 39(10).
%
% Copyright (C) - Cedric Cannard, 2023

function [RR, RR_t, peaks, signal, times, polarity, HR] = get_RR(signal, times, params)

fs       = params.fs;
sig_type = params.heart_signal;

% Sanity check: times should be in milliseconds, one sample every 1000/fs ms
% (the former test, max(times) < 1000, also refused recordings under 1 s)
if numel(times) > 1 && abs(median(diff(times(:)))*fs - 1) < 0.01
    error('get_RR: times appears to be in seconds (one sample every %.4f). Pass times in milliseconds (e.g. EEG.times).', median(diff(times(:))));
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

    % Flatline check, in any unit (uV, mV, V, raw amplifier units): the raw
    % signal barely varies around its offset. (After the band-pass, a
    % constant keeps a fraction of its offset: the FIR high-pass does not
    % reach zero gain at 0 Hz.)
    ok = isfinite(signal);
    if ~any(ok)
        error('get_RR: the ECG holds no finite sample (all NaN or Inf).');
    end
    dev = prctile(abs(signal(ok) - median(signal(ok))), 95);
    if isempty(dev) || dev <= 1e-5 * abs(median(signal(ok)))
        error('get_RR: ECG amplitude too small - likely a flat line.');
    end

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
        check_length(numel(signal), b_hp, fs, 'ECG', sprintf('%g-Hz high-pass', bp_cutoff(1)), 'ecg_bandpass');
        signal   = filtfilt(b_hp, 1, signal);

        lp_order = 5 * round(fs / bp_cutoff(2));
        b_lp     = fir1(lp_order, bp_cutoff(2) / (fs/2), 'low');
        check_length(numel(signal), b_lp, fs, 'ECG', sprintf('%g-Hz low-pass', bp_cutoff(2)), 'ecg_bandpass');
        signal   = filtfilt(b_lp, 1, signal);
        fprintf('  Bandpass filter: %.1f-%.1f Hz (order-%d/%d FIR, zero-phase)\n', ...
            bp_cutoff(1), bp_cutoff(2), hp_order, lp_order);
    end

    % P&T pipeline: differentiate, square, integrate
    dffecg      = [0; diff(signal)];
    sqrecg      = dffecg .^ 2;
    int_nb_coef = round(0.150 * fs);  % 150 ms standard window
    b_int       = hann(int_nb_coef) / sum(hann(int_nb_coef));
    signal_filt = filtfilt(b_int, 1, sqrecg);
    fprintf('  P&T integration filter: lowpass ~%.1f Hz (Hann window, zero-phase)\n', fs / int_nb_coef);

    % Energy threshold: 98th percentile of the integrated energy between
    % 1 s (skips the filter edge) and 90 s
    xs        = sort(signal_filt(round(fs) : min(round(fs*90), nSamp)));
    en_thres  = xs(ceil(0.98 * numel(xs)));
    poss_reg  = signal_filt > (peakThresh * en_thres);

    if ~any(poss_reg)
        warning('get_RR: no QRS candidates found - check signal quality.');
        [RR, RR_t, peaks, HR] = deal([]);
        return
    end

    % Search-back for missed beats: gaps between candidates longer than
    % 1.5x the median spacing are re-thresholded at half the level
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
    %   3. Default: single global vote across all detected QRS segments.
    %      Each segment votes for its dominant deflection (positive vs
    %      negative amplitude); the majority sets the polarity. Working on
    %      the filtered QRS windows makes it robust to baseline wander and
    %      T/P waves.
    % -----------------------------------------------------------------

    if ~isempty(ecg_polarity)
        % 1. Manual override
        pol_seg = repmat(ecg_polarity, 1, nb_peaks);
        fprintf('  Polarity: forced %s (params.ecg_polarity = %d)\n', ...
            ternary(ecg_polarity > 0, 'positive', 'negative'), ecg_polarity);

    elseif adaptive_pol
        % 2. Adaptive chunked polarity (legacy)
        chunk_sizes   = [floor(nSamp/fs/4), floor(nSamp/fs/2), Inf];
        chunk_sizes   = chunk_sizes(chunk_sizes >= 1);  % < 4 s of data gives 0-s chunks: use the whole recording
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
        % 3. Global amplitude vote (default): one vote per QRS segment.
        % Leads where S is deeper than R can vote negative; force the
        % polarity with params.ecg_polarity in that case.
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
    end

    % Flip an inverted lead so R-peaks are positive (all polarity modes).
    % Disable with params.ecg_flip_signal = false.
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
    pol      = sign(median(pol_seg));  % majority polarity (not used below)
    if pol == 0, pol = 1; end

    % Refractory period: of two peaks closer than ref_period, keep the larger |amplitude|
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

    % Micro-refinement: snap to the local extremum within +/-15 ms, using
    % the polarity of each kept peak's segment (pol_seg reindexed by sort
    % order and refractory mask).
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

    % RR intervals and HR from the sub-sample peak times (peaks stay sample
    % indices); drop the first peak so each interval is paired with the peak
    % that closes it
    RR_t    = subsample_times(signal, peaks, times, fs);
    RR      = diff(RR_t)';
    HR      = 60 ./ RR;
    RR_t(1)  = [];
    peaks(1) = [];

% =========================================================================
%  PPG - findpeaks with adaptive threshold and distance
% =========================================================================
elseif strcmpi(sig_type, 'ppg')

    % --- Default parameters -------------------------------------------
    detect_mode   = 'valleys';  % 'valleys' or 'peaks'
    height_method = 'prominence'; % or a global MinPeakHeight: 'mad', 'std', 'trimmean'

    if isfield(params, 'ppg_detect_mode'),   detect_mode   = params.ppg_detect_mode;   end
    if isfield(params, 'ppg_height_method'), height_method = params.ppg_height_method; end
    % ------------------------------------------------------------------

    % Polarity: in a blood-volume PPG the systolic upstroke is the fast edge,
    % so the slope (0.5-8 Hz) has a longer positive than negative tail. Raw
    % device outputs are often light intensity, inverted (all 70 PPGs of known
    % polarity we have: 61 ds003838 recordings and the sample dataset, by the
    % order of pulse foot and systolic peak after the ECG R-peak, and 8 Muse
    % optics channels): they are flipped, so that 'valleys' are the pulse
    % onsets and 'peaks' the systolic peaks whatever the device. Without the
    % flip, 'valleys' marked the systolic peak of these recordings (R + 384-424
    % ms instead of R + 92-112 ms), and the PPG-locked HEP matched the
    % ECG-locked one less well (r .52-.86 instead of .68-.92;
    % tests/compare_ppg_transit.m). The tails are compared at the 5th and 95th
    % percentiles of the slope (quantile skewness, -1 to 1): right on all 70,
    % clean and with injected steps, motion (up to 10% of the recording),
    % dropouts or clipping. The skewness used until 2026-11-10 was decided by
    % a few large transients: right in 52-67% with the larger artifacts, and
    % wrong on one clean PPG; the 1st and 99th percentiles in 85-88% with heavy
    % motion (tests/validate_ppg_polarity.m). Clean PPGs gave |q| >= 0.46,
    % three with scattered pulses (sub-052, 056, 080) |q| < 0.03: below 0.2,
    % a warning.
    if isfield(params, 'ppg_polarity') && ~isempty(params.ppg_polarity)
        polarity = sign(params.ppg_polarity);
        fprintf('  PPG polarity: forced %+d (params.ppg_polarity)\n', polarity);
    else
        [bP, aP] = butter(2, [0.5 min(8, 0.45*fs)]/(fs/2));
        q = quantile(diff(filtfilt(bP, aP, signal - mean(signal))), [0.05 0.5 0.95]);
        slopeSkew = (q(3) + q(1) - 2*q(2)) / (q(3) - q(1));
        polarity = 1 - 2*(slopeSkew < 0);
        fprintf('  PPG polarity: quantile skewness of the slope %+.2f (%s)\n', slopeSkew, ...
            ternary(polarity > 0, 'blood volume', 'light intensity, inverted'));
        if ~(abs(slopeSkew) >= 0.2)
            warning(['get_RR: PPG polarity unclear (quantile skewness of the slope %+.2f; ' ...
                'clean PPGs beyond +/-0.46): the PPG may hold few pulses. Check it, or set ''ppg_polarity''.'], slopeSkew);
        end
    end
    if polarity < 0
        signal = -signal;
        fprintf('  Signal flipped: PPG in blood-volume convention (systolic peak up)\n');
    end

    % Zero-phase FIR bandpass (skipped if pre-filtered externally)
    if ~isfield(params, 'ppg_bandpass') || params.ppg_bandpass
        hp_order = 2 * round(fs / 0.5);
        b_hp     = fir1(hp_order, 0.5 / (fs/2), 'high');
        check_length(numel(signal), b_hp, fs, 'PPG', '0.5-Hz high-pass', 'ppg_bandpass');
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

    if strcmpi(height_method, 'prominence')
        % Prominence relative to the local pulse amplitude (peak-to-peak over
        % 2 s), robust to slow amplitude changes (respiration, vasomotion)
        [~, tmp_peaks] = findpeaks(det_sig, 'MinPeakHeight', 0);
        if numel(tmp_peaks) < 2
            warning('get_RR: too few PPG peaks detected - check signal quality or detection mode.');
            [RR, RR_t, peaks, HR] = deal([]);
            return
        end
        minPeakDist = round(median(diff(tmp_peaks)) / 2);
        localRange = movmax(det_sig, round(2*fs)) - movmin(det_sig, round(2*fs));
        [~, peaks, ~, prom] = findpeaks(det_sig, 'MinPeakDistance', minPeakDist);
        peaks = peaks(prom >= 0.3 * localRange(peaks));
        fprintf('  PPG detection mode: %s | MinPeakDistance: %d samples (%.2f s) | prominence >= 30%% of the local amplitude\n', ...
            detect_mode, minPeakDist, minPeakDist/fs);
    else
        % Global MinPeakHeight (older behavior)
        minPeakHeight = estimateMinPeakHeight(det_sig, height_method);

        % MinPeakDistance = half the median interval of a first detection pass
        [~, tmp_peaks] = findpeaks(det_sig, 'MinPeakHeight', minPeakHeight);
        if numel(tmp_peaks) < 2
            warning('get_RR: too few PPG peaks detected - check signal quality or detection mode.');
            [RR, RR_t, peaks, HR] = deal([]);
            return
        end
        minPeakDist = round(median(diff(tmp_peaks)) / 2);

        [~, peaks] = findpeaks(det_sig, 'MinPeakHeight', minPeakHeight, 'MinPeakDistance', minPeakDist);
        fprintf('  PPG detection mode: %s | MinPeakHeight: %.2f | MinPeakDistance: %d samples (%.2f s)\n', ...
            detect_mode, minPeakHeight, minPeakDist, minPeakDist/fs);
    end

    % RR intervals and HR from the sub-sample pulse times (first peak
    % dropped, as for ECG)
    RR_t    = subsample_times(signal, peaks, times, fs);
    RR      = diff(RR_t);
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

function t = subsample_times(sig, pk, times, fs)
% Times (s, column) of the peaks or troughs pk of sig, refined to a
% fraction of a sample: vertex of the parabola through each extremum and its
% two neighbors (|shift| <= half a sample; samples at the edges, or that are
% not a local extremum, keep their sample time)
pk = pk(:);
t = times(pk);
in = pk > 1 & pk < numel(sig);
y = sig(pk(in) + (-1:1));
if size(y,1) ~= nnz(in), y = y.'; end      % a single peak gives a row
d = 0.5 * (y(:,1) - y(:,3)) ./ (y(:,1) - 2*y(:,2) + y(:,3));
d(~isfinite(d) | abs(d) > 0.5) = 0;
t(in) = t(in) + d/fs;
end

function check_length(n, b, fs, what, filt, opt)
% filtfilt needs more than 3 filter lengths of data: stop with the minimum
% duration instead of filtfilt's own message (ECG: 3 s with the default
% 3-Hz high-pass; PPG: 12 s with the 0.5-Hz high-pass)
need = 3*(numel(b) - 1);
if n <= need
    error(['get_RR: the %s is %.1f s long, and its %s filter needs more than %.1f s. ' ...
        'Use a longer recording, or a band-passed signal with ''%s'' set to false.'], what, n/fs, filt, need/fs, opt);
end
end
