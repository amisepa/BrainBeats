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
%   params      - params structure containing the following fields:
%                   fs (sample rate) and heart_signal ('ecg' or 'ppg')
%
% OUTPUTS:
%   RR          - RR intervals
%   RR_t        - time vector
%   peaks      - R peaks
%   sig         - ECG signal filtered
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
% Copyright (C), BrainBeats, Cedric Cannard, 2023

function [RR, RR_t, peaks, sig, tm, sign, HR] = get_RR(signal, params)

% Parameters
fs = params.fs;
sig_type = params.heart_signal;

if size(signal,1) < size(signal,2)
    signal = signal';
end

sign = [];
nSamp = size(signal,1);
tm = 1/fs:1/fs:nSamp/fs;

%% ECG
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
        ref_period = 0.25; % refractory period
    end

    % Constants
    med_smooth_nb_coef = round(fs/100);     % scales with fs
    int_nb_coef = round(7*fs/256);          % length is 7 for fs = 256 Hz

    % Use a consistent orientation: COLUMN vectors everywhere in ECG branch
    sig = signal(:);

    % If median of 20% of the samples are < min_amp, abort (flat line)
    min_amp = 0.1;
    if length(find(abs(sig)<min_amp)) / nSamp > 0.20
        error('ECG time series is a flat line')
    end

    % P&T operations (zero-phase where applicable)
    dffecg = diff(sig);                 % differentiate (one datum shorter)
    sqrecg = dffecg .* dffecg;          % square

    % integrate using zero-phase filtering (filtfilt)
    b_int = ones(1,int_nb_coef) / int_nb_coef;
    intecg = filtfilt(b_int, 1, sqrecg);

    % smooth using median filter (odd length)
    if mod(med_smooth_nb_coef, 2) == 0
        med_smooth_nb_coef = med_smooth_nb_coef + 1;
    end
    mdfint = medfilt1(intecg, med_smooth_nb_coef);

    % match length to sig
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

    max_force = []; % to force the energy threshold value
    if isempty(max_force)
        if nSamp/fs > 10
            ind_xs = ceil(98/100*length(xs));
            en_thres = xs(ind_xs);
        else
            ind_xs = ceil(99/100*length(xs));
            en_thres = xs(ind_xs);
        end
    else
        en_thres = max_force;
    end

    % candidate regions (energy)
    poss_reg = mdfint > (peakThresh * en_thres);
    if isempty(poss_reg)
        poss_reg(10) = true;
    end

    % search-back for missed beats
    if search_back
        indAboveThreshold = find(poss_reg);
        RRv = diff(tm(indAboveThreshold));
        medRRv = median(RRv(RRv > 0.01));
        indMissedBeat = find(RRv > 1.5*medRRv);

        indStart = indAboveThreshold(indMissedBeat);
        indEnd   = indAboveThreshold(indMissedBeat+1);

        for i = 1:length(indStart)
            poss_reg(indStart(i):indEnd(i)) = ...
                mdfint(indStart(i):indEnd(i)) > (0.5 * peakThresh * en_thres);
        end
    end


    % segment boundaries
    left  = find(diff([0; poss_reg])==1);
    right = find(diff([poss_reg; 0])==-1);

    nb_peaks = numel(left);
    if nb_peaks == 0
        peaks = [];
        RR = [];
        RR_t = [];
        HR = [];
        return
    end

    % ------------------------------------------------------------
    % FIX Jnauary 15, 2026: R-peak localization by QRS ordering, not amplitude
    % Inside each QRS segment, compute both max and min of raw ECG.
    % The R peak is the first dominant extremum of the QRS.
    % This prevents selecting S when S is the larger absolute deflection.
    % ------------------------------------------------------------

    % 1) Per-segment extrema and their order
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

    % 2) Decide global polarity from ordering across early beats
    % If max tends to come before min, treat R as positive, else negative.
    nb_pol = min(nb_peaks, max(3, round(30*fs))); % up to ~30 s worth of segments
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

    % 3) Pick peaks with polarity, but constrain search to occur before opposite extremum
    peaks = zeros(1, nb_peaks);
    pkval = zeros(1, nb_peaks);

    for i = 1:nb_peaks
        a = left(i);
        b = right(i);

        if pol > 0
            % R expected positive: force search before the min (S) when possible
            stop = min(idxMin(i), b);
            if stop <= a + 1
                stop = b;
            end
            [pkval(i), ii] = max(sig(a:stop));
            peaks(i) = a + ii - 1;

            % if the chosen max still occurs after the min (edge case), fall back to earliest max
            if peaks(i) > idxMin(i) && idxMax(i) >= a && idxMax(i) <= b
                peaks(i) = idxMax(i);
                pkval(i) = valMax(i);
            end
        else
            % R expected negative: force search before the max when possible
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

    % representative polarity value for reporting
    sign = median(pkval);

    % 4) Enforce refractory period (keep larger |amplitude| within ref period)
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

    % 5) Micro-refine around each peak to nearest local extremum of the chosen polarity
    Lsig  = numel(sig);
    micro = max(1, round(0.015 * fs)); % 15 ms

    for k = 1:numel(peaks)
        p  = peaks(k);
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

    if pol < 0
        fprintf(" - Peaks' polarity: negative \n");
    else
        fprintf(" - Peaks' polarity: positive \n");
    end
    fprintf(' - P&T energy threshold: %g \n', round(en_thres,2))

    % RR intervals and HR
    RR   = diff(peaks) ./ fs;
    RR_t = tm(peaks);
    HR   = 60 ./ diff(RR_t);

    %%  PPG
elseif strcmpi(sig_type, 'ppg')

    if ~exist('fs','var')
        error("You must provide your signal' sampling rate as 2nd input")
    end

    % PARAMETERS

    % Length of the buffer BUFLN to store a segment of the PPG signal for processing.
    % The length of this buffer (4096 samples by default) determines how much of the signal is held in memory for analysis at any given time.
    % A sufficiently large buffer size ensures that the algorithm has enough data to accurately detect pulse waves,
    % but it also must be balanced with computational efficiency. 4096 samples is a good compromise between these factors.
    if isfield(params,'ppg_buffer')
        BUFLN = params.ppg_buffer;
    else
        BUFLN = 4096; % must be a power of 2 (default = 4096).
    end

    % LPERIOD allows the algorithm to tailor its detection strategy based on the
    % initial segment of the PPG signal, improving the accuracy and reliability of pulse detection,
    % especially in the context of varying signal qualities or individual  differences.
    if isfield(params,'ppg_learnperiod')
        LPERIOD = fs*params.ppg_learnperiod;
    else
        LPERIOD  = fs*5;   % learning period in samples (default = 5 s).
    end

    % Minimum threshold value for the detection algorithm.
    % It serves as a baseline or lower limit for the algorithm to identify a pulse wave in the PPG signal.
    % The presence of a minimum threshold helps in preventing the algorithm from becoming overly sensitive to noise or artifacts in the PPG signal.
    % This is important for maintaining the reliability of pulse detection, ensuring that the algorithm does not mistake random signal fluctuations for actual heartbeats.
    if isfield(params,'ppg_learnthresh')
        minthresh = params.ppg_learnthresh;
    else
        minthresh = 5;  % default = 5
    end

    % The "Eye-Closing Period" refers to a specific time duration immediately following the detection of a pulse wave,
    % during which the algorithm refrains from detecting another pulse.
    % This is essentially a kind of refractory period specific to the PPG signal processing,
    % ensuring that the algorithm doesn't incorrectly identify multiple peaks (or pulse waves)
    % too close to each other, which could be artifacts or noise rather than actual physiological signals.
    if isfield(params,'ppg_eyeclosing')
        EyeClosing =  round(fs*params.ppg_eyeclosing);
    else
        EyeClosing = round(fs*0.65);     % default = 0.65 s (range: .4-.8)
    end

    % ExpectPeriod represents the maximum expected duration between two consecutive pulse waves (or heartbeats) in the PPG signal.
    % By dynamically adjusting the detection threshold based on the ExpectPeriod,
    % the algorithm can adapt to variations in pulse signal strength and quality, increasing the likelihood of detecting valid pulses even when the signal is weak or noisy.
    if isfield(params,'ppg_expctperiod')
        ExpectPeriod = round(fs*params.ppg_expctperiod);
    else
        ExpectPeriod = round(fs*5);    % default = 5 s
    end

    % SLPwindow defines the size of the window used for calculating the slope of the PPG signal.
    % The size of the SLPwindow affects how sensitive the algorithm is to changes in the signal.
    % A smaller window might make the algorithm more sensitive to rapid changes,
    % whereas a larger window might smooth out short-term fluctuations,
    % potentially improving detection stability but potentially missing rapid changes.
    if isfield(params,'ppg_slopewindow')
        SLPwindow = round(fs*params.ppg_slopewindow);
    else
        SLPwindow = round(fs*0.1);   % default = 0.1 s (range: .05-.4).
    end

    % INVALID signal (constant)
    INVALID_signal = -32758;

    % initiate variables
    timer = 0;
    peaks = [];
    beat_n = 1;
    from = 1;
    to = length(signal);

    % check signal
    if signal(1) <= INVALID_signal
        signal(1) = mean(signal);
    end
    inv = find(signal<=INVALID_signal);
    for i = 1:length(inv)
        signal(inv(i)) = signal(inv(i)-1);
    end

    % re-scale signal to +/-2000
    if length(signal) < 5*60*fs
        signal = (signal-min(signal))./(max(signal)-min(signal)).*4000-2000;
    else
        % find max/min every 5 minute for re-scaling signal
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
    T0 = T0/n; % T0=T0/(t1-from);
    Ta = 3*T0;

    learning = 1;  % turn learning mode ON

    % Main loop
    t = from;
    while t <= to

        if learning
            if t > from + LPERIOD  % end of learning period
                learning = 0;  % turn learning mode OFF
                T1 = T0;
                t = from;	% start over
            else
                T1 = 2*T0;
            end
        end

        [temp,ebuf,lbuf,tt_2] = slpsamp(t,signal,BUFLN,ebuf,lbuf,tt_2,SLPwindow);

        if temp > T1    % possible pulse near t
            timer = 0;
            % for counting the time after previous pulse
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

                % find valley from the original signal around 0.25 s of tpq
                valley_v = round(tpq);
                for valley_i = round(max(2,tpq-round(0.20*fs))):round(min(tpq+round(0.05*fs),length(signal)-1))

                    % If vally is too low, it cannot serve as an index, so move to the next time.
                    if valley_v <= 0
                        t = t + 1;
                        continue;
                    end

                    if signal(valley_v)>signal(valley_i) && signal(valley_i)<=signal(valley_i-1) && signal(valley_i)<=signal(valley_i+1)
                        valley_v = valley_i;
                    end
                end


                if ~learning

                    % If we are looking for the first peak
                    if beat_n == 1

                        % If the proposed peak index > 0
                        if round(valley_v) > 0
                            peaks(beat_n) = round(valley_v);
                            beat_n = beat_n + 1;
                        end
                    else
                        % Check if rounded valley_v is greater than the prior beat index
                        if round(valley_v) > peaks(beat_n-1)
                            peaks(beat_n) = round(valley_v);
                            beat_n = beat_n + 1;
                        end
                    end
                end


                % Adjust thresholds
                Ta = Ta + (maxd - Ta)/10;
                T1 = Ta / 3;

                % Lock out further detections during the eye-closing period
                t = tpq+EyeClosing;
            end
        else
            if ~learning
                % After learning period, decrease threshold if no pulse was detected recently
                timer = timer+1;
                if timer > ExpectPeriod && Ta > minthresh
                    Ta = Ta-1;
                    T1 = Ta / 3;
                end
            end
        end
        t=t+1;
    end

    sig = signal; % for plotting
    RR = diff(peaks) ./ fs;
    RR_t = peaks ./ fs;
    HR = 60 ./ diff(tm(peaks));   % heart rate (in bpm)

else
    error('Signal type must be ECG or PPG')
end

%% Subfunction

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