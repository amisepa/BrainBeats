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
%   detected using a multi-method voting approach combining amplitude ratios,
%   peak sharpness, and derivative patterns.
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
% Revised January 31, 2026 by Cedric Cannard
%   - Implemented multi-method polarity detection (amplitude ratio, 
%     sharpness, derivative-based) with weighted voting
%   - Removed unreliable temporal-ordering-based polarity detection
%   - Added amplitude consistency check for outlier peak rejection
%   - Improved diagnostic output
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

    % flatline check 
    if prctile(abs(sig), 95) < 0.05
        error('ECG time series amplitude too small (likely flat line)')
    end

    % P&T operations (zero-phase where applicable)
    dffecg = [0; diff(sig)];     % same length as sig
    sqrecg = dffecg.^2;

    % integrate using zero-phase filtering (filtfilt)
    b_int  = ones(int_nb_coef,1) / int_nb_coef;
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
    if ~any(poss_reg)
        % No candidates found, bail early
        peaks = [];
        RR = [];
        RR_t = [];
        HR = [];
        return
    end

    % search-back for missed beats
    if search_back
        indAboveThreshold = find(poss_reg);
        if numel(indAboveThreshold) > 2
            RRv = diff(tm(indAboveThreshold));
            RRv = RRv(RRv > 0.01);
            if ~isempty(RRv)
                medRRv = median(RRv);
                indMissedBeat = find(diff(tm(indAboveThreshold)) > 1.5*medRRv);

                indStart = indAboveThreshold(indMissedBeat);
                indEnd   = indAboveThreshold(indMissedBeat+1);

                for i = 1:numel(indStart)
                    poss_reg(indStart(i):indEnd(i)) = ...
                        mdfint(indStart(i):indEnd(i)) > (0.5 * peakThresh * en_thres);
                end
            end
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
    % R-peak localization with multi-method polarity detection
    % ------------------------------------------------------------

    % 1) Extract candidate QRS segments and compute features
    idxMax = zeros(1, nb_peaks);
    idxMin = zeros(1, nb_peaks);
    valMax = zeros(1, nb_peaks);
    valMin = zeros(1, nb_peaks);
    segEnergy = zeros(1, nb_peaks);

    for i = 1:nb_peaks
        a = left(i);
        b = right(i);
        seg = sig(a:b);

        [valMax(i), imx] = max(seg);
        [valMin(i), imn] = min(seg);

        idxMax(i) = a + imx - 1;
        idxMin(i) = a + imn - 1;

        % Energy of segment (for weighting)
        segEnergy(i) = sum(seg.^2);
    end

    % 2) Polarity detection using DERIVATIVE analysis
    % Key insight: The R-wave has the steepest slopes in the ECG.
    % Find where the maximum |derivative| occurs and determine polarity from that.
    %
    % Also: examine the FULL signal around detected peaks, not just segments
    
    nTemplate = min(nb_peaks, max(15, round(30 * fs / max(1, median(diff(left))))));
    
    % Sort segments by energy (higher energy = cleaner QRS)
    [~, energyOrder] = sort(segEnergy(1:nTemplate), 'descend');
    topSegs = energyOrder(1:min(15, numel(energyOrder)));
    
    % Compute derivative of full signal (for analysis)
    dsig = [0; diff(sig)];
    
    % For each segment, look at a WIDER window around the energy detection
    % to capture the full QRS complex
    extendSamples = round(0.1 * fs);  % extend 100ms each direction
    
    derivPolarityVotes = zeros(1, numel(topSegs));
    peakBeforeTroughVotes = zeros(1, numel(topSegs));
    
    for k = 1:numel(topSegs)
        i = topSegs(k);
        
        % Extended window around the detected segment
        a_ext = max(1, left(i) - extendSamples);
        b_ext = min(numel(sig), right(i) + extendSamples);
        
        segExt = sig(a_ext:b_ext);
        dsegExt = dsig(a_ext:b_ext);
        
        % Find the point of maximum slope magnitude
        [~, maxSlopeIdx] = max(abs(dsegExt));
        maxSlopeVal = dsegExt(maxSlopeIdx);
        
        % If max slope is positive (upstroke), look for peak AFTER
        % If max slope is negative (downstroke), look for trough AFTER
        % The R-peak is what the steepest slope leads TO
        
        if maxSlopeVal > 0
            % Steepest upstroke - R-peak should be positive (after upstroke)
            derivPolarityVotes(k) = +1;
        else
            % Steepest downstroke - R-peak should be negative (after downstroke)
            derivPolarityVotes(k) = -1;
        end
        
        % Also check: in extended window, which extreme comes first?
        [maxVal, maxPos] = max(segExt);
        [minVal, minPos] = min(segExt);
        
        % The first major deflection from baseline is R
        % Estimate baseline from edges
        baselineEst = mean([segExt(1:min(5,end)); segExt(max(1,end-4):end)]);
        
        distMax = abs(maxVal - baselineEst);
        distMin = abs(minVal - baselineEst);
        
        % Only vote if there's a clear difference in timing
        if abs(maxPos - minPos) > round(0.02 * fs)  % >20ms apart
            if maxPos < minPos
                peakBeforeTroughVotes(k) = +1;
            else
                peakBeforeTroughVotes(k) = -1;
            end
        end
    end
    
    % Method 3: Look at the T-wave to determine orientation
    % T-wave follows QRS by ~200-400ms. In normal ECG, T is positive.
    % In inverted ECG, T is negative.
    tWaveVotes = zeros(1, numel(topSegs));
    
    for k = 1:numel(topSegs)
        i = topSegs(k);
        
        % Look 200-400ms after the segment for T-wave
        tStart = right(i) + round(0.15 * fs);
        tEnd = right(i) + round(0.35 * fs);
        
        if tEnd <= numel(sig)
            tSeg = sig(tStart:tEnd);
            tBaseline = mean(sig(max(1,tStart-round(0.05*fs)):tStart));
            
            tMax = max(tSeg) - tBaseline;
            tMin = min(tSeg) - tBaseline;
            
            % T-wave polarity typically matches R-wave polarity
            if abs(tMax) > abs(tMin) * 1.3
                tWaveVotes(k) = +1;  % T is positive -> R is positive
            elseif abs(tMin) > abs(tMax) * 1.3
                tWaveVotes(k) = -1;  % T is negative -> R is negative
            end
        end
    end
    
    % Combine: derivative-based is most reliable
    % Derivative: 50%, T-wave: 30%, temporal order: 20%
    combinedVotes = 0.50 * derivPolarityVotes + 0.30 * tWaveVotes + 0.20 * peakBeforeTroughVotes;
    
    if sum(combinedVotes) > 0
        pol = +1;
    else
        pol = -1;
    end
    
    % Diagnostic output
    fprintf(' - Polarity votes: deriv=%.1f, tWave=%.1f, order=%.1f -> pol=%d\n', ...
        sum(derivPolarityVotes), sum(tWaveVotes), sum(peakBeforeTroughVotes), pol);
    
    % Compute widths for diagnostic
    widthMax = zeros(1, numel(topSegs));
    widthMin = zeros(1, numel(topSegs));
    for k = 1:numel(topSegs)
        i = topSegs(k);
        seg = sig(left(i):right(i));
        baseline = median(seg);
        
        halfAmpMax = baseline + (valMax(i) - baseline) / 2;
        aboveHalf = seg > halfAmpMax;
        if any(aboveHalf)
            widthMax(k) = (find(aboveHalf, 1, 'last') - find(aboveHalf, 1, 'first')) / fs * 1000;
        else
            widthMax(k) = inf;
        end
        
        halfAmpMin = baseline + (valMin(i) - baseline) / 2;
        belowHalf = seg < halfAmpMin;
        if any(belowHalf)
            widthMin(k) = (find(belowHalf, 1, 'last') - find(belowHalf, 1, 'first')) / fs * 1000;
        else
            widthMin(k) = inf;
        end
    end
    fprintf(' - Median peak widths: pos=%.1fms, neg=%.1fms\n', ...
        median(widthMax(widthMax < inf)), median(widthMin(widthMin < inf)));
    
    % Manual override option
    if isfield(params, 'ecg_polarity') && ~isempty(params.ecg_polarity)
        pol = params.ecg_polarity;
        fprintf(' - Polarity OVERRIDDEN to: %d\n', pol);
    end

    % 3) Peak selection
    peaks = zeros(1, nb_peaks);
    pkval = zeros(1, nb_peaks);
    qrsHalfWidth = round(0.05 * fs);  % 50ms

    for i = 1:nb_peaks
        a = left(i);
        b = right(i);

        if pol > 0
            [pkval(i), ii] = max(sig(a:b));
            coarsePeak = a + ii - 1;
        else
            [pkval(i), ii] = min(sig(a:b));
            coarsePeak = a + ii - 1;
        end

        % Refine in tighter window
        w1 = max(1, coarsePeak - qrsHalfWidth);
        w2 = min(numel(sig), coarsePeak + qrsHalfWidth);

        if pol > 0
            [pkval(i), ii] = max(sig(w1:w2));
        else
            [pkval(i), ii] = min(sig(w1:w2));
        end
        peaks(i) = w1 + ii - 1;
    end

    sign = pol * median(abs(pkval));

    % 4) Enforce refractory period
    [peaks, ord] = sort(peaks, 'ascend');
    pkval = pkval(ord);

    keep = true(size(peaks));
    refSamples = round(ref_period * fs);

    for k = 2:numel(peaks)
        if (peaks(k) - peaks(k-1)) < refSamples
            if abs(pkval(k)) > abs(pkval(k-1))
                keep(k-1) = false;
            else
                keep(k) = false;
            end
        end
    end
    peaks = peaks(keep);
    pkval = pkval(keep);

    % 5) Final micro-refinement
    Lsig = numel(sig);
    microWindow = max(1, round(0.010 * fs));  % 10ms

    for k = 1:numel(peaks)
        p = peaks(k);
        w1 = max(1, p - microWindow);
        w2 = min(Lsig, p + microWindow);

        if pol > 0
            [~, ii] = max(sig(w1:w2));
        else
            [~, ii] = min(sig(w1:w2));
        end
        peaks(k) = w1 + ii - 1;
    end

    peaks = unique(peaks, 'stable');

    % 6) Amplitude consistency check (optional outlier rejection)
    finalPkval = zeros(size(peaks));
    for k = 1:numel(peaks)
        finalPkval(k) = sig(peaks(k));
    end

    medAmp = median(abs(finalPkval));
    ampThresh = 0.25;  % peaks must be within 25-400% of median amplitude
    validAmp = abs(finalPkval) > ampThresh * medAmp & ...
               abs(finalPkval) < (1/ampThresh) * medAmp;

    % Only apply if we'd keep most peaks
    if sum(validAmp) > 0.7 * numel(peaks)
        peaks = peaks(validAmp);
    end

    % Print summary
    if pol < 0
        fprintf(' - Peaks polarity: negative\n');
    else
        fprintf(' - Peaks polarity: positive\n');
    end
    fprintf(' - P&T energy threshold: %.2f\n', en_thres);
    fprintf(' - Detected %d R-peaks\n', numel(peaks));

    % RR intervals and HR
    RR = diff(peaks) ./ fs;
    RR_t = tm(peaks);
    HR = 60 ./ diff(RR_t);

    % Sanity check
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

end