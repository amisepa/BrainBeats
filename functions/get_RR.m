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
% Original code from qppg from the Physionet Cardiovascular Signal Processing 
% Toolbox. Original authors: W. Zong (1998), revised by G. Moody (2010), Qiao Li
% (2010-2011), Adriana Vest (2011), and Giulia Da Poian (2018).
% 
% INPUTS:
%   signal      - raw ECG or PPG signal
%   fs          - sampling rate 
%   sig_type    - heart signal type: 'ecg' or 'ppg'
%
% OUTPUTS:
%   RR          - RR intervals
%   RR_t           - time vector
%   Rpeaks      - R peaks
%   sig         - ECG signal filtered
%
% Example:
% 
% When using this code, please cite:
%   Vest et al. (2018). An Open Source Benchmarked Toolbox for Cardiovascular
%   Waveform and Interval Analysis. Physiological measurement.
%
% Copyright (C) - Cedric Cannard, 2023, BrainBeats toolbox

function [RR, RR_t, Rpeaks, sig, tm, HR] = get_RR(signal, fs, sig_type)

nSamp = size(signal,1);
tm = 1/fs:1/fs:nSamp/fs;   % tm = 1/fs:1/fs:ceil(nSamp/fs);

%%%%%%%%%%%%%%%%%%%%%% ECG %%%%%%%%%%%%%%%%%%%%%%
if strcmpi(sig_type, 'ecg')

    % params
    peakThresh = .6;
    search_back = true;     % perform search back
    ref_period = 0.25;      % refractory period
    
    
    % Constants
    med_smooth_nb_coef = round(fs/100);
    int_nb_coef = round(7*fs/256);      % length is 7 for fs = 256Hz
    max_force = [];                     % to force the energy threshold value
    
    % Bandpass filter ECG signal. This sombrero hat has shown to give slightly
    % better results than a standard band-pass filter.
    b1 = [-7.757327341237223e-05  -2.357742589814283e-04 -6.689305101192819e-04 -0.001770119249103 ...
        -0.004364327211358 -0.010013251577232 -0.021344241245400 -0.042182820580118 -0.077080889653194...
        -0.129740392318591 -0.200064921294891 -0.280328573340852 -0.352139052257134 -0.386867664739069 ...
        -0.351974030208595 -0.223363323458050 0 0.286427448595213 0.574058766243311 ...
        0.788100265785590 0.867325070584078 0.788100265785590 0.574058766243311 0.286427448595213 0 ...
        -0.223363323458050 -0.351974030208595 -0.386867664739069 -0.352139052257134...
        -0.280328573340852 -0.200064921294891 -0.129740392318591 -0.077080889653194 -0.042182820580118 ...
        -0.021344241245400 -0.010013251577232 -0.004364327211358 -0.001770119249103 -6.689305101192819e-04...
        -2.357742589814283e-04 -7.757327341237223e-05];
    b1 = resample(b1,fs,250);
    sig = filtfilt(b1,1,double(signal))';
    
    % Plot spectra 
    % figure; 
    % subplot(2,1,1)
    % Fs = 250;
    % [pwr, f] = get_psd(sig,Fs*2,'hann',50,[],Fs,[0 60],'psd');
    % plot(f, pwr);
    % subplot(2,1,2)
    % spectrogram(sig,kaiser(256,5),220,512,Fs);
    % view(-45,65); colormap bone; xlim([0 50])
    
    % If median of 20% of the samples are < min_amp, abort (flat line)
    min_amp = 0.1;       
    % nSamp = length(signal);         % number of input samples
    if length(find(abs(sig)<min_amp))/nSamp > 0.20
        error('ECG time series is a flat line')
    end
    
        % P&T operations
        dffecg = diff(sig');  % (4) differentiate (one datum shorter)
        sqrecg = dffecg.*dffecg; % (5) square ecg
        intecg = filter(ones(1,int_nb_coef),1,sqrecg); % (6) integrate
        mdfint = medfilt1(intecg,med_smooth_nb_coef);  % (7) smooth
        delay  = ceil(int_nb_coef/2);
        mdfint = circshift(mdfint,-delay); % remove filter delay for scanning back through ECG
    
        % P&T threshold
        if nSamp/fs>90
            xs = sort(mdfint(fs:fs*90));
        else
            xs = sort(mdfint(fs:end));
        end
    
        if isempty(max_force)
            if nSamp/fs>10
                ind_xs = ceil(98/100*length(xs));
                en_thres = xs(ind_xs); % if more than ten seconds of ecg then 98% CI
            else
                ind_xs = ceil(99/100*length(xs));
                en_thres = xs(ind_xs); % else 99% CI
            end
        else
            en_thres = max_force;
        end
    
        % build an array of segments to look into
        poss_reg = mdfint>(peakThresh*en_thres);
    
        % in case empty because force threshold and crap in the signal
        if isempty(poss_reg)
            poss_reg(10) = 1;
        end
    
        % P&T QRS detection & search back
        if search_back
            indAboveThreshold = find(poss_reg); % ind of samples above threshold
            RRv = diff(tm(indAboveThreshold));  % compute RRv
            medRRv = median(RRv(RRv>0.01));
            indMissedBeat = find(RRv>1.5*medRRv); % missed a peak?
            % find interval onto which a beat might have been missed
            indStart = indAboveThreshold(indMissedBeat);
            indEnd = indAboveThreshold(indMissedBeat+1);
    
            % look for a peak on this interval by lowering the energy threshold
            for i = 1:length(indStart)
                poss_reg(indStart(i):indEnd(i)) = mdfint(indStart(i):indEnd(i))>(0.5*peakThresh*en_thres);
            end
        end
    
        % find indices into boudaries of each segment
        left  = find(diff([0 poss_reg'])==1);  % remember to zero pad at start
        right = find(diff([poss_reg' 0])==-1); % remember to zero pad at end
    
        % looking for max/min?
        nb_s = length(left<30*fs);
        loc  = zeros(1,nb_s);
        for j=1:nb_s
            [~,loc(j)] = max(abs(sig(left(j):right(j))));
            loc(j) = loc(j)-1+left(j);
        end
        sign = median(sig(loc));  % median was mean & sig was signal originally here
    
        % loop through all possibilities
        compt = 1;
        nb_peaks = length(left);
        maxval = zeros(1,nb_peaks);
        Rpeaks = zeros(1,nb_peaks);
        for i = 1:nb_peaks
            if sign > 0 % if sign is positive then look for positive peaks
                [maxval(compt), Rpeaks(compt)] = max(sig(left(i):right(i))); % sig was signal originally here
            else % if sign is negative then look for negative peaks
                [maxval(compt), Rpeaks(compt)] = min(sig(left(i):right(i))); % sig was signal originally here
            end
    
            % add offset of present location
            Rpeaks(compt) = Rpeaks(compt)-1+left(i);
    
            % refractory period - improve results
            if compt > 1
                if Rpeaks(compt)-Rpeaks(compt-1)<fs*ref_period && abs(maxval(compt))<abs(maxval(compt-1))
                    Rpeaks(compt)=[];
                    maxval(compt)=[];
                elseif Rpeaks(compt)-Rpeaks(compt-1)<fs*ref_period && abs(maxval(compt))>=abs(maxval(compt-1))
                    Rpeaks(compt-1)=[];
                    maxval(compt-1)=[];
                else
                    compt=compt+1;
                end
            else % if first peak then increment
                compt=compt+1;
            end
        end
    
        % r_pos = Rpeaks;       % QRS datapoint positions
        % r_t = tm(Rpeaks);     % QRS timestamp positions
        % r_amp = maxval;       % amplitude at QRS positions
        % hr = 60./diff(r_t);   % heart rate
        
    if sign < 0
        fprintf(" - Peaks' polarity: negative \n");
    else
        fprintf(" - Peaks' polarity: positive \n");
    end
    fprintf(' - P&T energy threshold: %g \n', en_thres)

    % RR intervals and HR
    RR = diff(Rpeaks) ./ fs;   % RR intervals in s
    RR_t = Rpeaks ./ fs;       % RR timestamps (always ignore 1st heartbeat)
    % RR_t = cumsum(Rpeaks);        % alternative method
    HR = 60 ./ diff(tm(Rpeaks));   % heart rate (in bpm)

    % Visualize
    % if params.vis
    %     figure('color','w');
    %     subplot(2,1,1);
    %     plot(tm,sig,'color','#0072BD'); hold on;
    %     plot(r_t,sig(r_pos),'.','MarkerSize',10,'color','#D95319');
    %     % plot(r_t,r_amp,'.','MarkerSize',10,'color','#D95319');
    %     title('Filtered ECG signal + R peaks');
    %     ylabel('mV'); xlim([0 tm(end)]); set(gca,'XTick',[])
    %
    %     subplot(2,1,2);
    %     plot(r_t(1:length(hr)),hr,'--','color','#A2142F','linewidth',1);
    %     xlim([0 tm(end)]);
    %     title('Heart Rate'); xlabel('Time (s)'); ylabel('bpm');
    %
    %     set(findall(gcf,'type','axes'),'fontSize',10,'fontweight','bold');
    % end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PPG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmpi(sig_type, 'ppg')

    
    if ~exist('fs','var')
        error("You must provide your signal' sampling rate as 2nd input")
    end

    BUFLN = 4096;                       % must be a power of 2
    LPERIOD  = fs*8;                    % learning period in samples (default = 8 s)
    thresh = 5;                         % minimum threshold value (default = 5)
    EyeClosing = round(fs* .3);        % eye-closing period is set to 0.3 sec (300 ms)
    ExpectPeriod = round(fs* 2.5);	    % threshold in s (default = 2.5) -> adjust if no pulse found
    SLPwindow = round(fs* .17);        % Slope window size (deafult = 170 ms)
    timer = 0;
    
    Rpeaks = [];  % these are actually onsets of heartbeats (or pulse wave) 
    % but this makes it easier with outputs of the function
    beat_n = 1;

    % start and end of time series
    from = 1;
    to = length(signal);
    
    % check signal
    INVALID_signal = -32758;    
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
    aet = 0;
    
    t1 = 8*fs;
    t1 = t1+from;
    T0 = 0;
    n = 0;
    for t = from:t1
        [temp,ebuf,lbuf,tt_2, aet] = slpsamp(t,signal,BUFLN,ebuf,lbuf,tt_2, aet,SLPwindow);
        if temp > INVALID_signal
            T0 = T0+temp;
            n=n+1;
        end
    end
    T0 = T0/n; % T0=T0/(t1-from);
    Ta = 3 * T0;
    
    learning=1;
    
    % Main loop
    t = from;
    while t <= to
    
        if (learning)
            if (t > from + LPERIOD)
                learning = 0;
                T1 = T0;
                t = from;	% start over
            else
                T1 = 2*T0;
            end
        end
    
        [temp,ebuf,lbuf,tt_2, aet] = slpsamp(t,signal,BUFLN,ebuf,lbuf,tt_2, aet,SLPwindow);
    
        if temp > T1    % possible pulse near t
            timer = 0;
            % for counting the time after previous pulse
            maxd = temp;
            mind = maxd;
            tmax = t;
            for tt = t + 1: t + EyeClosing-1
                [temp2 ,ebuf,lbuf,tt_2, aet] = slpsamp(tt,signal,BUFLN,ebuf,lbuf,tt_2, aet,SLPwindow);
                if temp2 > maxd
                    maxd=temp2;
                    tmax=tt;
                end
            end
            if maxd == temp
                t = t+1;
                continue;
            end
    
            for tt = tmax :-1: t-EyeClosing/2+1
                [temp2 ,ebuf,lbuf,tt_2, aet] = slpsamp(tt,signal,BUFLN,ebuf,lbuf,tt_2, aet,SLPwindow);
                if temp2< mind
                    mind=temp2;
                end
            end
            if maxd > mind+10
                onset = (maxd-mind)/100+2;
                tpq = t-round(0.04*fs);
                maxmin_2_3_threshold=(maxd-mind)*2.0/3;
                for tt = tmax:-1:t-EyeClosing/2+1
                    [temp2, ebuf,lbuf,tt_2, aet] = slpsamp(tt,signal,BUFLN,ebuf,lbuf,tt_2, aet,SLPwindow);
                    if temp2 < maxmin_2_3_threshold
                        break
                    end
                end
                for tt = tt:-1:t - EyeClosing / 2 + round(0.024*fs)
                    [temp2 ,ebuf,lbuf,tt_2, aet] = slpsamp(tt,signal,BUFLN,ebuf,lbuf,tt_2, aet,SLPwindow);
                    [temp3 ,ebuf,lbuf,tt_2, aet] = slpsamp(tt-round(0.024*fs),signal,BUFLN,ebuf,lbuf,tt_2, aet,SLPwindow);
                    if temp2 - temp3<onset
                        tpq = tt-round(0.016*fs);
                        break
                    end
                end
    
                % find valley from the original signal around 0.25s of tpq
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
                            Rpeaks(beat_n) = round(valley_v);
                            beat_n = beat_n + 1;
                        end
                    else
                        % Check if rounded valley_v is greater than the prior beat index
                        if round(valley_v) > Rpeaks(beat_n-1)
                            Rpeaks(beat_n) = round(valley_v);
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
                if timer > ExpectPeriod && Ta > thresh
                    Ta = Ta-1;
                    T1 = Ta / 3;
                end
            end
        end
        t=t+1;
    end
    
    % Discard first beat because algorithm always finds first minimum value, so trace-back logic
    % will find a fir
    % Rpeaks(1) = [];

    RR = diff(Rpeaks) ./ fs;
    RR_t = Rpeaks ./ fs;
    HR = 60 ./ diff(tm(Rpeaks));   % heart rate (in bpm)
    
    sig = signal; % for output
else
    error('Signal type must be ECG or PPG')
end

%% Subfunction

function [beat1,ebuf,lbuf,tt_2, aet] = slpsamp(t,signal,BUFLN,ebuf,lbuf,tt_2, aet,SLPwindow)


while t > tt_2
    prevVal = 0;

    if tt_2>0 && tt_2-1>0 && tt_2<length(signal) && tt_2-1<length(signal)
        val2 = signal(tt_2 - 1);
        val1 = signal(tt_2);
    else
        val2 = prevVal;
        val1 = val2;
    end
    prevVal = val2;
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
