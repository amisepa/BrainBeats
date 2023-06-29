%% Process PPG signal
% Detect heart beats (pulse waveforms) in a continuous PPG signal. Works
% best with adult human signals sampled at 125 Hz, but works with any sampling
% rate using on-the-fly resampling. The output marks the heart beats onsets
%
% They are remove/interpolate artifacts
%
% Original code from the Physionet Cardiovascular Signal Processing Toolbox
% function called qppg
% Original authors: W. Zong (1998), revised by G. Moody (2010), Qiao Li
% (2010-2011), Adriana Vest (2011), and Giulia Da Poian (2018).
%
% When using this code, please cite:
% Vest et al. (2018). An Open Source Benchmarked Toolbox for Cardiovascular
%   Waveform and Interval Analysis. Physiological measurement.
%
% Copyright (C) - Cedric Cannard, 2023, BrainBeats toolbox

function [heartbeats, signal] = process_ppg_fast(signal,fs)

% start and end of time series
from = 1;
to = length(signal);

if ~exist('fs','var')
    error("You must provide your signal' sampling rate as 2nd input")
end

% params
BUFLN = 4096;                       % must be a power of 2
LPERIOD  = fs*8;                    % learning period in samples (default = 8 s)
Tm = 5;                             % minimum threshold value (default = 5)
EyeClosing = round(fs * .34);       % eye-closing period is set to 0.34 sec (340 ms)
ExpectPeriod = round(fs * 2.5);	    % threshold in s (default = 2.5) -> adjust if no pulse found
SLPwindow = round(fs * .17);        % Slope window size (deafult = 170 ms)
timer = 0;

heartbeats = [];
beat_n = 1;

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
    n=1;
    for i=1:5*60*fs:length(signal)
        max_signal(n)=max(signal(i:min(i+5*60*fs-1,length(signal))));
        min_signal(n)=min(signal(i:min(i+5*60*fs-1,length(signal))));
        n=n+1;
    end
    signal = (signal-median(min_signal))./(median(max_signal)-median(min_signal)).*4000-2000;
end


ebuf(1:BUFLN)=0;
lbuf=ebuf;
if from>BUFLN
    tt_2=from-BUFLN;
else
    tt_2=0;
end
aet=0;

t1=8*fs;
t1 = t1+from;
T0 = 0;
n=0;
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

%    /* Main loop */
t = from;
while t <= to

    if (learning)
        if (t > from + LPERIOD)
            learning = 0;
            T1 = T0;
            t = from;	% /* start over */
        else
            T1 = 2*T0;
        end
    end

    [temp,ebuf,lbuf,tt_2, aet] = slpsamp(t,signal,BUFLN,ebuf,lbuf,tt_2, aet,SLPwindow);

    if (temp > T1)    % /* found a possible ABP pulse near t */
        timer = 0;
        % /* used for counting the time after previous ABP pulse */
        maxd = temp;
        mind = maxd;
        tmax=t;
        for (tt = t + 1: t + EyeClosing-1)
            [temp2 ,ebuf,lbuf,tt_2, aet] = slpsamp(tt,signal,BUFLN,ebuf,lbuf,tt_2, aet,SLPwindow);
            if temp2 > maxd
                maxd=temp2;
                tmax=tt;
            end
        end
        if (maxd == temp)
            t=t+1;
            continue;
        end

        for tt = tmax :-1: (t - EyeClosing / 2 +1)
            [temp2 ,ebuf,lbuf,tt_2, aet] = slpsamp(tt,signal,BUFLN,ebuf,lbuf,tt_2, aet,SLPwindow);
            if temp2< mind
                mind=temp2;
            end
        end
        if maxd>mind+10
            onset=(maxd-mind)/100+2;
            tpq=t-round(0.04*fs);
            maxmin_2_3_threshold=(maxd-mind)*2.0/3;
            for tt=tmax:-1:t-EyeClosing/2+1
                [temp2, ebuf,lbuf,tt_2, aet] = slpsamp(tt,signal,BUFLN,ebuf,lbuf,tt_2, aet,SLPwindow);
                if temp2<maxmin_2_3_threshold
                    break;
                end
            end
            for tt=tt:-1:t - EyeClosing / 2 + round(0.024*fs)
                [temp2 ,ebuf,lbuf,tt_2, aet] = slpsamp(tt,signal,BUFLN,ebuf,lbuf,tt_2, aet,SLPwindow);
                [temp3 ,ebuf,lbuf,tt_2, aet] = slpsamp(tt-round(0.024*fs),signal,BUFLN,ebuf,lbuf,tt_2, aet,SLPwindow);
                if temp2-temp3<onset
                    tpq=tt-round(0.016*fs);
                    break;
                end
            end

            % find valley from the original signal around 0.25s of tpq
            valley_v = round(tpq);
            for valley_i=round(max(2,tpq-round(0.20*fs))):round(min(tpq+round(0.05*fs),length(signal)-1))

                % If vally is too low, it cannot serve as an index, so move to the next time.
                if valley_v <= 0
                    t = t + 1;
                    continue;
                end

                if signal(valley_v)>signal(valley_i) && signal(valley_i)<=signal(valley_i-1) && signal(valley_i)<=signal(valley_i+1)
                    valley_v=valley_i;
                end
            end


            if (~learning)

                % If we are looking for the first peak
                if beat_n == 1

                    % If the proposed peak index > 0
                    if round(valley_v) > 0
                        heartbeats(beat_n) = round(valley_v);
                        beat_n = beat_n + 1;
                    end
                else
                    % Check if rounded valley_v is greater than the prior beat index
                    if round(valley_v) > heartbeats(beat_n-1)
                        heartbeats(beat_n) = round(valley_v);
                        beat_n = beat_n + 1;
                    end
                end
            end


            % /* Adjust thresholds */
            Ta = Ta + (maxd - Ta)/10;
            T1 = Ta / 3;

            % /* Lock out further detections during the eye-closing period */
            t = tpq+EyeClosing;
        end
    else
        if (~learning)
    	    % /* Once past the learning period, decrease threshold if no pulse
    	    %   was detected recently. */
            timer = timer+1;
            if (timer > ExpectPeriod && Ta > Tm)
                Ta=Ta-1;
                T1 = Ta / 3;
            end
        end
    end

    t=t+1;

end

% Discard first beat because algorithm always finds first minimum value, so trace-back logic
% will find a fir
heartbeats(1) = [];

%% Subfunction
function [beat1,ebuf,lbuf,tt_2, aet] = slpsamp(t,signal,BUFLN,ebuf,lbuf,tt_2, aet,SLPwindow)


while (t > tt_2)
    prevVal=0;

    if (tt_2>0) && (tt_2-1>0) && (tt_2<length(signal)) && (tt_2-1<length(signal))
        val2=signal(tt_2 - 1);
        val1=signal(tt_2);
    else
        val2=prevVal;
        val1=val2;
    end
    prevVal=val2;
    dy =  val1-val2;
    if (dy < 0)
        dy = 0;
    end
    tt_2=tt_2+1;
    M=round(mod(tt_2,(BUFLN-1))+1);
    et=dy;
    ebuf(M)=et;
    %         M2=round(mod(tt_2-SLPwindow,(BUFLN-1))+1);
    %         aet=aet+et-ebuf(M2);
    aet=0;
    for i=0:SLPwindow-1
        p=M-i;
        if p<=0
            p=p+BUFLN;
        end
        aet=aet+ebuf(p);
    end
    lbuf(M) = aet;

end
M3=round(mod(t,(BUFLN-1))+1);
beat1=lbuf(M3);

