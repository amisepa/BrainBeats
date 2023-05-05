% Detect QRS onsets and J-points using a nonlinearly-scaled ECG curve
% length feature. This version has been optimized for adult human ECG
% sampled at 125 Hz, but it can analyze ECGs sampled at any frequency
% with on-the-fly resampling.
%
% INPUTS:
%   signal:   ECG signal, the swqrsm.c used the signal in raw adus, if the input
%           signal is in physical units, when the function found the peak-peak
%           value < 10, it will multiply the value by WFDB_DEFGAIN (200)
%   t:      time index for resampling
%   srate:     sampling frequency (default: 125)
%   lineFreq: power line frequency
%   TmDEF:  minimum threshold value (default: 100)
%   jflag:  annotate J-points (ends of QRS complexes) (default: 0)
%
% OUTPUT:
%   qrs: QRS fiducial mark in samples
%   jpoints: J-points annotation, if jflag==1
%
% Code adapted from the Physionet Cardiovascular toolbox.
% 
% Cedric Cannard, 2022

function [r, jpoints] = get_rwave2(signal,srate,lineFreq,TmDEF,jflag)

if nargin < 5
    jflag = 0;
end
if nargin < 4
    TmDEF = 100; % minimum threshold value (100 = default)
end

% Warning if srate is no 125 hz
if srate ~= 125
    warning('Sample rate is not 125 Hz. This is strongly recommended to get an accurate Signal quality index (SQI)!')
%     [N,D] = rat(125/srate);                % rational fraction approximation
% %     tmp = [125/srate, N/D]                 % approximation accuracy check
%     newData = resample(signal, N, D);     % resampled signal
%     newTimes = resample(t, N, D);     % Resampled times (for plotting only)
%     figure; plot(t(1:srate),signal(1:srate)); hold on; plot(newTimes(1:125),newData(1:125)); 
%     legend(['Previous sample rate: ' num2str(srate) ' Hz'], 'Resampled to 125 Hz')
%     set(gca,'XTick',[]); xlabel('1 second of signal')
%     signal = newData;
%     oldFs = srate;
%     srate = 125;
end
swqrsm.bufln = 16384;   % must be a power of 2, see ltsamp()
baseline = 0.25;         % eye-closing period is set to 0.25 sec (250 ms)
MaxQRSw = 0.13;         % maximum QRS width (130ms)
NDP	= 2.5;              % adjust threshold if no QRS found in NDP seconds
timer_d = 0;
gain = 200.0;
swqrsm.lfsc = fix(1.25*gain*gain/srate);	% /* length function scale constant */
Tm = fix(TmDEF / 5.0);

% test signal is physical units (mV) or raw adus units
datatest = signal(1:fix(length(signal)/srate)*srate);
if length(datatest)>srate
    datatest = reshape(datatest,srate,[]);
end
test_ap = median(max(datatest)-min(datatest));
if test_ap < 10 % peak-peak < 10 mV, may be physical units
    signal = signal*gain;
end

swqrsm.signal = signal;
swqrsm.lbuf = zeros(swqrsm.bufln,1);
swqrsm.ebuf = zeros(swqrsm.bufln,1);
swqrsm.ebuf(1:end)=fix(sqrt(swqrsm.lfsc));
swqrsm.lt_tt = 0;
swqrsm.aet = 0;
swqrsm.Yn = 0;
swqrsm.Yn1 = 0;
swqrsm.Yn2 = 0;

r = [];
jpoints = [];


% The LP filter will have a notch at the power line frequency
swqrsm.LPn = fix(srate/lineFreq);
if swqrsm.LPn > 8
    swqrsm.LPn = 8;	% avoid filtering too agressively
end
swqrsm.LP2n = 2 * swqrsm.LPn;
EyeClosing = fix(srate * baseline); % set eye-closing period
ExpectPeriod = fix(srate * NDP);   % maximum expected RR interval
swqrsm.LTwindow = fix(srate * MaxQRSw);   % length transform window size

% Average the first 8 seconds of the length-transformed samples
% to determine the initial thresholds Ta and T0. The number of samples
% in the average is limited to half of the ltsamp buffer if the sampling
% frequency exceeds about 2 KHz.
% for i = 1:2000
%   ltdata(i) = ltsamp(i);
% end

t1 = srate*8;
if t1 > fix(swqrsm.bufln*0.9)
    t1 = swqrsm.bufln/2;
end

T0 = 0;
for t = 1:t1
    swqrsm = ltsamp(t,swqrsm);
    T0 = T0 + swqrsm.lt_data;
end

T0 = T0 / t1;
Ta = 3 * T0;
t = 1;
learning = 1;
while t < length(signal)
    if learning == 1
        if t > t1
            learning = 0;
            T1 = T0;
            t = 1;	% /* start over */
        else
            T1 = 2*T0;
        end
    end

    % Compare a length-transformed sample against T1.
    swqrsm = ltsamp(t,swqrsm);
    if swqrsm.lt_data > T1   % found a possible QRS near t
        timer_d = 0;            % used for counting the time after previous QRS
        maxd = swqrsm.lt_data;
        mind = maxd;
        for tt = t+1:t + fix(EyeClosing/2)
            swqrsm = ltsamp(tt,swqrsm);
            if (swqrsm.lt_data > maxd)
                maxd = swqrsm.lt_data;
            end
        end
        for tt = t-1:-1:t - fix(EyeClosing/2)
            swqrsm = ltsamp(tt,swqrsm);
            if (swqrsm.lt_data < mind)
                mind = swqrsm.lt_data;
            end
        end

        % if there is a QRS near tt, find QRS onset (PQ junction)
        if maxd > mind+10
            onset = fix(maxd/100) + 2;
            tpq = t - 5;
            for tt = t:-1:t-fix(EyeClosing/2)
                swqrsm = ltsamp(tt,swqrsm);
                swqrsm_1 = ltsamp(tt-1,swqrsm);
                swqrsm_2 = ltsamp(tt-2,swqrsm_1);
                swqrsm_3 = ltsamp(tt-3,swqrsm_2);
                swqrsm_4 = ltsamp(tt-4,swqrsm_3);
                if (swqrsm.lt_data - swqrsm_1.lt_data < onset && swqrsm_1.lt_data - swqrsm_2.lt_data < onset && swqrsm_2.lt_data - swqrsm_3.lt_data < onset && swqrsm_3.lt_data - swqrsm_4.lt_data < onset)
                    swqrsm = swqrsm_4;
                    tpq = tt - swqrsm.LP2n;  % account for phase shift
                    break;
                end
            end

            if learning ~= 1
                if tpq > length(signal), break; end

                % Record an annotation at the QRS onset
                r = [r tpq];

                % J-points processing
                if jflag
                    tj = t+5;
                    for tt = t:t + fix(EyeClosing/2)
                        swqrsm = ltsamp(tt, swqrsm);
                        if swqrsm.lt_data > maxd - fix(maxd/10)
                            tj = tt;
                            break;
                        end
                    end
                    if tj > length(signal), break; end

                    % Record an annotation at the J-point
                    jpoints = [jpoints tj];
                end
            end

            % Adjust thresholds
            Ta = Ta + (maxd - Ta)/10;
            T1 = Ta / 3;

            % Lock out further detections during the eye-closing period
            t = t + EyeClosing;
        end

        % Once past the learning period, decrease threshold if no QRS was detected
    elseif learning ~= 1
        timer_d = timer_d + 1;
        if timer_d > ExpectPeriod && Ta > Tm
            Ta = Ta -1;
            T1 = Ta / 3;
        end
    end
    t = t+1;
end

%-------------- SUBFUNCTION --------------%
% ltsamp() returns a sample of the length transform of the input at time t.
% Since this program analyzes only one signal, ltsamp() does not have an
% input argument for specifying a signal number; rather, it always filters
% and returns samples from the signal designated by the global variable 'sig'.
% The caller must never "rewind" by more than swqrsm.bufln samples (the
% length of ltsamp()'s buffers).
function swqrsm = ltsamp(t,swqrsm)

while (t > swqrsm.lt_tt)
    swqrsm.Yn2 = swqrsm.Yn1;
    swqrsm.Yn1 = swqrsm.Yn;
    v0=swqrsm.signal(1);
    v1=swqrsm.signal(1);
    v2=swqrsm.signal(1);
    if swqrsm.lt_tt>0 && swqrsm.lt_tt<=length(swqrsm.signal)
        v0 = swqrsm.signal(swqrsm.lt_tt);
    end
    if swqrsm.lt_tt-swqrsm.LPn>0 && (swqrsm.lt_tt-swqrsm.LPn)<=length(swqrsm.signal)
        v1 = swqrsm.signal(swqrsm.lt_tt-swqrsm.LPn);
    end
    if swqrsm.lt_tt-swqrsm.LP2n>0 && (swqrsm.lt_tt-swqrsm.LP2n)<=length(swqrsm.signal)
        v2 = swqrsm.signal(swqrsm.lt_tt-swqrsm.LP2n);
    end
    if v0~=-32768 && v1~=-32768 && v2~=-32768
        swqrsm.Yn = 2*swqrsm.Yn1 - swqrsm.Yn2 + v0 - 2*v1 + v2;
    end
    dy = fix((swqrsm.Yn - swqrsm.Yn1) / swqrsm.LP2n);	%	/* lowpass derivative of input */
    swqrsm.lt_tt=swqrsm.lt_tt+1;
    swqrsm.et = fix(sqrt(swqrsm.lfsc +dy*dy)); % /* length transform */
    id = mod(swqrsm.lt_tt,swqrsm.bufln);
    if id == 0
        id = swqrsm.bufln;
    end
    swqrsm.ebuf(id) = swqrsm.et;
    id2 = mod(swqrsm.lt_tt-swqrsm.LTwindow,swqrsm.bufln);
    if id2 == 0
        id2 = swqrsm.bufln;
    end
    swqrsm.aet = swqrsm.aet + (swqrsm.et - swqrsm.ebuf(id2));
    swqrsm.lbuf(id) = swqrsm.aet;
    % swqrsm.lbuf contains the average of the length-transformed samples over
    % the interval from tt-swqrsm.LTwindow+1 to tt
end

id3 = mod(t,swqrsm.bufln);
if id3 == 0
    id3 = swqrsm.bufln;
end
swqrsm.lt_data = swqrsm.lbuf(id3);


