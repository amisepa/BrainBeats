% GET_RWAVE3 - Fast R-peak detector, used as the independent second
% detector for the ECG signal quality index (get_sqi_ecg).
%
% Usage:
%   rwave = get_rwave3(signal, srate)
%
% Inputs:
%   signal - raw ECG (vector)
%   srate  - sampling rate (Hz)
% Outputs:
%   rwave  - R-peak sample indices (column), or -1 if no QRS is found
%
% Method: ~17.5 Hz resonant bandpass, derivative, 30 Hz lowpass, squaring
% and 150 ms integration; marks where the energy exceeds 15% of the local
% (2-s) maximum, at least 350 ms apart. Each R-peak is the largest
% |deflection| of the 1-Hz highpassed ECG within 70 ms of the mark
% (shifted back 30 ms for the causal filter delay). Polarity-independent.
%
% Adapted from heplab_fastdetect (HEPLAB toolbox).
% Adapted for BrainBeats by Cedric Cannard, 2023

function rwave = get_rwave3(signal,srate)

signal = double(signal);
signal = signal(:);

% Resonant bandpass at ~17.5 Hz (QRS energy), then first difference
Q = 3;
gain = 1.2;
w0 = 17.5625*2*pi;
NUM = gain*w0^2;
DEN = [1,(w0/Q),w0^2];
[B,A] = bilinear(NUM,DEN,srate);
ecg_flt = filtfilt(B,A,signal);
ecg_flt = filter([1 -1],1,ecg_flt);

% low-pass 30 Hz (causal), normalize and square
[B,A] = butter(8,30/(srate/2));
ecg_flt = filter(B,A,ecg_flt);
ecg_flt = ecg_flt /max(abs(ecg_flt));
ecg_flt = ecg_flt.^2;

% moving-average integration over 150 ms
N = round(0.150*srate); ecg_flt=1/N*filter(ones(1,N),1,ecg_flt);

% search window for the R-peak after each mark (70 ms)
area = round(0.070*srate);

% threshold: fraction of the local maximum energy
gain = 0.15;

% window for the local maximum (2 s)
comp = round(2*srate);

% minimum interval between marks (350 ms)
step = round(0.350*srate);

% scan step (10 ms)
step_10ms = round(0.01*srate);

% marks are moved back 30 ms before searching for the R-peak
ret = round(0.030*srate);

% scan the energy signal for threshold crossings
sz = length(signal);
n = 1;
Rwave = [];
rwave = [];
while n < sz
    
    if (n+comp) <= sz
        lmt = gain*max(abs(ecg_flt(n:n+comp)));
    else
        lmt = gain*max(abs(ecg_flt(sz-comp:sz)));
    end
    
    if (ecg_flt(n)>lmt && n<sz)
        Rwave = [Rwave ; n];
        n=n+step;
    else
        n=n+step_10ms;
    end
end

mark_count = length(Rwave);

% zero-phase 1-Hz highpass of the raw ECG, to locate the peaks
[B,A] = butter(4,1/(srate/2),'high'); % cf = 1 Hz;
ecg_flt = filtfilt(B,A,signal);

% no QRS found (e.g. flat or disconnected channel)
if isempty(Rwave)
    rwave = -1;
    return
end

% move marks back to compensate for the causal filter delay
Rwave = Rwave-ret;
if Rwave(1)<1
    Rwave(1)=1;
end

% locate R peaks: largest |deflection| within 70 ms after each mark
for i = 1:mark_count
    
    if sz>=Rwave(i) + area
        [~,mark] = max(abs(ecg_flt(Rwave(i):Rwave(i)+area)));
    else
        [~,mark] = max(abs(ecg_flt(Rwave(i):sz)));
    end
    
    % calculate and save mark
    mark = mark+Rwave(i)-1;
    rwave = [rwave ; mark];    
end

if isempty(Rwave)
    rwave = -1;
end
