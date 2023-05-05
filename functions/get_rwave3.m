% From heplab_fastdetect

function rwave = get_rwave3(signal,srate)

signal = double(signal);
signal = signal(:);

% Filter at 17 Hz
Q = 3;
gain = 1.2;
w0 = 17.5625*2*pi;
NUM = gain*w0^2;
DEN = [1,(w0/Q),w0^2];
[B,A] = bilinear(NUM,DEN,srate);
ecg_flt = filtfilt(B,A,signal);
ecg_flt = filter([1 -1],1,ecg_flt);

% low-pass 30 Hz
[B,A] = butter(8,30/(srate/2));
ecg_flt = filter(B,A,ecg_flt);
ecg_flt = ecg_flt /max(abs(ecg_flt));
ecg_flt = ecg_flt.^2;

% integration
N = round(0.150*srate); ecg_flt=1/N*filter(ones(1,N),1,ecg_flt);

% start fast algorithm, area to look for R
area = round(0.070*srate);

% gain
gain = 0.15;

% comparison limit: 2 sec
comp = round(2*srate);

% minimum interval between marks (350ms)
step = round(0.350*srate);

% 10ms stet
step_10ms = round(0.01*srate);

ret = round(0.030*srate);

% repeat for every cardiac cycle
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

% n = 1;
mark_count = length(Rwave);

% filter
[B,A] = butter(4,1/(srate/2),'high'); % cf = 1 Hz;
ecg_flt = filtfilt(B,A,signal);

% return to signal
Rwave = Rwave-ret;
if Rwave(1)<1
    Rwave(1)=1;
end

% locate R peaks
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

