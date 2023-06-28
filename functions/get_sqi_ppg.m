% calculate_ppgsqi.m - run PPG SQI based on beat template correlation on
% 30 s PPG segments in loop
%
% Inputs:
%     beats     - PPG annotation time (samples), read from ple annot file,
%                 But the ann-time is the OFFSET based on wave(1)
%     signal    - PPG signal
%     fs        - sampling frequency (default fs = 125Hz)
%     winlength - window length in s (default = 30)
%
% Outputs:
%     annot     - ppg sqi annotation
%                     E - excellent beat;
%                     A - acceptable beat;
%                     Q - unacceptable beat
%     sqi       - ppg sqi matrix
%                     [N,1]: SQI based on Direct compare
%                     [N,2]: SQI based on Linear resampling
%                     [N,3]: SQI based on Dynamic time warping
%                     [N,4]: SQI based on Clipping detection
%     template  - Current PPG beat template
%     valid     - valid (>=1) or invalid (0) template
%
% Original code by Qiao Li (2011)
%
% Copyright (C) - Cedric Cannard, 2023

function [sqi, sqi_mu, annot] = get_sqi_ppg(beats,signal,fs, winlength)

% Params
winlength = winlength*fs; % default = 30-s window

% Initialize
template = [];
beat_i = 1;
annot = [];
for j = 1:length(beats)
    annot(j) = 'Q';
end

sqimatrix_all = zeros(length(beats),3);
sqi_mu = zeros(1,ceil(length(signal)/fs/30));

% loop every 30-sec
for j = 1:ceil(length(signal)/winlength)
    databegin = (j-1)*winlength+1;
    dataend = min(length(signal),j*winlength);
    annf = find(beats<=dataend);

    if length(annf)<=1, continue; end
    if length(annf) < length(beats)
        annf(length(annf)+1) = annf(length(annf))+1;
    end
    annf = find(beats(annf) >= databegin);
    if length(annf) <= 1, continue;  end
    anntime = beats(annf);

    % prolong the data window an extra 3s
    wave = signal( databegin:min(length(signal), max(dataend, anntime(length(anntime))+3*fs)) ); 

    % set the anntime to be the 0 offset of the selected data
    anntime = beats(annf)-databegin+1;

    % PPG SQI analysis
    [annot, sqimatrix, template, valid] = PPG_SQI_buf(wave,anntime,template,30*fs,fs);
    for k=1:length(annot)
        annot(annf(k))=annot{k};
        sqimatrix_all(annf(k),:)=sqimatrix(k,1:3); % 1:4
        beat_i=beat_i+1;
    end
    sqi_mu(j) = mean(mean(sqimatrix(:,1:3)'));

end
sqi = round(mean(sqimatrix_all(:,1:3),2)');




%% Subfunction 1

% [annot sqimatrix template valid] = PPG_SQI_buf(wave,anntime,template_ahead,windowlen,Fs)
% 
% PPG_SQI.m - PPG SQI based on beat template correlation.
% (as an advice, the algorithm get 30 beats at each call and run in loop)
% by Qiao Li 30 Mar 2011
% 
% input: 
%     wave:       PPG data; 
%     anntime:    PPG annotation time (samples), read from ple annot file,
%                 But the ann-time is the OFFSET based on wave(1)
%     template:   Last PPG beat template 
%     windowlen:  length of window to calculate template(default: 30s)
%     Fs       :  sampling frequency (default Fs=125Hz)
% output:
%     annot:      ppg sqi annotation
%                     E - excellent beat; 
%                     A - acceptable beat; 
%                     Q - unacceptable beat
%     sqimatrix:  ppg sqi matrix   
%                     [N,1]: SQI based on Direct compare
%                     [N,2]: SQI based on Linear resampling
%                     [N,3]: SQI based on Dynamic time warping
%                     [N,4]: SQI based on Clipping detection
%     template:   Current PPG beat template
%     valid:      1 or greater for valid template, 
%                 0 for invalid template
%	
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%
% 03-03-2017
% Edits by Adriana Vest
% - Changed output variable annot from numeric to cell to preserve
%   characters
% - Style changes to align loops and conditional statements
%
% 12-01-2017 Modified by Giulia Da Poian: sampling frequency as input
% parameter instead of fixed fs = 125
% 
% 12-19-2017 Modified by Giulia Da Poian: replaced dp_dtw with dpfast 

function [annot, sqimatrix, template, valid] = PPG_SQI_buf(wave,anntime,template,windowlen,Fs)

    if nargin < 5
        Fs = 125; 
    end
    if nargin < 4 || isempty(windowlen)
        windowlen=30*Fs;
    end
    if nargin < 3 || isempty(template)
        template=[];
    end
    if nargin < 2
        sprintf('Error: must provide wave, anntime');
        annot=[];
        sqimatrix=[];
        template=[];
        valid=0;
        return;
    end

    annot=[];
    sqimatrix=[];
    
    % get PPG template
    [t t2 v]=template_pleth(wave(1:min(windowlen,length(wave))), anntime(find(anntime<min(windowlen,length(wave)))),0, Fs);

    if v<1 && length(template)<1 % Current template invalid && no previous template available
        for j=1:length(anntime)
            annot{j}='Q'; % Unacceptable
        end
        t=[];
    else
        % Using previous template as the template
        if v<1
            t=template;
        end
        % Using t2 if available
        if v>1
            t=t2;
        end
        
        % Calculate the PLA of template for dynamic time warping
        d1=t;
        d1=(d1-min(d1))/(max(d1)-min(d1)).*100;
        [y1 pla1]=PLA(d1,1,1);

        % Main Loop
        for j=1:length(anntime)-1

            % SQI1: Direct compare
            % Calculate correlation coefficients based on the template
            % length
            beatbegin=anntime(j);
            beatend=anntime(j+1);
            % 07/11/2011 ADD max beat length <= 3s detection 
            if beatend-beatbegin>3*Fs
                beatend=beatbegin+3*Fs;
            end
            templatelength=length(t);
            if beatbegin+templatelength-1 > length(wave) || beatend > length(wave) || beatbegin < 1
                continue;
            end
            currentb = j;
            cc = corrcoef(t,wave(beatbegin:beatbegin+templatelength-1));
            c1(j) = cc(1,2);
            if (c1(j)<0)
                c1(j)=0;
            end
            sqimatrix(currentb,1)=int8(c1(j)*100);

            % SQI2: Linear resampling
            % Calculate correlation coefficients based on the 
            % linear resampling (interp1)
            
            y=interp1(1:beatend-beatbegin, wave(beatbegin:beatend-1),1:(beatend-beatbegin-1)/(templatelength-1):(beatend-beatbegin),'spline');
            y(isnan(y))=0;
            cc=corrcoef(t,y);
            c2(j)=cc(1,2);
            if (c2(j)<0)
                c2(j)=0;
            end
            sqimatrix(currentb,2)=int8(c2(j)*100);

            % SQI3: Dynamic Time Warping                
            % Calculate correlation coefficients based on the dynamic time
            % warping
            d2=wave(beatbegin:beatend-1);
            
            % if beat too long, set SQI = 0;
            if (length(d2)>length(d1)*10)
                c3(j)=0;
            else
                d2=(d2-min(d2))/(max(d2)-min(d2)).*100;
               [y2 pla2]=PLA(d2,1,1);

               [w ta tb] = simmx_dtw(y1,pla1,y2,pla2);
               try % try to use the fast version if possible
                   [p,q,Dm] = dpfast(w);
               catch
                   [p,q,Dm] = dp_dtw(w); 
               end
               [ym1, ym2, yout1] = draw_dtw(y1,pla1,p,y2,pla2,q); 
                cc=corrcoef(y1,ym2);
                c3(j)=cc(1,2);
                if (c3(j)<0)
                    c3(j)=0;
                end
            end
            sqimatrix(currentb,3)=int8(c3(j)*100);
                
            % SQI4: Clipping detection   
            d2=wave(beatbegin:beatend-1);
            y=diff(d2);
            clipthreshold=0.5;
            c4(j)=int8(length(find(abs(y)>clipthreshold))/length(y)*100);
            sqimatrix(currentb,4)=c4(j);

            % SQI: Combined

            sqibuf=[sqimatrix(currentb,1) sqimatrix(currentb,2) sqimatrix(currentb,3) sqimatrix(currentb,4)];
            if min(sqibuf)>=90
                annot{currentb}='E'; % Excellent
            else
                if length(find(sqibuf>=90))>=3 || (median(sqibuf(1:3))>=80 && sqibuf(1)>=50 && sqibuf(4)>=70) || min(sqibuf)>=70
                    annot{currentb}='A'; % Acceptable
                else
                    annot{currentb}='Q'; % Unacceptable
                end
            end
                
        end
    end
    
    template=t;
valid=v;



%% subfunction 2

% [template t2 valid] = template_pleth(wave,anntime,temp_ahead,samp_freq)
% 
% template_pleth.m - PPG waveform template creation.
% by Qiao Li 21 Feb 2011
% 
% input: 
%     wave:       PPG data; 
%     anntime:    PPG annotation time (sample), read from ple annot file,
%                 the time is a offset based on wave(1)
%     temp_ahead: N samples before the beginning of PPG waveform mark,
%                 default is 0
%     samp_freq:  sampling frequency, default is 125Hz
%     
% output:
%     template:   PPG waveform template based on normal-length beats
%     t2:         PPG waveform template2 based on the beats which have 
%                 good correlation with the template
%     valid:      1 for valid template, 2 for valid t2, 3 for both valid,
%                 0 for invalid template and t2
%	
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%

function [template, t2, valid] = template_pleth(wave,anntime,temp_ahead,samp_freq)

if nargin < 4
    samp_freq = 125;
end

if nargin < 3
    temp_ahead = 0;
end

% according to heart rate max(300bpm) and min(20bpm) to get max and min
% beat-by-beat interval
hr_max=300;
bb_interval_min=samp_freq*60/hr_max;
hr_min=20;
bb_interval_max=samp_freq*60/hr_min;

normal_beat_length_min=0.7;
normal_beat_lentth_max=1.5;
normal_beat_percent_threshold=0.5;

% using xcorr to get the basic period of the PPG as the length of template
y=xcorr(detrend(wave));

len=length(wave);
lena=length(anntime);
i=len+1;
n=1;
peaki=[];
while i<len*2-1
    if y(i)>y(i-1) && y(i)>y(i+1)
        peaki(n)=y(i);
        peakp(n)=i;
        n=n+1;
    end
    i=i+1;
end
if (length(peaki)<1) || (lena<1)
    template=[];
    t2=[];
    valid=0;
    return;
end
[c i]=max(peaki);
peakp=peakp-len;
i=peakp(i);

cycle=samp_freq;
if i<len-1
    cycle=i;
end

% cumulate the beats with reasonable length to get template
p0=1;
i=anntime(p0);
while i-temp_ahead <1 
    p0=p0+1;
    if (p0>lena)
        template=wave;
        valid=0;
        return;
    end
    i=anntime(p0);
end

if p0+1>=lena
    template=[];
    t2=[];
    valid=0;
    return;
end

beat_interval=diff(anntime(p0:length(anntime)));
median_bi=median(beat_interval);
if ~isnan(median_bi)
    temp_peak=abs(peakp-median_bi);
    [m i]=min(temp_peak);
    cycle=peakp(i);
else
    template=[];
    t2=[];
    valid=0;
    return;
end
   
% the length of template valid detection
valid=1;
if cycle > bb_interval_max || cycle < bb_interval_min
    valid=0;
    template=zeros(1,cycle);
    t2=zeros(1,cycle);
    return;
end
    
n=0;
d1=0;
invalidn=0;
currentbeatlength=anntime(p0+1)-anntime(p0);
if currentbeatlength>0 %cycle*normal_beat_length_min %% && currentbeatlength < cycle*normal_beat_lentth_max
    d1=wave(i-temp_ahead:i+cycle-1);%-temp_ahead);
    n=1;
else
    invalidn=invalidn+1;
    d1=zeros(cycle+temp_ahead,1);
end
    
p0=p0+1;
if p0<lena-1
    i=anntime(p0);
    n=1;
    invalidn=0;
    while i<len-cycle && p0<lena-1
        currentbeatlength=anntime(p0+1)-anntime(p0);
        if currentbeatlength>0 %cycle*normal_beat_length_min %% && currentbeatlength < cycle*normal_beat_lentth_max
            d1=d1+wave(i-temp_ahead:i+cycle-1);%-temp_ahead);
            n=n+1;
        else
            invalidn=invalidn+1;
        end
        p0=p0+1;
        i=anntime(p0);
    end
    d1=d1./n;
    % normal beat is less than the reasonable percentage of all beats
    if (n/(n+invalidn))<normal_beat_percent_threshold
        valid=0;
    end
else
    valid=0;
end

% Compare each beat to the template to get template 2
d2=0;
if (valid)
    p0=2;
    i=anntime(p0);
    n=0;
    while i<len-cycle && p0<lena-1
        cc=corrcoef(d1,wave(i-temp_ahead:i+cycle-1));
        if cc(1,2)>0.8
            d2=d2+wave(i-temp_ahead:i+cycle-1);
            n=n+1;
        end
        p0=p0+1;
        i=anntime(p0);
    end
    d2=d2./n;
    valid = uint8(valid);
    if n>length(anntime)*normal_beat_percent_threshold
        valid = bitor(valid, uint8(2));
    end
end
template=d1;
t2=d2;

%% subfunction 3

% Piecewise linear approximation (PLA) 
% See: HJLM. Vullings, Automated ECG segmentation with Dynamic Time Warping 

function [y1 PLA] = PLA(input,s,th)

% input:    source data
% s:        step
% th:       threshold

if nargin<3
    th=10;
end
if nargin<2
    s=10;
end
    
n=length(input);
% Low pass filter 
% [b,a] = cheby1(3,0.5,40/62.5);
% y1=filter(b,a,input);
% b = fir1(31,40/62.5);
% y=conv(b,input);
% y1=y(int16((length(y)-length(input))/2):(int16((length(y)-length(input))/2+length(input)-1)));

s1=s;
pp=1;
PLA(pp)=1;
pp=pp+1;
i=1;
while i<n
    if (i+s1>=n)
        i_plus_s=n;
    else
        i_plus_s=i+s1;
    end
    interrupt=0;
    while (interrupt==0)
    j=i+1;
    while j<=i_plus_s
        distance=input(i_plus_s)-input(i);
        dcur=input(j)-input(i)-(distance*(j-i)/(i_plus_s-i));
        % distance is larger than threshold, truncate old line
        if abs(dcur)>th
            s1=j-i;
            i_plus_s=i+s1;
            j=i+1;
            interrupt=1;
            continue;
        end
        j=j+1;
    end
    if interrupt==1
        PLA(pp)=j-1;
        pp=pp+1;
        i=j-2;
        s1=s;
    else
        % distance is smaller than threshold, expand old line
        if (i_plus_s>=n)
            i_plus_s=n;
            break;
        else
            if ((i_plus_s+s1)>=n)
                i_plus_s=n;
            else
                i_plus_s=i_plus_s+s1;
            end
        end
    end
    end
    i=i+1;
end
PLA(pp)=n;
y1=input;    
% PLA=y1(int16((length(y1)-length(input))/2):(int16((length(y1)-length(input))/2+length(input)-1)));

%% subfunction 4

% function [w ta tb] = simmx_dtw(y1,pla1,y2,pla2)
%
% calculate a sim matrix between y1 and y2.
%

function [w, ta, tb] = simmx_dtw(y1,pla1,y2,pla2)

% slope1(1)=y1(pla1(1));
% slope2(1)=y2(pla2(1));

slope1(1)=0;
slope2(1)=0;

ta(1)=1;
tb(1)=1;
for i=1:length(pla1)-1
    slope1(i+1)=(y1(pla1(i+1))-y1(pla1(i)))/(pla1(i+1)-pla1(i));
    ta(i+1)=pla1(i+1)-pla1(i);
end
for i=1:length(pla2)-1
    slope2(i+1)=(y2(pla2(i+1))-y2(pla2(i)))/(pla2(i+1)-pla2(i));
    tb(i+1)=pla2(i+1)-pla2(i);
end
al=length(slope1);
bl=length(slope2);
for i=1:bl
    A1(:,i)=slope1;
end
for j=1:al
    B1(j,:)=slope2';
end

w = abs(A1-B1);

%% subfunction 5

function [p,q,D] = dp_dtw(M)
% [p,q] = dp(M) 
%    Use dynamic programming to find a min-cost path through matrix M.
%    Return state sequence in p,q
% 2003-03-15 dpwe@ee.columbia.edu

% Copyright (c) 2003 Dan Ellis <dpwe@ee.columbia.edu>
% released under GPL - see file COPYRIGHT

[r,c] = size(M);

% costs
D = zeros(r+1, c+1);
D(1,:) = NaN;
D(:,1) = NaN;
D(1,1) = 0;
D(2:(r+1), 2:(c+1)) = M;

% traceback
phi = zeros(r,c);

for i = 1:r
  for j = 1:c
    [dmax, tb] = min([D(i, j), D(i, j+1), D(i+1, j)]);
    D(i+1,j+1) = D(i+1,j+1)+dmax;
    phi(i,j) = tb;
  end
end

% Traceback from top left
i = r; 
j = c;
p = i;
q = j;
while i > 1 && j > 1
  tb = phi(i,j);
  if (tb == 1)
    i = i-1;
    j = j-1;
  elseif (tb == 2)
    i = i-1;
  elseif (tb == 3)
    j = j-1;
  else    
    error
  end
  p = [i,p];
  q = [j,q];
end
p=[1,p];
q=[1,q];
% Strip off the edges of the D matrix before returning
D = D(2:(r+1),2:(c+1));

%% subfunction 6

function [ymodify, y2modify, r] = draw_dtw(y1,pla1,p,y2,pla2,q)

l=length(p);
i=1;
j=1;
la=pla1(p(1));
lb=pla2(q(1));
while p(i+1)==p(1)
    i=i+1;
end
while q(j+1)==q(1)
    j=j+1;
end
point=max(i,j)+1;
outputx=[];
outputy=[];
while (point<=l)
    tp=1;
    tq=1;
    while (point+1<=l && ((p(point+1)==p(point)) || (q(point+1)==q(point))))
        point=point+1;
    end
    if point<=l
        laold=la;
        lbold=lb;
        la=pla1(p(point));
        lb=pla2(q(point));
        %intvb=(la-laold)/(lb-lbold);
        intvb=(lb-lbold)/(la-laold)/10;
        xx=pla1(p(i)):intvb:pla1(p(point));
        yy=y2(pla2(q(j)):pla2(q(point)));
        x1=xx(1):(xx(length(xx))-xx(1))/(length(yy)-1):xx(length(xx));
        if length(xx)<=1
            y=yy(end);
        else
            y = griddedInterpolant(x1,yy,'spline');
            y = y(xx);
            % y = interp1(x1,yy,xx,'spline'); % 12-19-2017 Modified by Giulia Da Poian
            % replaced interp1 with griddedInterpolant for speed
        end
%        plot(xx,y,'k');
        outputx=[outputx,xx];
        outputy=[outputy,y];
        i=point;
        j=point;
    end
    point=point+1;
end
% outputxx=unique(outputx);
j=1;
i=1;    
while i<=length(outputx)-1 
    while i<length(outputx)-1 && outputx(i)==outputx(i+1)
        i=i+1;
    end
    outputxx(j)=outputx(i);
    outputyy(j)=outputy(i);
    i=i+1;    
    j=j+1;
end
outputxx(j)=outputx(i);
outputyy(j)=outputy(i);
% 12-01-2017 added by Giulia Da Poian to solve Error: The grid vectors must contain unique points
[~, uidx] = unique(outputxx);
y=interp1(outputxx(uidx),outputyy(uidx),1:(outputxx(length(outputxx(uidx)))-1)/(length(y1)-1):outputxx(length(outputxx(uidx))),'spline');

y(isnan(y))=0;
%plot(y,'m');
y2modify=y;
% ylength=min(y,y1);
[S,R] = size(y);
[R2,Q] = size(y1);
if R ~= R2
    difference=dist(y,y1');
else
    difference=dist(y,y1);
end
    
meany1=sqrt(sum(y1.^2));
r=difference/meany1;
if r<0.0010
    ymodify=y1.*0.9+y'.*0.1;
else
    ymodify=y1;
end