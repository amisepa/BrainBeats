% Find individualized frequency bounds based on power spectra distribution.
% 
% Input power in uv^2 (not normalized) is best!
% 
% Method developed by Corcoran et al. (2017) in the resting IAF toolbox.
% 
% Copyright (C) - Cedric Cannard, 2023, BrainBeats toolbox

function [bounds, peak] = get_freqBounds(pwr, f, fs, w, winSize, mpow)

% Parameters
Fw = 11;            % SGF frame width (11 corresponding to a frequency span of ~2.69 Hz @ ~.24Hz frequency resolution)
k = 5;              % SGF polynomial order (default = 5)
mdiff = .2;    % minimal height difference distinguishing a primary peak from
                % competing peaks (default = 0.2; i.e. 20% peak height)
% mpow = 1;       % error bound (SD) threshold to differentiate peaks from background spectral noise (default = 1)


% fit 1st order poly (regression line) to normalised spectra (log-scaled)
[pfit, sig] = polyfit(f, 10*log10(pwr), 1);     

% derive yval coefficients of fitted polynomial and delta (std dev) error estimate
[yval, del] = polyval(pfit, f, sig);         

% take [minPowThresh * Std dev] as upper error bound on background spectral noise
minPow = yval + (mpow * del);                

% apply Savitzky-Golay filter to fit curves to spectra & estimate 1st and 2nd derivatives
% [d0, d1, d2] = sgfDiff(pwr, Fw, k, fs, winSize);
% [d0, d1, d2] = sgfDiff(x, Fw, poly, Fs, tlen)
[~, g] = sgolay(k, Fw);      
dt = fs/winSize;
dx = zeros(length(pwr),3);
for p = 0:2                 % p determines order of estimated derivatives
    dx(:,p+1) = conv(pwr, factorial(p)/(-dt)^p * g(:,p+1), 'same');
end
d0 = dx(:,1);        % smoothed signal post S-G diff filt
d1 = dx(:,2);        % 1st derivative
d2 = dx(:,3);        % 2nd derivative

% frequency resolution
% exp_tlen = nextpow2(winSize); 
% fres = fs/2.^exp_tlen; 
fres = f(2)-f(1);

% take derivatives, find peak(s) and boundaries of the frequency band
% [peak, pos1, pos2]  = peakBounds(d0, d1, d2, f, w, minPow, mdiff, fres);
% [peakF, posZ1, posZ2] = peakBounds(d0, d1, d2, f, w, minPow, minDiff, fres)

% evaluate derivative for zero-crossings
[~, lower_lim] = min(abs(f-w(1)));      % set lower bound for alpha band
[~, upper_lim] = min(abs(f-w(2)));      % set upper bound for alpha band

negZ = zeros(1,4);  % zero-crossing count & frequency bin
cnt = 0;     
for k = lower_lim-1:upper_lim+1             % step through frequency bins in alpha band (start/end at bound -/+ 1 to make sure don't miss switch)
    if sign(d1(k)) > sign(d1(k+1))              % look for switch from positive to negative derivative values (i.e. downward zero-crossing)
       	[~, maxk] = max([d0(k), d0(k+1)]);      % ensure correct frequency bin is picked out (find larger of two values either side of crossing (in the smoothed signal))
       	if maxk == 1
        	maxim = k;
        elseif maxk == 2
            maxim = k+1;
        end
        cnt = cnt+1;                % advance counter by 1
      	negZ(cnt,1) = cnt;          % zero-crossing (i.e. peak) count
        negZ(cnt,2) = maxim;        % keep bin index for later
      	negZ(cnt,3) = f(maxim);     % zero-crossing frequency            
     	negZ(cnt,4) = d0(maxim);    % power estimate
    end
end
    
% sort out appropriate estimates for output
if negZ(1,1) == 0                   % if no zero-crossing detected --> report NaNs
    peak = NaN;
    subBin = NaN;
elseif size(negZ, 1) == 1           % if singular crossing...
    if log10(negZ(1, 4)) > minPow(negZ(1,2))      % ...and peak power is > minimum threshold --> report frequency
        peakBin = negZ(1, 2);
        peak = negZ(1, 3);
    else
        peak = NaN;                % ...otherwise, report NaNs
        subBin = NaN;        
    end
else 
    negZ = sortrows(negZ, -4);     % if >1 crossing, re-sort from largest to smallest peak...
    if log10(negZ(1, 4)) > minPow(negZ(1,2))        % ...if highest peak exceeds min threshold...
        if negZ(1, 4)*(1-mdiff) > negZ(2, 4)      % ...report frequency of this peak.
        	peakBin = negZ(1, 2);
            peak = negZ(1, 3); 
        else                        % ...if not...
            peak = NaN;
            subBin = negZ(1, 2);                    % ... index as a subpeak for starting alpha bound search.
        end
    else
        peak = NaN;                % ...otherwise, report NaNs
        subBin = NaN;
    end
end


% search for positive (upward going) zero-crossings (minima / valleys) either side of peak/subpeak(s)
slen = round(1/fres);               % define number of bins included in shollow slope search (approximate span = 1 Hz)
if isnan(peak) && isnan(subBin)       % if no evidence of peak activity, no parameter estimation indicated
    
    pos1 = NaN;
    pos2 = NaN;
    % f1 = NaN;
    % f2 = NaN;
    % inf1 = NaN;
    % inf2 = NaN;
    % Q = NaN;
    % Qf = NaN;

elseif isnan(peak)     % deal with spectra lacking a clear primary peak (similar strategy to peak; take highest subpeak as start point, look for minima)
    
    [f1, pos1] = findF1(f, d0, d1, negZ, minPow, slen, subBin); 
    [f2, pos2] = findF2(f, d0, d1, negZ, minPow, slen, subBin); 
        
    % inflections / Q values not calculated as these spectra won't be included in averaged channel peak analyses 
    % inf1 = NaN;     
    % inf2 = NaN;
    % Q = NaN;
    % Qf = NaN;

else            % now for the primary peak spectra
    
    [f1, pos1] = findF1(f, d0, d1, negZ, minPow, slen, peakBin);  
    [f2, pos2] = findF2(f, d0, d1, negZ, minPow, slen, peakBin); 
    
    % define boundaries by inflection points (requires 2nd derivative of smoothed signal)
    % inf1 = zeros(1,2);                  % initialise for zero-crossing count & frequency
    % cnt = 0;                            % start counter at 0
    % for k = 1:peakBin-1                 % step through frequency bins prior peak
    %     if sign(d2(k)) > sign(d2(k+1))                  % look for switch from positive to negative derivative values (i.e. downward zero-crossing)
    %         [~, mink] = min(abs([d2(k), d2(k+1)]));     % ensure correct frequency bin is picked out (find smaller of two values either side of crossing)
    %         if mink == 1
    %             min1 = k;
    %         else
    %             min1 = k+1;
    %         end
    %         cnt = cnt+1;                % advance counter by 1
    %         inf1(cnt,1) = cnt;          % zero-crossing count
    %         inf1(cnt,2) = f(min1);      % zero-crossing frequency
    %     end
    % end
    % 
    % % sort out appropriate estimates for output
    % if size(inf1, 1) == 1               % if singular crossing --> report frequency
    %     inf1 = inf1(1, 2);
    % else
    %     inf1 = sortrows(inf1, -2);      % sort by frequency values (descending)...
    %     inf1 = inf1(1, 2);              % take highest frequency (bin nearest to peak)
    % end
    % 
    % for k = peakBin+1:length(d2)-1                      % step through frequency bins post peak
    %     if sign(d2(k)) < sign(d2(k+1))                  % look for upward zero-crossing
    %         [~, mink] = min(abs([d2(k), d2(k+1)]));     % ensure frequency bin nearest zero-crossing point picked out (find smaller of two values either side of crossing)
    %             if mink == 1
    %                 min2 = k;
    %             else
    %                 min2 = k+1;
    %             end
    %             inf2 = f(min2);         % zero-crossing frequency
    %             break                   % break loop (only need to record first crossing)
    % 
    %     end
    % 
    % end
    % 
    % % estimate approx. area under curve between inflection points either
    % % side of peak, scale by inflection band width 
    % Q = trapz(f(min1:min2), d0(min1:min2));
    % Qf = Q / (min2-min1);

end 

bounds = [pos1 pos2];
