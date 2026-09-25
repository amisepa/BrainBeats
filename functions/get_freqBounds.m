% GET_FREQBOUNDS - Individualized bounds of a frequency band from the shape
% of the power spectrum (peak and surrounding minima within a search window).
%
% Adapts the restingIAF alpha-band method (Corcoran et al. 2018) to any band:
% the spectrum is normalized by its mean (as in restingIAF), smoothed and
% differentiated with a Savitzky-Golay filter, the highest peak above the
% log10 1/f fit (+ mpow SD) is searched within w, and the band bounds are
% the minima (or flattening of the slope) on either side of that peak.
% Used by the 'individualized' option of get_eeg_features.
%
% Usage:
%   [bounds, peak] = get_freqBounds(pwr, f, fs, w, winSize, mpow)
%
% Inputs:
%   pwr     - power spectrum of one channel, in uV^2/Hz (not normalized or in dB)
%   f       - frequencies of pwr (Hz)
%   fs      - sample rate (Hz)
%   w       - peak search window [low high] (Hz)
%   winSize - pwelch window length in samples (scales the S-G derivatives, dt = fs/winSize)
%   mpow    - threshold (SDs above the 1/f regression fit) for a peak to be
%             distinguished from background noise (e.g. 1 for alpha, 0.25 for other bands)
%
% Outputs:
%   bounds  - [lower upper] band bounds (Hz); NaN when no peak (or subpeak)
%             is detected in w, or when a bound is not found
%   peak    - peak frequency (Hz); NaN if no clear primary peak
%
% Reference:
%   Corcoran, A.W., Alday, P.M., Schlesewsky, M., & Bornkessel-Schlesewsky, I.
%   (2018). Toward a reliable, automated method of individual alpha frequency
%   (IAF) quantification. Psychophysiology, 55(7), e13064.
%
% Copyright (C) - Cedric Cannard, 2023

function [bounds, peak] = get_freqBounds(pwr, f, fs, w, winSize, mpow)

% normalize the spectrum by its mean (as restingIAF), so that the slope
% criteria in findF1/findF2 (|d1| < 1) do not depend on the power units
pwr = pwr(:) / mean(pwr);
f = f(:);

% frequency resolution
fres = f(2)-f(1);

% Parameters
Fw = 2*floor(2.69/fres/2) + 1;  % SGF frame width spanning ~2.69 Hz (11 bins @ ~.24Hz frequency resolution, restingIAF default)
Fw = max(Fw, 5);
k = min(5, Fw-2);   % SGF polynomial order (default = 5, must be < Fw)
mdiff = .2;    % minimal height difference distinguishing a primary peak from
                % competing peaks (default = 0.2; i.e. 20% peak height)

% fit 1st order poly (regression line) to the log10 spectrum
[pfit, sig] = polyfit(f, log10(pwr), 1);

% derive yval coefficients of fitted polynomial and delta (std dev) error estimate
[yval, del] = polyval(pfit, f, sig);

% take [minPowThresh * Std dev] as upper error bound on background spectral
% noise (log10 units, compared with log10 of the peak power below)
minPow = yval + (mpow * del);

% apply Savitzky-Golay filter to fit curves to spectra & estimate 1st and 2nd derivatives
% (sgfDiff in restingIAF)
[~, g] = sgolay(k, Fw);      
dt = fs/winSize;
dx = zeros(length(pwr),3);
for p = 0:2                 % p determines order of estimated derivatives
    dx(:,p+1) = conv(pwr, factorial(p)/(-dt)^p * g(:,p+1), 'same');
end
d0 = dx(:,1);        % smoothed signal post S-G diff filt
d1 = dx(:,2);        % 1st derivative
d2 = dx(:,3);        % 2nd derivative

% find peak(s) and boundaries of the frequency band from the derivatives
% (peakBounds in restingIAF, without the inflection points and Q)

% evaluate derivative for zero-crossings
[~, lower_lim] = min(abs(f-w(1)));      % set lower bound of the search window
[~, upper_lim] = min(abs(f-w(2)));      % set upper bound of the search window

negZ = zeros(1,4);  % zero-crossing count & frequency bin
cnt = 0;
for k = max(lower_lim-1,1):min(upper_lim+1,length(d1)-1)   % step through frequency bins in the window (start/end at bound -/+ 1 to make sure don't miss switch, within the spectrum)
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
            subBin = negZ(1, 2);                    % ... index as a subpeak for starting bound search.
        end
    else
        peak = NaN;                % ...otherwise, report NaNs
        subBin = NaN;
    end
end


% search for positive (upward going) zero-crossings (minima / valleys) either side of peak/subpeak(s)
slen = round(1/fres);               % define number of bins included in shallow slope search (approximate span = 1 Hz)
if isnan(peak) && isnan(subBin)       % if no evidence of peak activity, no parameter estimation indicated

    pos1 = NaN;
    pos2 = NaN;

elseif isnan(peak)     % deal with spectra lacking a clear primary peak (similar strategy to peak; take highest subpeak as start point, look for minima)

    [f1, pos1] = findF1(f, d0, d1, negZ, minPow, slen, subBin);
    [f2, pos2] = findF2(f, d0, d1, negZ, minPow, slen, subBin);

else            % now for the primary peak spectra

    [f1, pos1] = findF1(f, d0, d1, negZ, minPow, slen, peakBin);
    [f2, pos2] = findF2(f, d0, d1, negZ, minPow, slen, peakBin);

end

bounds = [pos1 pos2];
