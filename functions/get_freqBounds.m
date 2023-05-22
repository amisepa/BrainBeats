function [bounds, peak] = get_freqBounds(pwr, f, fs, w, winSize, norm)

% Parameters
% cmin = 1;     % min number of channel for cross-channel averages (1-6)
% fRange = [1 40];    % spectral range (set to filter passband)
% w = [13 30];         % alpha peak search window (Hz)
Fw = 11;            % SGF frame width (11 corresponding to a frequency span of ~2.69 Hz @ ~.24Hz frequency resolution)
k = 5;              % SGF polynomial order (default = 5)
mdiff = .20;    % minimal height difference distinguishing a primary peak from
                % competing peaks (default = 0.20; i.e. 20% peak height)
mpow = 1;       % error bound (SD) threshold to differentiate peaks from background spectral noise (default = 1)

% normalize power (default = true)
if norm
    pwr = pwr / mean(pwr); 
end

% fit 1st order poly (regression line) to normalised spectra (log-scaled)
[pfit, sig] = polyfit(f, log(pwr), 1);     

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
[peak, pos1, pos2, f1, f2, inf1, inf2, Q, Qf] = peakBounds(d0, d1, d2, f, w, minPow, mdiff, fres);
bounds = [pos1 pos2];

