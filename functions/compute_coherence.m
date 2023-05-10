%% get_coherence
% First compute power-spectrum for each channel and the cross-spectrum
% on each epoch using fft, and then computes the coherence. 
%  
% Inputs:
%   - data1 and data2: continuous EEG data (e.g., 2 different electrodes)
%   - sRate: data sampling rate in Hz
%   - wSize: window size in s (default = 2)
%   - fInt: vector of frequencies of interest for coherence output; default = up to Nyquist frequency)
% 
% Outputs:
%   - c_mu: mean coherence score (from 0 to 1)
%   - c: coherence for each sliding window (from 0 to 1)
%   - f: frequency vector
% 
% Example: 
%   [c_mu, c, f] = get_coherence(data1, data2, EEG.srate, 2, 1:100)
%   figure; plot(f,c); xlabel('Frequencies (Hz)'); ylabel('Coherence');
%   or 
%   c = get_coherence(data1, data2, EEG.srate, 2, 10)
% 
% Cedric Cannard, January 2022

function [c_mu, c, f] = get_coherence(data1, data2, sRate, wSize)

dt = 1/sRate;   % discrete sampling in s (0.002 = 2 ms)
df = 1/wSize;   % frequency resolution
fNQ = 1/dt/2;   % Nyquist frequency
f = df:df:fNQ;  % whole frequency vector availabe in data

% %check frequency vector    
% if fInt(end)> f(end)
%     warning('Upper frequency bound selected is above the Nyquist limit. Changing to Nyquist frequency.');
%     fInt = fInt(1):f(end);
% end

% Initiate variables
nWind = floor((size(data1,2)/sRate)/wSize);
Sxx = zeros(1,sRate*wSize);
Syy = zeros(1,sRate*wSize);
Sxy = zeros(1,sRate*wSize);

% Compute power spectra and cross-spectrum for each sliding window
window = 1;
for iWind = 1:wSize*sRate:nWind*sRate%-wSize*sRate
    Sxx(window,:) = 2*dt^2/wSize * fft(data1(iWind:iWind+wSize*sRate-1)).*conj(fft(data1(iWind:iWind+wSize*sRate-1)));
    Syy(window,:) = 2*dt^2/wSize * fft(data2(iWind:iWind+wSize*sRate-1)).*conj(fft(data2(iWind:iWind+wSize*sRate-1)));
    Sxy(window,:) = 2*dt^2/wSize * fft(data1(iWind:iWind+wSize*sRate-1)).*conj(fft(data2(iWind:iWind+wSize*sRate-1)));
    window = window+1;
end

% Keep only positive frequencies
Sxx = Sxx(:,1:size(Sxx,2)/2);
Syy = Syy(:,1:size(Syy,2)/2);
Sxy = Sxy(:,1:size(Sxy,2)/2);

% Compute coherence
c = abs( Sxy ./ ( sqrt(Sxx) .* sqrt(Syy) ) );
% c = Sxy ./ ( sqrt(Sxx) .* sqrt(Syy) );

% Mean coherence across windows
Sxx_mu = mean(Sxx,1);
Syy_mu = mean(Syy,1);
Sxy_mu = mean(Sxy,1);
c_mu = abs( Sxy_mu ./ (sqrt(Sxx_mu) .* sqrt(Syy_mu)) );

% % Coherence for frequencies of interest
% if size(fInt,2) > 1
%     f = find(f==fInt(1)):find(f==fInt(end));
%     c = c(:,f);
%     f = fInt(1):1/wSize:fInt(end);
% else
%     f = find(f==fInt); %if only one frequency was selected
%     c = c(:,f);
% end
