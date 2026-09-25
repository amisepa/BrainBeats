% FINDF1 - Lower bound (f1) of the individual alpha band: searches the 1st
% derivative for a local minimum or near-horizontal slope before the alpha peak.
%
% Standalone copy of the findF1 subfunction of restingIAF, used by
% get_freqBounds (restingIAF uses its own copy), with
% lessThan1 as a local function. Returns NaN when no bound is found.
%
% Usage:
%   [f1, posZ1] = findF1(f, d0, d1, negZ, minPow, slen, bin)
%
% Inputs:
%   f      - frequency bin vector (Hz)
%   d0     - smoothed PSD estimate vector
%   d1     - 1st derivative of d0
%   negZ   - negative zero-crossings (peaks) found in the search window
%            [count, bin, frequency, power] (one row per peak)
%   minPow - minimum power threshold per bin (log10) defining candidate peaks
%   slen   - number of bins examined for a shallow slope (~1 Hz)
%   bin    - frequency bin of the peak / subpeak
%
% Outputs:
%   f1     - frequency bin of the lower bound (NaN if not found)
%   posZ1  - frequency of the lower bound (Hz; NaN if not found)
%
% Part of the restingIAF package, (c) Andrew W. Corcoran, 2016-2018
% (Corcoran et al. 2018, Psychophysiology). See github.com/corcorana/restingIAF.
function [f1, posZ1] = findF1(f, d0, d1, negZ, minPow, slen, bin)

posZ1 = zeros(1,4);

% contingency for multiple peaks - aim to identify left-most peak in range for upper bound of k in next loop (avoid falling into local minimum)
if size(negZ, 1) >1
    negZ = sortrows(negZ, 3);       % sort by frequency (ascending)
    for z = 1:size(negZ, 1)
        if log10(negZ(z, 4)) > minPow(negZ(1,2)) || negZ(z, 4) > (0.5* d0(bin)) % relax power constraint, as may result in excessively narrow alpha window in noisy spectra with shallow peakF (i.e. precisely where we want CoG to take breadth into account)
            leftPeak = negZ(z, 2);
            break                   % break off search when conditions satisfied
        else 
            leftPeak = bin;        % if fail to satisfy search conditions, default to bin (sub)peak
        end
    end
else 
    leftPeak = bin;                % if no other peaks were identified, take bin (sub)peak as boundary
end

cnt = 0;                        % start counter at 0
for k = 2:leftPeak-1            % step through frequency bins up to left-most peak in search window
    if sign(d1(k)) < sign(d1(k+1))        % look for switch from negative to positive derivative values (i.e. upward/positive zero-crossing)
        [~, mink] = min(abs([d0(k-1), d0(k), d0(k+1)]));    % search around crossing for local minimum in d0 (indexing 1st derivative sometimes results in small errors)
        if mink == 1
            minim = k-1;
        elseif mink == 2
            minim = k;
        else
            minim = k+1;
        end

        cnt = cnt+1;                % advance counter by 1
        posZ1(cnt,1) = cnt;         % zero-crossing count
        posZ1(cnt,2) = minim;       % zero-crossing frequency bin
        posZ1(cnt,3) = f(minim);    % zero-crossing frequency

       	% look for consistent low d1 values for signs of shallow slope (levelling off)
    elseif abs(d1(k)) < 1 && lessThan1(d1(k+1:min(k+slen,length(d1))))
        minim = k;
        cnt = cnt+1;                % advance counter by 1
        posZ1(cnt,1) = cnt;         % zero-crossing count
        posZ1(cnt,2) = minim;       % zero-crossing frequency bin
        posZ1(cnt,3) = f(minim);    % zero-crossing frequency
    end

end

% sort out appropriate estimates for output
if cnt == 0                         % if no crossing or shallow slope found --> report NaNs
    f1 = NaN;
    posZ1 = NaN;
elseif size(posZ1, 1) == 1          % if singular crossing --> report frequency
    f1 = posZ1(1, 2);
    posZ1 = posZ1(1, 3);
else                                % else sort by frequency values (descending), take highest frequency (bin nearest to peak)
    posZ1 = sortrows(posZ1, -3);
    f1 = posZ1(1, 2);
    posZ1 = posZ1(1, 3);
end

%% Subfunction

function tval = lessThan1(d1)
% True if all values of a 1st derivative segment (~1 Hz) are within +/- 1
% (shallow slope). Copy of the lessThan1 subfunction of restingIAF.
if length(d1) < 2
    error('Length of 1st derivative segment < 2');
end
tval = all(abs(d1) < 1);
