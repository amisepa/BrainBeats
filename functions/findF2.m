% FINDF2 - Upper bound (f2) of the individual alpha band: searches the 1st
% derivative for a local minimum or near-horizontal slope after the alpha peak.
%
% Standalone copy of the findF2 subfunction of restingIAF, used by
% get_freqBounds (restingIAF uses its own copy), with
% lessThan1 as a local function. Returns NaN when no bound is found.
%
% Usage:
%   [f2, posZ2] = findF2(f, d0, d1, negZ, minPow, slen, bin)
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
%   f2     - frequency bin of the upper bound (first one found; NaN if none)
%   posZ2  - frequency of the upper bound (Hz; NaN if none)
%
% Part of the restingIAF package, (c) Andrew W. Corcoran, 2016-2018
% (Corcoran et al. 2018, Psychophysiology). See github.com/corcorana/restingIAF.
function [f2, posZ2] = findF2(f, d0, d1, negZ, minPow, slen, bin)

posZ2 = zeros(1,4);

% contingency for multiple peaks - try to identify right-most peak in range for upper bound of k in next loop (avoid falling into local minima)
if size(negZ, 1) >1
    negZ = sortrows(negZ, -3);      % sort by frequency (descending)
    for z = 1:size(negZ, 1)
        if log10(negZ(z, 4)) > minPow(negZ(1,2)) || negZ(z, 4) > (0.5* d0(bin))
            rightPeak = negZ(z, 2);
            break                 	% break search when conditions satisfied
        else 
            rightPeak = bin;       % if fail to satisfy search conditions, default to bin (sub)peak
        end
    end
else 
    rightPeak = bin;               % if no other peaks identified, take bin (sub)peak as boundary
end

cnt = 0;                            % start counter at 0
for k = rightPeak+1:length(d1) - slen     % step through frequency bins following right-most peak (trim end of range to allow for following conditional search of d1 values < 1)
    if sign(d1(k)) < sign(d1(k+1))            % look for switch from negative to positive derivative values (i.e. upward/positive zero-crossing)
        [~, mink] = min(abs([d0(k-1), d0(k), d0(k+1)]));    % search around crossing for local minimum in d0 (indexing 1st derivative sometimes results in small errors)
        if mink == 1
            minim = k-1;
        elseif mink == 2
            minim = k;
        else
            minim = k+1;
        end

        cnt = cnt+1;                % advance counter by 1
        posZ2(cnt,1) = cnt;         % zero-crossing count
        posZ2(cnt,2) = minim;       % zero-crossing frequency bin
        posZ2(cnt,3) = f(minim);    % zero-crossing frequency

        % look for consistent low d1 values for signs of shallow slope (levelling off)
    elseif abs(d1(k)) < 1 && lessThan1(d1(k+1:k+slen))
        minim = k;
        cnt = cnt+1;                % advance counter by 1
        posZ2(cnt,1) = cnt;         % zero-crossing count
        posZ2(cnt,2) = minim;       % zero-crossing frequency bin
        posZ2(cnt,3) = f(minim);    % zero-crossing frequency
    end
end

if cnt == 0                         % if no crossing or shallow slope found --> report NaNs
    f2 = NaN;
    posZ2 = NaN;
else
    f2 = posZ2(1, 2);
    posZ2 = posZ2(1, 3);            % can simply take first estimate for output
end

%% Subfunction

function tval = lessThan1(d1)
% True if all values of a 1st derivative segment (~1 Hz) are within +/- 1
% (shallow slope). Copy of the lessThan1 subfunction of restingIAF.
if length(d1) < 2
    error('Length of 1st derivative segment < 2');
end
tval = all(abs(d1) < 1);
