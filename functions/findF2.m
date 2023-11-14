% Search 1st derivative for evidence of local minima or near horizontal
% function post alpha peak. This location will be taken as the lower
% bound of the individual alpha band used to calculate CoG (f2).
function [f2, posZ2] = findF2(f, d0, d1, negZ, minPow, slen, bin)

posZ2 = zeros(1,4);

% contingency for multiple peakF - try to identify right-most peak in range for upper bound of k in next loop (avoid falling into local minima)
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
    rightPeak = bin;               % if no other peakF identified, take bin (sub)peak as boundary
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

f2 = posZ2(1, 2);
posZ2 = posZ2(1, 3);                % can simply take first estimate for output
