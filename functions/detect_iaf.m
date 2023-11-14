% Detect the individual alpha frequency (IAF) using 2 methods:
% 1) the peak alpha frequency (PAF) or 2) the alpha center of gravity (CoG).
%
% Please cite and see more detail here:
%   Corcoran et al., (2018). Toward a reliable, automated method of
%   individual alpha frequency (IAF) quantification.
%
% Cedric Cannard, 2023

function iaf = detect_iaf(pwr, freqs, winLength, params)
% [pSpec.sums, pSpec.chans, freqs] = restingIAF(pwr, nchan, cmin, fRange, fs, alphaRange, Fw, k);

fs = params.fs;
winLength = winLength*fs;   % window lenth in sample --> was tlen
alphaRange = [7 14];        % alpha range to search (Hz) --> was w
minChan = 3;                % min number of channels for cross-channel averages (between 1 and 6) --> was cmin
Fw = 11;                    % SGF frame width (11 corresponding to a frequency span of ~2.69 Hz @ ~.24Hz frequency resolution)
poly = 5;                   % SGF polynomial order --> was k
mpow = 1;
fRes = freqs(2)-freqs(1);   % Frequency resolution
% exp_tlen = nextpow2(winLength); fRes = fs/2.^exp_tlen;  % freq resolution

% PSD
% fh = str2func(taper);
% [pxx, f] = pwelch(data(kx,:), fh(tlen), tover, nfft, Fs);
% frex = dsearchn(f, fRange(1)):dsearchn(f, fRange(2));
% f = f(frex); pxx = pxx(frex);

% Calculate minPower vector
[pfit, sig] = polyfit(freqs, log10(pwr), 1);    % fit 1st order poly (regression line) to normalised spectra (log-scaled)
[yval, del] = polyval(pfit, freqs, sig);        % derive yval coefficients of fitted polynomial and delta (std dev) error estimate
minPow = yval + (mpow * del);                   % takes [minPowThresh * Std dev] as upper error bound on background spectral noise

%% Savitzky-Golay filter to fit curves to spectra & estimate 1st and 2nd
% derivatives Savitzky-Golay Smoothing and differentiation filter for
% fitting curves to spectra and extracting estimates of the 0th (smoothed),
% 1st, & 2nd derivative function of an input PSD. Savitzky-Golay filters
% are optimal in the sense that they minimize the least-squares error in
% fitting a polynomial to frames of noisy data.
[~, g] = sgolay(poly, Fw);
dt = fs/winLength;
dx = zeros(length(pwr),3);
for p = 0:2           % order of estimated derivatives
    dx(:,p+1) = conv(pwr, factorial(p)/(-dt)^p * g(:,p+1), 'same');
end
d0 = dx(:,1);        % smoothed signal post S-G diff filt
d1 = dx(:,2);        % 1st derivative
d2 = dx(:,3);        % 2nd derivative

%% Find peak(s) and boundaries of alpha band from derivatives
% Take derivatives from Savitzky-Golay curve-fitting and differentiation
% function sgfDiff, pump out estimates of alpha-band peak & bounds.
% Also calculates primary peak area Qf via integration between inflections.
%
% Depends on `findF1`, `findF2`, and `lessThan1` functions to locate
% bounds of individual alpha band.
%
% Outputs:
%   peakF = peak frequency estimate
%   posZ1 = freq of 1st positive zero-crossing (lower bound alpha interval)
%   posZ2 = freq of 2nd positive zero-crossing (upper bound alpha interval)
%   f1 = freq bin for posZ1
%   f2 = freq bin for posZ2
%   inf1 = inflection point, ascending edge
%   inf2 = inflection point, descending edge
%   Q = area under peak between inf1 & inf2
%   Qf = Q divided by bandwidth of Q

% [peakF, pos1, pos2, f1, f2, inf1, inf2, Q, Qf] = peakBounds(d0, d1, d2, freqs,fRange,minPow,mdiff,fRes);
% function [peakF, posZ1, posZ2, f1, f2, inf1, inf2, Q, Qf] = peakBounds(d0, d1, d2, freqs, alphaRange, minPow, minDiff, fRes)

% evaluate derivative for zero-crossings
% [~, lower_alpha] = min(abs(freqs-alphaRange(1)));      % set lower bound for alpha band
% [~, upper_alpha] = min(abs(freqs-alphaRange(2)));      % set upper bound for alpha band

negZ = zeros(1,4);                              % initialise for zero-crossing count & frequency bin
cnt = 0;                                        % start counter at 0
for k = alphaRange(1)-1:alphaRange(2)+1             % step through frequency bins in alpha band (start/end at bound -/+ 1 to make sure don't miss switch)
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
        negZ(cnt,3) = freqs(maxim);     % zero-crossing frequency
        negZ(cnt,4) = d0(maxim);    % power estimate
    end
end

% Sort out appropriate estimates for output
% if no zero-crossing detected --> report NaNs
if negZ(1,1) == 0
    peakF = NaN;
    subBin = NaN;

    % if singular crossing
elseif size(negZ, 1) == 1

    % if peak power is > minimum threshold --> report frequency
    if log10(negZ(1, 4)) > minPow(negZ(1,2))
        peakBin = negZ(1, 2);
        peakF = negZ(1, 3);
    else
        peakF = NaN;
        subBin = NaN;
    end

    % if >1 crossing, re-sort from largest to smallest peak
else
    negZ = sortrows(negZ, -4);

    % if highest peak exceeds min threshold
    if log10(negZ(1, 4)) > minPow(negZ(1,2))

        % report frequency of this peak.
        if negZ(1, 4)*(1-minDiff) > negZ(2, 4)
            peakBin = negZ(1, 2);
            peakF = negZ(1, 3);
        else
            peakF = NaN;
            subBin = negZ(1, 2);  % index as a subpeak for starting alpha bound search.
        end
    else
        peakF = NaN;
        subBin = NaN;
    end
end

%% Search for positive (upward going) zero-crossings (minima / valleys) either side of peak/subpeak(s)

% define number of bins included in shollow slope search (approximate span = 1 Hz)
slen = round(1/fRes);

% if no evidence of peak activity, no parameter estimation indicated
if isnan(peakF) && isnan(subBin)
    posZ1 = NaN;
    posZ2 = NaN;
    f1 = NaN;
    f2 = NaN;
    inf1 = NaN;
    inf2 = NaN;
    Q = NaN;
    Qf = NaN;

    % Deal with spectra lacking a clear primary peak (similar strategy to peak; take highest subpeak as start point, look for minima)
elseif isnan(peakF)

    [f1, posZ1] = findF1(freqs, d0, d1, negZ, minPow, slen, subBin);
    [f2, posZ2] = findF2(freqs, d0, d1, negZ, minPow, slen, subBin);

    % inflections / Q values not calculated as these spectra won't be included in averaged channel peak analyses
    inf1 = NaN;
    inf2 = NaN;
    Q = NaN;
    Qf = NaN;

    % now for the primary peak spectra
else

    [f1, posZ1] = findF1(freqs, d0, d1, negZ, minPow, slen, peakBin);
    [f2, posZ2] = findF2(freqs, d0, d1, negZ, minPow, slen, peakBin);

    % define boundaries by inflection points (requires 2nd derivative of smoothed signal)
    inf1 = zeros(1,2);                  % initialise for zero-crossing count & frequency
    cnt = 0;                            % start counter at 0

    % step through frequency bins prior peak
    for k = 1:peakBin-1

        % look for switch from positive to negative derivative values (i.e. downward zero-crossing)
        if sign(d2(k)) > sign(d2(k+1))

            % ensure correct frequency bin is picked out (find smaller of two values either side of crossing)
            [~, mink] = min(abs([d2(k), d2(k+1)]));
            if mink == 1
                min1 = k;
            else
                min1 = k+1;
            end
            cnt = cnt+1;                    % counter
            inf1(cnt,1) = cnt;              % zero-crossing count
            inf1(cnt,2) = freqs(min1);      % zero-crossing frequency
        end
    end

    % sort out appropriate estimates for output
    if size(inf1, 1) == 1               % if singular crossing --> report frequency
        inf1 = inf1(1, 2);
    else
        inf1 = sortrows(inf1, -2);      % sort by frequency values (descending)...
        inf1 = inf1(1, 2);              % take highest frequency (bin nearest to peak)
    end

    for k = peakBin+1:length(d2)-1                      % step through frequency bins post peak
        if sign(d2(k)) < sign(d2(k+1))                  % look for upward zero-crossing
            [~, mink] = min(abs([d2(k), d2(k+1)]));     % ensure frequency bin nearest zero-crossing point picked out (find smaller of two values either side of crossing)
            if mink == 1
                min2 = k;
            else
                min2 = k+1;
            end
            inf2 = freqs(min2);         % zero-crossing frequency
            break                   % break loop (only need to record first crossing)

        end

    end

    % estimate approx. area under curve between inflection points either
    % side of peak, scale by inflection band width
    Q = trapz(freqs(min1:min2), d0(min1:min2));
    Qf = Q / (min2-min1);

end



%% Estimate gravities for smoothed spectra (average IAF window across channels)
% Takes smoothed channel spectra and associated estimates of individual
% alpha bandwidth [f1:f2], calculate mean bandwidth, estimate CoG across
% all channels (as per Klimesch's group; e.g, 1990, 1993, & 1997 papers).
% Outputs:
%   cogs = centre of gravity derived from averaged f1:f2 frequency window
%   sel = channels contributing estimates of alpha window bandwidth
%   iaw = bounds of individual alpha window

% [ gravs, selG, iaw ] = chanGravs(d0, freqs, f1, f2);
% function [cogs, sel, iaw] = chanGravs(d0, freqs, f1, f2)
% trim off any NaNs in f1/f2 vectors
trim_f1 = f1(~isnan(f1));
trim_f2 = f2(~isnan(f2));

% derive average f1 & f2 values across chans, then look for nearest freq bin
mean_f1 = dsearchn(freqs, mean(freqs(trim_f1)));
mean_f2 = dsearchn(freqs, mean(freqs(trim_f2)));
iaw = [mean_f1, mean_f2];

% calculate CoG for each channel spectra on basis of averaged alpha window
cogs = zeros(1,size(d0,2));
for d = 1:length(cogs)
    if isempty(trim_f1) || isempty(trim_f2)
        cogs(d) = NaN;
    else
        cogs(d) = nansum(d0(mean_f1:mean_f2,d).*freqs(mean_f1:mean_f2)) / sum(d0(mean_f1:mean_f2,d));
    end
end

% report which channels contribute to averaged window
sel = ~isnan(f1);


%% calculate average pt estimates/spectra across k-th channels for each j-th recording
% Takes channel-wise estimates of peak alpha frequency (PAF) / centre of
% gravity (CoG) and calculates mean and standard deviation if cmin
% condition satisfied. PAFs are weighted in accordance with qf, which aims
% to quantify the relative strength of each channel peak.
%
% Outputs:
%   selP = channels contributing peak estimates to calculation of mean PAF
%   sums = structure containing summary estimates (m, std) for PAF and CoG

% [ selP, pSum ] = chanMeans(gravs, selG, [pChans(:).peakF], [pChans(:).d0], [pChans(:).Qf], cmin);
% function [selP, sums] = chanMeans(chanCogs, selG, peakF, specs, qf, cmin)

% function [selP, sums] = chanMeans(chanCogs, selG, peakF, specs, qf, cmin)
selP = ~isnan(peakF);               % evaluate whether channel provides estimate of PAF

% qWt = nansum(qf)/sum(selP);       % average area under peak (Qf) across viable channels (depricated: was used when calculating cross-recording comparisons)
chanWts = qf/max(qf);               % channel weightings scaled in proportion to Qf value of channel manifesting highest Qf

% average across peakF
if sum(selP) < cmin                 % if number of viable channels < cmin threshold, don't calculate cross-channel mean & std
    sums.paf = NaN;
    sums.pafStd = NaN;
    sums.muSpec = NaN;
else                                % else compute (weighted) cross-channel average PAF estimate and corresponding std of channel PAFs
    sums.paf = nansum(bsxfun(@times, peakF, chanWts))/nansum(chanWts);
    sums.pafStd = nanstd(peakF);
    % estimate averaged spectra for plotting
    wtSpec = bsxfun(@times, specs, chanWts);
    sums.muSpec = nansum(wtSpec, 2)/nansum(chanWts);
end

% now for the gravs (no channel weighting, all channels included if cmin satisfied)
if sum(selG) < cmin
    sums.cog = NaN;
    sums.cogStd = NaN;
else
    sums.cog = nanmean(chanCogs, 2);
    sums.cogStd = nanstd(chanCogs);
end

%% retain gravity estimates and selected channels in channel data struct
% (only loop through trimmed channels)
for cx = 1:nchan
    pChans(cx).gravs = gravs(cx);
    pChans(cx).selP = selP(cx);
    pChans(cx).selG = selG(cx);
end

%% get total number of chans that contributed to PAF/CoG estimation
pSum.pSel = sum(selP);
pSum.gSel = sum(selG);
pSum.iaw = iaw;


