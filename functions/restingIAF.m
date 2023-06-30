function [pSum, pChans, f] = restingIAF(data, nchan, cmin, fRange, Fs, w, Fw, k, varargin)
% Primary function for running `restingIAF` analysis routine for estimating
% two indices of individual alpha frequency (IAF): Peak alpha frequency
% (PAF) and the alpha centre of gravity (CoG) or mean frequency.
%
% Calls on the Signal Processing Toolbox function `pwelch` to derive power
% spectral density estimates of one or more vectors of EEG channel data,
% which are subsequently smoothed by `sgfDiff` (depends on `sgolay`).
% Alpha peak activity is parameterised by `peakBounds` and CoG estimated by
% `chanGravs`. Channel-wise peak and gravity estimates are averaged across
% channels by `chanMeans`.
%
% This function and all custom-designed dependencies are part of the
% `restingIAF` package, (c) Andrew W. Corcoran, 2016-2018.
%
% Please consult our methods paper for a more detailed exposition of the
% analysis routine, factors to consider when selecting parameter settings,
% and a study of its performance on empirical and simulated EEG signals:
%
% Corcoran, A.W., Alday, P.M., Schlesewsky, M., & Bornkessel-Schlesewsky,
%   I. (2018). Toward a reliable, automated method of individual alpha
%   frequency (IAF) quantification. Psychophysiology, 55(7), e13064.
%   doi: 10.1111/psyp.13064.
%
% Each release version is also citable and archived on GitHub, making it
% easier for others to fully replicate your analysis (see README.md).
%
% Visit github.com/corcorana/restingIAF for further info on licencing and
% updates on package development.
%
% Outputs:
%   pSum    = structure containing summary statistics of alpha-band parameters
%   pChans  = structure containing channel-wise spectral and alpha parameter data
%   f       = trimmed vector of frequency bins resolved by `pwelch`
%
% Inputs:
%   data    = vector or matrix containing continuous EEG channel data
%             (matrix rows = channels, cols = sample points)
%   nchan   = number of channels in data array
%   cmin    = minimum number of channel estimtes that must be resolved in
%             order to calculate average PAF/CoG estimates
%   fRange  = frequency range to be included in analysis (e.g., [1, 40] Hz)
%   Fs      = EEG sampling rate
%   w       = bounds of alpha peak search window (e.g., [7 13])
%   Fw      = frame width, Savitzky-Golay filter (corresponds to number of
%             freq. bins spanned by filter; must be odd)
%   k       = polynomial order, Savitzky-Golay filter (must be < Fw)
%
% Optional inputs:
%   mpow    = error bound (s.d.) used to determine threshold differentiating
%             substantive peaks from background spectral noise (default = 1)
%   mdiff   = minimal height difference distinguishing a primary peak from
%             competing peaks (default = 0.20; i.e. 20% peak height)
%   taper   = taper window function applied by `pwelch` (default = 'hamming')
%   tlen    = length of taper window applied by `pwelch` (default = 4 sec)
%   tover   = length of taper window overlap in samples (default = 50% window length)
%   nfft    = specify number of FFT points used to calculate PSD (default =
%             next power of 2 above window length)
%   norm    = normalise power spectra (default = true)
%
% setup inputParser
p = inputParser;
p.addRequired('data',...
    @(x) validateattributes(x, {'numeric'}, ...
    {'2d', 'nonempty'}));
p.addRequired('nchan',...
    @(x) validateattributes(x, {'numeric'}, ...
    {'scalar', 'integer', 'positive'}));
p.addRequired('cmin',...
    @(x) validateattributes(x, {'numeric'}, ...
    {'scalar', 'integer', 'positive', '<=', size(data, 1)}));
p.addRequired('fRange',...
    @(x) validateattributes(x, {'numeric'}, ...
    {'integer', 'nonnegative', 'increasing', 'size', [1,2]}));
p.addRequired('Fs',...
    @(x) validateattributes(x, {'numeric'}, ...
    {'scalar', 'integer', 'positive', '>=', 2*fRange(2)}));
p.addRequired('w',...
    @(x) validateattributes(x, {'numeric'}, ...
    {'nonnegative', 'increasing', 'size', [1,2], '>', fRange(1), '<', fRange(2)}));
p.addRequired('Fw',...
    @(x) validateattributes(x, {'numeric'}, ...
    {'scalar', 'integer', 'positive', 'odd'}));
p.addRequired('k',...
    @(x) validateattributes(x, {'numeric'}, ...
    {'scalar', 'integer', 'positive', '<', Fw }));

p.addOptional('mpow', 1,...
    @(x) validateattributes(x, {'numeric'}, ...
    {'scalar', 'positive'}));
p.addOptional('mdiff', .20,...
    @(x) validateattributes(x, {'numeric'}, ...
    {'scalar', '>=', 0, '<=', 1}));
p.addOptional('taper', 'hamming',...
    @(x) validateattributes(x, {'char'}, ...
    {}));
p.addOptional('tlen', (Fs*4),...
    @(x) validateattributes(x, {'numeric'}, ...
    {'scalar', 'integer', 'positive'}));
p.addOptional('tover', [],...
    @(x) validateattributes(x, {'numeric'}, ...
    {'integer', 'nonnegative'}));
p.addOptional('nfft', [],...
    @(x) validateattributes(x, {'numeric'}, ...
    {'integer', 'positive'}));
p.addOptional('norm', true,...
    @(x) validateattributes(x, {'logical'}, ...
    {'scalar'}));

p.parse(data, nchan, cmin, fRange, Fs, w, Fw, k, varargin{:})

mpow    = p.Results.mpow;
mdiff   = p.Results.mdiff;
taper   = p.Results.taper;
tlen    = p.Results.tlen;
tover   = p.Results.tover;
nfft    = p.Results.nfft;
norm    = p.Results.norm;

% struct for channel data (PSD estimates & derivatives, some additional info)
pChans = struct('pxx', [], 'minPow', [], 'd0', [], 'd1', [], 'd2', [],...
    'peaks', [], 'pos1', [], 'pos2', [], 'f1', [], 'f2', [], 'inf1', [], 'inf2', [],...
    'Q', [],'Qf', [], 'gravs', [], 'selP', [], 'selG', [] );


for kx = 1:nchan
    if sum(isnan(data(kx,:)))==0      % ensure no channel NaNs

        % perform pwelch routine to extract PSD estimates by channel
        fh = str2func(taper);
        [pxx, f] = pwelch(data(kx,:), fh(tlen), tover, nfft, Fs);

        % delimit range of freq bins to be included in analysis
        frex = dsearchn(f, fRange(1)):dsearchn(f, fRange(2));
        f = f(frex);
        pxx = pxx(frex);      % truncate PSD to frex range

        % normalise truncated PSD
        if norm == true
            pChans(kx).pxx = pxx / mean(pxx);
        else
            pChans(kx).pxx = pxx;
        end

        % calculate minPower vector
        [pfit, sig] = polyfit(f, log10(pChans(kx).pxx), 1);     % fit 1st order poly (regression line) to normalised spectra (log-scaled)
        [yval, del] = polyval(pfit, f, sig);                    % derive yval coefficients of fitted polynomial and delta (std dev) error estimate
        pChans(kx).minPow = yval + (mpow * del);                % takes [minPowThresh * Std dev] as upper error bound on background spectral noise

        % apply Savitzky-Golay filter to fit curves to spectra & estimate 1st and 2nd derivatives
        [pChans(kx).d0, pChans(kx).d1, pChans(kx).d2] = sgfDiff(pChans(kx).pxx, Fw, k, Fs, tlen);

        % calculate frequency resolution
        exp_tlen = nextpow2(tlen);
        fres = Fs/2.^exp_tlen;

        % take derivatives, find peak(s) and boundaries of alpha band
        [pChans(kx).peaks, pChans(kx).pos1, pChans(kx).pos2, pChans(kx).f1, pChans(kx).f2,...
            pChans(kx).inf1, pChans(kx).inf2, pChans(kx).Q, pChans(kx).Qf]...
            = peakBounds(pChans(kx).d0, pChans(kx).d1, pChans(kx).d2, f, w,...
            pChans(kx).minPow, mdiff, fres);

    else
        warning('Row #%s contains NaNs, skipping channel...', num2str(kx))
        nchan = nchan-1;    % trim nchans for cx loop later on
    end

end

% estimate gravities for smoothed spectra (average IAF window across channels)
[ gravs, selG, iaw ] = chanGravs([pChans(:).d0], f, [pChans(:).f1], [pChans(:).f2] );

% calculate average pt estimates/spectra across k-th channels for each j-th recording
[ selP, pSum ] = chanMeans(gravs, selG, [pChans(:).peaks], [pChans(:).d0], [pChans(:).Qf], cmin);

% retain gravity estimates and selected channels in channel data struct
% (only loop through trimmed channels)
for cx = 1:nchan
    pChans(cx).gravs = gravs(cx);
    pChans(cx).selP = selP(cx);
    pChans(cx).selG = selG(cx);
end

% get total number of chans that contributed to PAF/CoG estimation
pSum.pSel = sum(selP);
pSum.gSel = sum(selG);
pSum.iaw = iaw;

%% Subfunction 1

function [d0, d1, d2] = sgfDiff(x, Fw, poly, Fs, tlen)
% Savitzky-Golay Smoothing and Differentiation Filter for extracting
% estimates of the 0th (smoothed), 1st, & 2nd derivative function of an
% input PSD.
%
% "Savitzky-Golay filters are optimal in the sense that they minimize the
% least-squares error in fitting a polynomial to frames of noisy data."
%
% Depends on MATLAB Signal Processing Toolbox function to implement S-G
% filter.
%
% Part of the `restingIAF` package, (c) Andrew W. Corcoran, 2016-2017.
% Visit github.com/corcorana/restingIAF for licence, updates, and further
% information about package development and testing.
%
% Outputs:
%   d0 = smoothed PSD estimates
%   d1 = 1st derivative of d0
%   d2 = 2nd derivative of d0
%
% Inputs:
%   x = spectral data to be filtered/differentiated (vector)
%   poly = polynomial order (integer, must be < Fw)
%   Fw = frame width (i.e. number of samples, must be an odd integer)
%   Fs = sampling rate (integer)
%   tlen = taper length (i.e. number of samples of pwelch window, integer)

% setup inputParser
% p = inputParser;
% p.addRequired('x',...
%                 @(x) validateattributes(x, {'numeric'}, ...
%                 {'vector'}));
% p.addRequired('Fw',...
%                 @(x) validateattributes(x, {'numeric'}, ...
%                 {'scalar', 'integer', 'positive', 'odd'}));
% p.addRequired('poly',...
%                 @(x) validateattributes(x, {'numeric'}, ...
%                 {'scalar', 'integer', 'positive', '<', Fw }));
% p.addRequired('Fs',...
%                 @(x) validateattributes(x, {'numeric'}, ...
%                 {'scalar', 'integer', 'positive'}));
% p.addRequired('tlen',...
%                 @(x) validateattributes(x, {'numeric'}, ...
%                 {'scalar', 'integer', 'positive'}));
% p.parse(x, Fw, poly, Fs, tlen)

[~, g] = sgolay(poly, Fw);
dt = Fs/tlen;
dx = zeros(length(x),3);
for p = 0:2                 % p determines order of estimated derivatives
    dx(:,p+1) = conv(x, factorial(p)/(-dt)^p * g(:,p+1), 'same');
end
d0 = dx(:,1);        % smoothed signal post S-G diff filt
d1 = dx(:,2);        % 1st derivative
d2 = dx(:,3);        % 2nd derivative

%% Subfunction 2

function [peakF, posZ1, posZ2, f1, f2, inf1, inf2, Q, Qf] = peakBounds(d0, d1, d2, f, w, minPow, minDiff, fres)
% Take derivatives from Savitzky-Golay curve-fitting and differentiation
% function sgfDiff, pump out estimates of alpha-band peak & bounds.
% Also calculates primary peak area Qf via integration between inflections.
%
% Depends on `findF1`, `findF2`, and `lessThan1` functions to locate
% bounds of individual alpha band.
%
% Part of the `restingIAF` package, (c) Andrew W. Corcoran, 2016-2017.
% Visit github.com/corcorana/restingIAF for licence, updates, and further
% information about package development and testing.
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
%
% Required inputs:
%   d0 = smoothed PSD estimate vector
%   d1 = 1st derivative vector
%   d2 = 2nd derivative vector
%   f = frequency bin vector
%   w = bounds of initial alpha window
%   minPow = vector of minimum power threshold values defining candidate peaks (regression fit of background spectral activity)
%   minDiff = minimum difference required to distinguish peak as dominant (proportion of primary peak height)
%   fres = frequency resolution (determine how many bins to search to establish shallow rolloff in d1)

% evaluate derivative for zero-crossings
[~, lower_alpha] = min(abs(f-w(1)));      % set lower bound for alpha band
[~, upper_alpha] = min(abs(f-w(2)));      % set upper bound for alpha band

negZ = zeros(1,4);                              % initialise for zero-crossing count & frequency bin
cnt = 0;                                        % start counter at 0
for k = lower_alpha-1:upper_alpha+1             % step through frequency bins in alpha band (start/end at bound -/+ 1 to make sure don't miss switch)
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
    peakF = NaN;
    subBin = NaN;
elseif size(negZ, 1) == 1           % if singular crossing...
    if log10(negZ(1, 4)) > minPow(negZ(1,2))      % ...and peak power is > minimum threshold --> report frequency
        peakBin = negZ(1, 2);
        peakF = negZ(1, 3);
    else
        peakF = NaN;                % ...otherwise, report NaNs
        subBin = NaN;
    end
else
    negZ = sortrows(negZ, -4);     % if >1 crossing, re-sort from largest to smallest peak...
    if log10(negZ(1, 4)) > minPow(negZ(1,2))        % ...if highest peak exceeds min threshold...
        if negZ(1, 4)*(1-minDiff) > negZ(2, 4)      % ...report frequency of this peak.
            peakBin = negZ(1, 2);
            peakF = negZ(1, 3);
        else                        % ...if not...
            peakF = NaN;
            subBin = negZ(1, 2);                    % ... index as a subpeak for starting alpha bound search.
        end
    else
        peakF = NaN;                % ...otherwise, report NaNs
        subBin = NaN;
    end
end


% search for positive (upward going) zero-crossings (minima / valleys) either side of peak/subpeak(s)
slen = round(1/fres);    % define number of bins included in shollow slope search (approximate span = 1 Hz)

if isnan(peakF) && isnan(subBin)       % if no evidence of peak activity, no parameter estimation indicated

    posZ1 = NaN;
    posZ2 = NaN;
    f1 = NaN;
    f2 = NaN;
    inf1 = NaN;
    inf2 = NaN;
    Q = NaN;
    Qf = NaN;

elseif isnan(peakF)     % deal with spectra lacking a clear primary peak (similar strategy to peak; take highest subpeak as start point, look for minima)

    [f1, posZ1] = findF1(f, d0, d1, negZ, minPow, slen, subBin);
    [f2, posZ2] = findF2(f, d0, d1, negZ, minPow, slen, subBin);

    % inflections / Q values not calculated as these spectra won't be included in averaged channel peak analyses
    inf1 = NaN;
    inf2 = NaN;
    Q = NaN;
    Qf = NaN;

else            % now for the primary peak spectra

    [f1, posZ1] = findF1(f, d0, d1, negZ, minPow, slen, peakBin);
    [f2, posZ2] = findF2(f, d0, d1, negZ, minPow, slen, peakBin);

    % define boundaries by inflection points (requires 2nd derivative of smoothed signal)
    inf1 = zeros(1,2);                  % initialise for zero-crossing count & frequency
    cnt = 0;                            % start counter at 0
    for k = 1:peakBin-1                 % step through frequency bins prior peak
        if sign(d2(k)) > sign(d2(k+1))                  % look for switch from positive to negative derivative values (i.e. downward zero-crossing)
            [~, mink] = min(abs([d2(k), d2(k+1)]));     % ensure correct frequency bin is picked out (find smaller of two values either side of crossing)
            if mink == 1
                min1 = k;
            else
                min1 = k+1;
            end
            cnt = cnt+1;                % advance counter by 1
            inf1(cnt,1) = cnt;          % zero-crossing count
            inf1(cnt,2) = f(min1);      % zero-crossing frequency
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
            inf2 = f(min2);         % zero-crossing frequency
            break                   % break loop (only need to record first crossing)

        end

    end

    % estimate approx. area under curve between inflection points either
    % side of peak, scale by inflection band width
    Q = trapz(f(min1:min2), d0(min1:min2));
    Qf = Q / (min2-min1);

end

%% Subfunction 3

function [cogs, sel, iaw] = chanGravs(d0, f, f1, f2)
% Takes smoothed channel spectra and associated estimates of individual
% alpha bandwidth [f1:f2], calculate mean bandwidth, estimate CoG across
% all channels (as per Klimesch's group; e.g, 1990, 1993, & 1997 papers).
%
% Part of the `restingIAF` package, (c) Andrew W. Corcoran, 2016-2017.
% Visit github.com/corcorana/restingIAF for licence, updates, and further
% information about package development and testing.
%
% Outputs:
%   cogs = centre of gravity derived from averaged f1:f2 frequency window
%   sel = channels contributing estimates of alpha window bandwidth
%   iaw = bounds of individual alpha window
%
% Inputs:
%   d0 = vector / matrix of smoothed PSDs
%   f = frequency bin vector
%   f1 = vector of f1 bin indices (lower bound, individual alpha window)
%   f2 = vector of f2 bin indices (upper bound, individual alpha window)

% trim off any NaNs in f1/f2 vectors
trim_f1 = f1(~isnan(f1));
trim_f2 = f2(~isnan(f2));

% derive average f1 & f2 values across chans, then look for nearest freq bin
mean_f1 = dsearchn(f, mean(f(trim_f1)));
mean_f2 = dsearchn(f, mean(f(trim_f2)));
iaw = [mean_f1, mean_f2];

% calculate CoG for each channel spectra on basis of averaged alpha window
cogs = zeros(1,size(d0,2));
for d = 1:length(cogs)
    if isempty(trim_f1) || isempty(trim_f2)
        cogs(d) = NaN;
    else
        cogs(d) = nansum(d0(mean_f1:mean_f2,d).*f(mean_f1:mean_f2)) / sum(d0(mean_f1:mean_f2,d));
    end
end

% report which channels contribute to averaged window
sel = ~isnan(f1);

%% subfunction 4

function [selP, sums] = chanMeans(chanCogs, selG, peaks, specs, qf, cmin)
% Takes channel-wise estimates of peak alpha frequency (PAF) / centre of
% gravity (CoG) and calculates mean and standard deviation if cmin
% condition satisfied. PAFs are weighted in accordance with qf, which aims
% to quantify the relative strength of each channel peak.
%
% Part of the `restingIAF` package, (c) Andrew W. Corcoran, 2016-2017.
% Visit github.com/corcorana/restingIAF for licence, updates, and further
% information about package development and testing.
%
% Outputs:
%   selP = channels contributing peak estimates to calculation of mean PAF
%   sums = structure containing summary estimates (m, std) for PAF and CoG
%
% Inputs:
%   chanCogs = vector of channel-wise CoG estimates
%   selG = vector (logical) of channels selected for individual alpha band estimation
%   peaks = vector of channel-wise PAF estimates
%   specs = matrix of smoothed spectra (helpful for plotting weighted estimates)
%   qf = vector of peak quality estimates (area bounded by inflection points)
%   cmin = min number of channels required to generate cross-channel mean

% channel selection and weights
selP = ~isnan(peaks);               % evaluate whether channel provides estimate of PAF

% qWt = nansum(qf)/sum(selP);       % average area under peak (Qf) across viable channels (depricated: was used when calculating cross-recording comparisons)
chanWts = qf/max(qf);               % channel weightings scaled in proportion to Qf value of channel manifesting highest Qf

% average across peaks
if sum(selP) < cmin                 % if number of viable channels < cmin threshold, don't calculate cross-channel mean & std
    sums.paf = NaN;
    sums.pafStd = NaN;
    sums.muSpec = NaN;
else                                % else compute (weighted) cross-channel average PAF estimate and corresponding std of channel PAFs
    sums.paf = nansum(bsxfun(@times, peaks, chanWts))/nansum(chanWts);
    sums.pafStd = std(peaks,'omitnan');
    % estimate averaged spectra for plotting
    wtSpec = bsxfun(@times, specs, chanWts);
    sums.muSpec = nansum(wtSpec, 2)/nansum(chanWts);
end

% now for the gravs (no channel weighting, all channels included if cmin satisfied)
if sum(selG) < cmin
    sums.cog = NaN;
    sums.cogStd = NaN;
else
    sums.cog = mean(chanCogs, 2, 'omitnan');
    sums.cogStd = std(chanCogs,'omitnan');
end

%% subfunction 5

function [f1, posZ1] = findF1(f, d0, d1, negZ, minPow, slen, bin)
% Searches 1st derivative for evidence of local minima or near horizontal
% function prior to alpha peak. This location will be taken as the lower
% bound of the individual alpha band used to calculate CoG (f1).
%
% Part of the `restingIAF` package, (c) Andrew W. Corcoran, 2016-2017.
% Visit github.com/corcorana/restingIAF for licence, updates, and further
% information about package development and testing.
%
% Outputs:
%   f1 = frequency bin indexing f1
%   posZ1 = frequency of f1
%
% inputs:
%   f = frequency bin vector
%   d0 = smoothed PSD estimate vector
%   d1 = 1st derivative vector
%   negZ = vector / matrix of negative zero-crossings detected within alpha window
%   minPow = vector of minimum power threshold values defining candidate peaks
%   slen = number of bins examined for evidence of shallow slope (~1 Hz interval)
%   bin = frequency bin indexing peak / subpeak

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
    elseif abs(d1(k)) < 1 && lessThan1(d1(k+1:k+slen))
        minim = k;
        cnt = cnt+1;                % advance counter by 1
        posZ1(cnt,1) = cnt;         % zero-crossing count
        posZ1(cnt,2) = minim;       % zero-crossing frequency bin
        posZ1(cnt,3) = f(minim);    % zero-crossing frequency
    end

end

% sort out appropriate estimates for output
if size(posZ1, 1) == 1              % if singular crossing --> report frequency
    f1 = posZ1(1, 2);
    posZ1 = posZ1(1, 3);
else                                % else sort by frequency values (descending), take highest frequency (bin nearest to peak)
    posZ1 = sortrows(posZ1, -3);
    f1 = posZ1(1, 2);
    posZ1 = posZ1(1, 3);
end

%% Subfunction 6

function [f2, posZ2] = findF2(f, d0, d1, negZ, minPow, slen, bin)
% Searches 1st derivative for evidence of local minima or near horizontal
% function post alpha peak. This location will be taken as the lower
% bound of the individual alpha band used to calculate CoG (f2).
%
% Part of the `restingIAF` package, (c) Andrew W. Corcoran, 2016-2017.
% Visit github.com/corcorana/restingIAF for licence, updates, and further
% information about package development and testing.
%
% Outputs:
%   f2 = frequency bin indexing f2
%   posZ2 = frequency of f2
%
% Required inputs:
%   f = frequency bin vector
%   d0 = smoothed PSD estimate vector
%   d1 = 1st derivative of d0
%   negZ = vector / matrix of negative zero-crossings detected within alpha window
%   minPow = vector of minimum power threshold values defining candidate peaks
%   slen = number of bins examined for evidence of shallow slope (~1 Hz interval)
%   bin = frequency bin indexing peak / subpeak

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

f2 = posZ2(1, 2);
posZ2 = posZ2(1, 3);                % can simply take first estimate for output

%% subfunction 7

function tval = lessThan1(d1)
% When encounter 1st derivative absolute value < 1, eval whether following
% values (within segment) remain < +/- 1 
%
% Part of the `restingIAF` package, (c) Andrew W. Corcoran, 2016-2017.
% Visit github.com/corcorana/restingIAF for licence, updates, and further
% information about package development and testing.
%
% Output:
%   tval = logical
%
% Required input:
%   d1 = segment of 1st derivative spanning approx. 1 Hz

% setup variable check
if ~exist('d1', 'var')
    error('Provide vector containing series of 1st derivative values for slope analysis')
elseif length(d1) < 2
    error('Length of 1st derivative segment < 2');
end

t = zeros(1, length(d1));
for kx = 1:length(d1)
    t(kx) = abs(d1(kx) < 1);
end

if all(t == 1)
    tval = 1;
else
    tval = 0;
end
