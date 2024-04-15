%% Computes multiscale fuzzy entropy (MFE)
%
% INPUTS:
%   signal: univariate signal - a vector of size 1 x N (the number of sample points)
%   m: embedding dimension (default = 2)
%   r: threshold (it is usually equal to 0.15 of the standard deviation of a signal - 
%       because we normalize signals to have a standard deviation of 1, here, r is usually equal to 0.15)
%   tau: time lag (it is usually equal to 1)
%   coarseType: 'Mean', 'Standard deviation' (default), 'Variance'
%   nScales: the number of scale factors (default = 15)
%   filtData: bandpass filter each scale factor to control for spectral 
%       bias (1) or not (0; default)
%   fs: sample rate (Hz)
%   n: fuzzy power (default = 2)
%
% OUTPUTS:
%   mfe: entropy values for each scale factor
%   scales: scale factors
%   scales_bounds: lower and upper frequency bounds of each time scale
% 
% EXAMPLE:
%   [mfe, scales, scale_bounds] = compute_mfe(signal, m, r, tau, coarseType, nScales, filtData, fs, n, useGPU)
%   [mfe, scales, scale_bounds] = compute_mfe(EEG.data(iChan,:), [], [], [], 'SD', 30, 0, EEG.srate, [], 0)
% 
% Please cite:
%   [1] Azami & Escudero (2016), "Refined Multiscale Fuzzy Entropy based on
%   Standard Deviation for Biomedical Signal Analysis", Medical & Biological
%   Engineering & Computing.
% 
%   [2] Costa, Goldberger, Peng (2002). Multiscale entropy analysis of 
%   complex physiologic time series. Phys Rev Lett. 
% 
%   [3] Kosciessa, Kloosterman, Garrett (2020). Standard multiscale entropy 
%   reflects neural dynamics at mismatched temporal scales: What's signal 
%   irregularity got to do with it? Plos Comput Biol.
% 
%   [4] Grandy, Garrett, Schmiedek, Werkle-Bergner (2016). On the estimation 
%   of brain signal entropy from sparse neuroimaging data. Sci Rep.
%
% Cedric Cannard, 2022

function [mfe, scales, scale_bounds] = compute_mfe(signal, m, r, tau, coarseType, nScales, filtData, fs, n, useGPU)

% FIXME: Determine max number of scales using file length. As a rough guideline, 
% some researchers suggest having at least 10 times as many data points as 
% the embedding dimension m used in the entropy calculation. For instance, 
% if you are using an embedding dimension of 2, you might want to have at 
% least 20 data points in each coarse-grained time series. So, the maximum 
% scale factor τ_max that you could use would be approximately N/20 when 
% using an embedding dimension of 2.

% Lowest_freq = 1 / (length(signal)/1000)
% highest_freq = fs / (2*nScales)
scale_bounds = {};

if ~exist('m','var') || isempty(m)
    m = 2; 
end
if ~exist('r','var')|| isempty(r)
    r = .15; 
end
if ~exist('tau','var') || isempty(tau) 
    tau = 1; 
end
if ~exist('coarseType','var') || isempty(coarseType)
    coarseType = 'SD'; 
end
if ~exist('nScales','var') || isempty(nScales)
    nScales = 20; 
end
if ~exist('filtData','var') || isempty(filtData)
    filtData = 0; 
end
if ~exist('n','var') || isempty(n)
    n = 2; 
end
if ~exist('useGPU','var') || isempty(useGPU)
    useGPU = 0; 
end

% Downsample signal if tau is > 1
if tau > 1
    signal = downsample(signal, tau); 
end

% Max scale factor cannot be greater than Nyquist frequency (for EEG)
if exist('fs','var')
    nf = fs/2;
    if nScales >= nf
        warning("Scale factor cannot be as high as the Nyquist frequency. Lowering it to %g", nf-1);
    end
end

% Center and normalize signal to SD of 1
signal = zscore(signal);

% Simplify SD coarse-graining name
if contains(lower(coarseType), 'standard')
    coarseType = 'SD';
end

if filtData
    disp('Applying bandpasss-filter to remove spectral bias at each scale (please reference Kosciessa et al. (2020).')
end

mfe = nan(1,nScales);
parfor iScale = 1:nScales

    fprintf('   scale %d \n', iScale)
  
    % Make copy of signal in case it is bandpass-filtered at each scale
    sig = signal;

    % Move it to GPU if applicable
    if useGPU && length(sig) > 1000
        try
            sig = gpuArray(sig);
            disp('Using GPU computing')
        catch
            disp('Could not use GPU computing.')
        end
    end

    % Corresponding scale frequency bounds
    if exist('fs','var')
        if iScale == 1
            upperBound = nf;
        else
            upperBound = (1/iScale).*nf + .05*((1./iScale).*nf);
        end
        lowerBound = (1/(iScale+1)).*nf - .05*((1./(iScale+1)).*nf);
        
        % scale_bounds(iScale,:) = {round(lowerBound,3) round(upperBound,3) };
        scale_bounds(iScale,:) = {sprintf('scale %g: %g - %g',iScale, round(lowerBound,1), round(upperBound,1))} ;
    end

    % Control for broadband spectral contributions using a bandpass filter at each scale factor (see Kosciessa et al 2020)
    if filtData
		
		% FIXME: can filters be designed outside loops to increase speed?
        if iScale == 1
            [b,a] = butter(10,lowerBound/nf,'high');   % use Butterworth highpass
            sig = filtfilt(b,a,sig);
       else
            if upperBound-lowerBound > .05*nf    % for broad passbands: Chebyshev Type I zero-phase bandpass
                [b,a] = cheby1(10,1,upperBound/nf,'low');   % lowpass
                % figure; freqz(b,a);
                sig = filtfilt(b,a,sig);
                [b,a] = cheby1(10,1,lowerBound/nf,'high');  % highpass
                sig = filtfilt(b,a,sig);
            else                                 % Butterworth zero-phase bandpass
                [b,a] = butter(4,upperBound/nf,'low');
                sig = filtfilt(b,a,sig);
                [b,a] = butter(4,lowerBound/nf,'high');
                sig = filtfilt(b,a,sig);
            end
        end
        
%         % Visualize filter effect on the power spectrum
%         [psd,f] = pwelch(sig,[],[],[],fs);
%         plot(f,psd); hold on;
%         title(num2str(iScale)); legend([num2str(lowcutoff) '-' num2str(highcutoff)]);        
    end

    % Because of filtering, the scale-wise SD decreases relative to the global scale-
    % invariant similarity bound r [29]. Thus, r is recalculated for each scale,
    % thereby normalizing MSE with respect to changes in overall time series variation at each scale
    if filtData
        % r_adj = 0.5*std(sig);
        r_adj = r*std(sig);
    else
        r_adj = r;
    end

    y = reshape(sig(1:floor(length(sig)/iScale)*iScale),iScale,[]);

    switch coarseType
        case 'Mean'
            sig = mean(y,'omitnan');
        case 'SD'
            sig = std(y,'omitnan');
        case 'Variance'
            sig = var(y,'omitnan');
    end

    mfe(:,iScale) = compute_fe(sig, m, r_adj, n, tau, useGPU);

end

% scale factors
scales = 1:nScales;

% Remove NaN scales
scales(isnan(mfe)) = [];
mfe(isnan(mfe)) = [];
if exist('fs','var')
    scale_bounds(isnan(mfe),:) = [];
end

% Plot
figure; plot(scales,mfe)
