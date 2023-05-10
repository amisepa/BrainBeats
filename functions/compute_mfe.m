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
%   scales: lower and upper frequency bounds of each time scale
% 
% Ref:
%   [1] Azami & Escudero (2016), "Refined Multiscale Fuzzy Entropy based on
%   Standard Deviation for Biomedical Signal Analysis", Medical & Biological
%   Engineering & Computing.
% 
%   [2] Costa, Goldberger, Peng (2002). Multiscale entropy analysis of 
%   complex physiologic time series. Phys Rev Lett. 
%
% Cedric Cannard, 2022

function [mfe, scales] = compute_mfe(signal, m, r, tau, coarseType, nScales, filtData, fs, n, useGPU)

% FIXME: Determine max number of scales using file length. As a rough guideline, 
% some researchers suggest having at least 10 times as many data points as 
% the embedding dimension m used in the entropy calculation. For instance, 
% if you are using an embedding dimension of 2, you might want to have at 
% least 20 data points in each coarse-grained time series. So, the maximum 
% scale factor τ_max that you could use would be approximately N/20 when 
% using an embedding dimension of 2.

% Lowest_freq = 1 / (length(signal)/1000)
% highest_freq = fs / (2*nScales)

if ~exist('m','var'), m = 2; end
if ~exist('r','var'), r = .15; end
if ~exist('n','var'), n = 2; end
if ~exist('tau','var'), tau = 1; end
if tau > 1, signal = downsample(signal, tau); end

% Max scale factor cannot be greater than Nyquist frequency (for EEG)
if exist('fs','var')
    nf = fs/2;
    if nScales >= nf
        warning("Scale factor cannot be as high as the Nyquist frequency. Lowering it to %g", nf-1);
    end
end

% Move it to GPU if applicable
if useGPU && length(signal) > 1000
    try
        signal = gpuArray(signal);
        disp('Using GPU computing')
    catch
        disp('Could not use GPU computing.')
    end
end

% Signal is centered and normalised to standard deviation 1
signal = zscore(signal);

% Simplify SD coarse-graining name
if contains(lower(coarseType), 'standard')
    coarseType = 'SD';
end

mfe = nan(1,nScales);
scales = nan(1,nScales);
for iScale = 1:nScales

    fprintf('   scale %d \n', iScale)
  
    % Make copy of signal in case it is bandpass-filtered at each scale
    sig = signal;


    % Bandpass filter outside these bounds to control for spectral bias (see Kosciessa et al 2020)
    if filtData
        disp('Applying bandpasss-filter to remove spectral bias (please reference Kosciessa et al. (2020).')

        % Scale frequency bounds
        upperBound = (1/iScale).*nf + .05*((1./iScale).*nf);
        lowerBound = (1/(iScale+1)).*nf - .05*((1./(iScale+1)).*nf);
        scales(:,iScale) = [round(lowerBound,3) round(upperBound,3) ];


        if iScale == 1
            [D,C] = butter(10,lowerBound/nf,'high');   % use Butterworth highpass
            upperBound = nf;
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
        
        % FIXME: add this from filedtrip code and move filtering below this
        % padlength = ceil(size(data.trial{1},2)./2); % use half the length of trial 1 as padding (JQK)
        % x_pad = cellfun(@(a) ft_preproc_padding(a, 'mean', padlength), data.trial, 'UniformOutput', false );    % add padding
        % x_pad = cellfun(@transpose, x_pad, 'UniformOutput', false);                                                 % transpose for filtfilt: time x chan
        % if sc == 1 % only HPF
        %     resamp_x_pad = cellfun(@(x_pad) filtfilt(D,C,x_pad), x_pad, 'UniformOutput', false );  % high-pass filter data
        % else
        %     resamp_x_pad = cellfun(@(x_pad) filtfilt(B,A,x_pad), x_pad, 'UniformOutput', false );                       % low-pass filter data
        %     resamp_x_pad = cellfun(@(resamp_x_pad) filtfilt(D,C,resamp_x_pad), resamp_x_pad, 'UniformOutput', false );  % high-pass filter data
        % end
        % resamp_x_pad = cellfun(@transpose, resamp_x_pad, 'UniformOutput', false);                                   % transpose back : chan x time again
        % resamp_x = cellfun(@(resamp_x_pad) ft_preproc_padding(resamp_x_pad, 'remove', padlength), ...                % remove padding
        %     resamp_x_pad, 'UniformOutput', false );
        % %figure; hold on; plot(resamp_x{1}(1,:)); plot(data.trial{1}(1,:))
        % % create data_filt structure
        % data_filt = data;
        % data_filt.trial = resamp_x;
        % clear resamp_* x_pad;

%         % Visualize filter effect on the power spectrum
%         [psd,f] = pwelch(sig,[],[],[],fs);
%         plot(f,psd); hold on;
%         title(num2str(iScale)); legend([num2str(lowcutoff) '-' num2str(highcutoff)]);        
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

    mfe(:,iScale) = compute_fe(sig, m, r, n, tau,useGPU);

end

% scale numbers instead of freqs for NN series
if ~filtData
    scales = 1:nScales;
end

% Remove NaN scales
idx = isnan(mfe);    
mfe(idx) = [];
scales(idx) = [];

