% Calculates the refined composite multiscale fuzzy entropy (RCMFE) 
% whose coarse-graining uses standard deviation.
% 
% INPUTS:
%   signal: univariate signal - a vector of size 1 x N (the number of sample points)
%   m: embedding dimension
%   r: threshold (it is usually equal to 0.15 of the standard deviation of a signal - 
%       because we normalize signals to have a standard deviation of 1, here, r is usually equal to 0.15)
%   tau: time lag (it is usually equal to 1)
%   coarseType: 'Mean', 'Standard deviation' (default), 'Variance'
%   nScales: the number of scale factors (default = 15)
%   fs: sample rate (Hz)
%   n: fuzzy power
%   usegpu: use GPU computing (1), or not (0)
%
% Outputs:
%   rcmfe: entropy values for each scale factor
%   scales: lower and upper frequency bounds of each time scale
%
% 
% Ref:
%   [1] H. Azami and J. Escudero, "Refined Multiscale Fuzzy Entropy based on
%   Standard Deviation for Biomedical Signal Analysis", Medical & Biological
%   Engineering & Computing, 2016.
%
% Cedric Cannard, 2022

function [rcmfe, scales] = compute_rcmfe(signal,m,r,tau,coarseType,nScales,fs,n,useGPU)

% if not srate provided somehow, try to interpolate
% if ~exist('fs', 'var')
% 	fs = (1/(EEG.times(end)) - EEG.times(1));
% 	warning(fprintf('No sample rate inputted, sample rate estimated: %g', fs))
%     errordlg('You need to input the sample rate.'); return
% end

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

% Center and normalize signal to SD of 1
signal = zscore(signal);

% Simplify SD coarse-graining name
if contains(lower(coarseType), 'standard')
    coarseType = 'SD';
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

rcmfe = nan(1,nScales);
for iScale = 1:nScales
    
    fprintf('Scale %g \n', iScale);

    temp_A = [];
    temp_B = [];
    
    parfor ii = 1:iScale
        
        if useGPU
            sig = gpuArray(signal(ii:end));
        else
            sig = signal(ii:end);
        end
    
        % Coarse-graining process 
        y = reshape(sig(1:floor(length(sig)/iScale)*iScale),iScale,[]);
        switch coarseType
            case 'Mean'
                sigCoarse = mean(y,'omitnan');
            case 'SD'
                sigCoarse = std(y,'omitnan');
            case 'Variance'
                sigCoarse = var(y,'omitnan');
        end
        
        % Compute fuzzy entropy on the coearsed signal
        [~, p] = compute_fe(sigCoarse,m,r,n,tau,useGPU);

        temp_A = [temp_A p(1)];
        temp_B = [temp_B p(2)];
        
    end
    
    % output
    A = sum(temp_A);
    B = sum(temp_B);
    rcmfe(iScale) = log(A/B);
    
end

% Scale factors to export
scales = 1:nScales;

% Remove NaN scales
scales(isnan(rcmfe)) = [];
rcmfe(isnan(rcmfe)) = [];

