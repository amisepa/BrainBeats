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

function [rcmfe, scales] = compute_rcmfe(signal,m,r,tau,coarseType,nScales,fs,n,usegpu)

% if not srate provided somehow, try to interpolate
% if ~exist('fs', 'var')
% 	fs = (1/(EEG.times(end)) - EEG.times(1));
% 	warning(fprintf('No sample rate inputted, sample rate estimated: %g', fs))
%     errordlg('You need to input the sample rate.'); return
% end

% Simplify SD coarse-graining name
if contains(lower(coarseType), 'standard')
    coarseType = 'SD';
end

% max scale factor
nf = fs/2;  % Nyquist frequency
if nScales >= nf
	warning("Scale factor cannot be as high as signal's Nyquist frequency. Lowering it to %g", nf-1);
end

% Signal is centered and normalized to standard deviation 1
disp('Normalizing signal to SD = 1.')
signal = signal - mean(signal);
signal = signal ./ std(signal);

rcmfe = nan(1,nScales);
for iScale = 1:nScales
    
    fprintf('Scale %g \n', iScale);

    temp_A = [];
    temp_B = [];
    
    for ii = 1:iScale
        
        if usegpu
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
        [~, p] = compute_fe(sigCoarse,m,r,n,tau,usegpu);

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

