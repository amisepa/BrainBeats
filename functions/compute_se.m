%% Computes sample entropy (SE) of a univariate signal. 
% Sample entropy quantifies the likelihood that a sequence of m consecutive 
% data points that matches another sequence of the same length (match within
% a tolerance of r) will still match the other sequence when their length 
% is increased of one sample (sequences of length m + 1).
% 
% Inputs:
%   x   - univariate signal - a vector of size 1 x n (the number of sample points)
%   m   - embedding dimension (default = 2)
%   r   - threshold/tolerance (.15 of the signal' sd)
%   tau - time lag (it is usually equal to 1)
%
% Outputs:
%   se  - sample entropy
%   p   - a vector of length 2 (the number of template matches of length m 
%           and number of forward matches of length m+1).
%
% Please cite:
%   Azami & Escudero, "Refined Multiscale Fuzzy Entropy based on Standard Deviation for Biomedical Signal Analysis", 
%       Medical & Biological Engineering & Computing, 2016.
%   Richman & Moorman, "Physiological time-series analysis using approximate entropy and sample entropy"
%       American Journal of Physiology-Heart and Circulatory Physiology, vol. 278, no. 6, pp.H2039-H2049, 2000.
%
% Cedric Cannard, August 2022

function [entropy,p] = compute_se(signal,m,r,tau)

% Defaults
if ~exist('m', 'var'), m = 2; end
if ~exist('r', 'var'), r = .15*std(signal); end
if ~exist('tau', 'var'), tau = 1; end

% Downsample
if tau > 1, signal = downsamp(signal, tau); end

n = length(signal);
p = zeros(1,2);
sMat = zeros(m+1,n-m);
parfor i = 1:m+1
    sMat(i,:) = signal(i:n-m+i-1);
end

for k = m:m+1
    count = zeros(1,n-m);
    tempMat = sMat(1:k,:);
    
    parfor i = 1:n-k
        % calculate Chebyshev distance without counting self-matches
        dist = max(abs(tempMat(:,i+1:n-m) - repmat(tempMat(:,i),1,n-m-i)));
        
        % calculate the number of distance values that are less than the threshold r
        D = (dist < r);
        count(i) = sum(D)/(n-m);
    end
    
    p(k-m+1) = sum(count)/(n-m);
end
entropy = log(p(1)/p(2));


% Downsample subfunction
function y = downsamp(x, n, phase)

if nargin<2 || nargin>3, print_usage; end

if nargin<3
  phase = 0;
end

if phase > n - 1
  warning('This is incompatible with Matlab (phase = 0:n-1). See octave-forge signal package release notes for details.')
end

if isvector(x)
  y = x(phase + 1:n:end);
else
  y = x(phase + 1:n:end,:);
end


