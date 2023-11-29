%% Fast method to compute sample entropy.
% Sample entropy quantifies the likelihood that a sequence of m consecutive 
% data points that matches another sequence of the same length (match within
% a tolerance of r) will still match the other sequence when their length 
% is increased of one sample (sequences of length m + 1).
% 
% Developed by Shamim Nemati and Giulia Da Poian (Physionet toolbox)
% Please cite: 
%   Vest A, Da Poian G, Li Q, Liu C, Nemati S, Shah A, Clifford GD, 
%   An Open Source Benchmarked Toolbox for Cardiovascular Waveform and Interval Analysis", 
%   Physiological measurement 39, no. 10 (2018): 105004.
% 
% Inputs:
%   y   - num_var x num_samples
%   m   - template length
%   r   - radius
%
% Output:
%   se  - sample Entropy
%
% Cedric Cannard, 2022

function se = compute_se_fast(signal,m,r)

if size(signal, 1) > size(signal,2), signal = signal'; end

xx = convert_to_lagged_form(signal, m)';
Dxx = pdist(xx,'chebychev');

yy = convert_to_lagged_form(signal, m+1)';
Dyy = pdist(yy,'chebychev');

A = mean( Dxx < r ) ;
B = mean( Dyy < r );

se = -log(B/A);

%% Subfunctions

% Create an observation vector yy(:,t) containing the last k values of y, newest first
% e.g., k=2, y = (a1 a2 a3)     yy  = a2 a3
%                (b1 b2 b3)           b2 b2
%                                     a1 a2
%                                     b1 b2
function yy = convert_to_lagged_form(y, k)
[s, T] = size(y);
bs = s*ones(1,k);
yy = zeros(k*s, T-k+1);
for i = 1:k
    yy(block(i,bs), :) = y(:, k-i+1:end-i+1); 
end

% BLOCK Return a vector of subscripts corresponding to the specified blocks.
% sub = block(blocks, block_sizes)
% e.g., block([2 5], [2 1 2 1 2]) = [3 7 8].
function sub = block(blocks, block_sizes)
blocks = blocks(:)';
block_sizes = block_sizes(:)';
skip = [0 cumsum(block_sizes)];
start = skip(blocks)+1;
fin = start + block_sizes(blocks) - 1;
sub = [];
for j = 1:length(blocks)
    sub = [sub start(j):fin(j)];
end
