% VAR_SIM - Simulate a VAR (or VARX) model with Gaussian innovations.
%
% Usage:
%   X = var_sim(A, V, N, nTrials, burn)
%   X = var_sim(A, V, N, nTrials, burn, B, U)
%
% INPUTS:
%   A, V    - VAR coefficients [A1 ... Ap] (n x n*p) and innovation
%             covariance (n x n)
%   N       - samples per trial
%   nTrials - independent realizations (default 1)
%   burn    - samples discarded at the start of each (default 1000)
%   B, U    - optional exogenous part: B = [B0 ... Bq] (n x nIn*(q+1)) and
%             the input U (nIn x N, the same for every trial; NaN read as
%             0). The burn-in runs with zero input.
%
% OUTPUT:
%   X - n x N x nTrials
%
% Draw the random stream with rng() before calling for reproducibility.
%
% Copyright (C) - Cedric Cannard, BrainBeats 2026

function X = var_sim(A, V, N, nTrials, burn, B, U)

if nargin < 4 || isempty(nTrials), nTrials = 1; end
if nargin < 5 || isempty(burn), burn = 1000; end
if nargin < 7, B = []; U = []; end
n = size(A,1); p = size(A,2)/n;
T = N + burn;
L = chol(V, 'lower');
X = reshape(L*randn(n, T*nTrials), n, T, nTrials);   % innovations, overwritten in place
if ~isempty(B)
    U(~isfinite(U)) = 0;
    nu = size(U,1); q = size(B,2)/nu - 1;
    Ub = [zeros(nu, burn+q) U];                       % input with zero burn-in
    D = zeros(n, T);                                  % exogenous drive, same for all trials
    for k = 0:q
        D = D + B(:, k*nu+(1:nu)) * Ub(:, (q+1:q+T) - k);
    end
    X = X + D;
end
X = permute(X, [1 3 2]);                              % n x nTrials x T: one page per time
Ak = reshape(A, n, n, p);
for t = p+1:T
    acc = X(:,:,t);
    for k = 1:p
        acc = acc + Ak(:,:,k)*X(:,:,t-k);
    end
    X(:,:,t) = acc;
end
X = permute(X(:,:,burn+1:end), [1 3 2]);
