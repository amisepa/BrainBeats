% VAR_DESIGN - Lagged regressors of a VAR(X) model and the regressors of
% each equation.
%
% Usage:
%   [Z, cols, P, qx] = var_design(X, U, rows, p, q)
%
% INPUTS:
%   X    - nVar x nSamples series
%   U    - nIn x nSamples exogenous inputs ([] = none)
%   rows - time points (indices) to build regressors for
%   p    - lag orders: a scalar (all equations and variables), an nVar
%          vector (one order per equation), or an nVar x nVar matrix P,
%          P(i,j) = lags of variable j in equation i (e.g. long own lags and
%          short cross lags)
%   q    - lags of U per equation (0..q): scalar or nVar vector (default:
%          the largest cross-lag order of each equation, or p)
%
% OUTPUTS:
%   Z    - (nVar*pm + nIn*(qm+1)) x numel(rows): X at lags 1..pm (block k
%          holds all variables at lag k), then U at lags 0..qm
%   cols - nVar x size(Z,1) logical: regressors of each equation
%   P, qx - the lag-order matrix and exogenous lags actually used
%
% Copyright (C) - Cedric Cannard, BrainBeats 2026

function [Z, cols, P, qx] = var_design(X, U, rows, p, q)

n = size(X,1); nu = size(U,1);
if isscalar(p)
    P = p*ones(n);
elseif isvector(p)
    P = repmat(p(:), 1, n);
else
    P = p;
end
if nargin < 5 || isempty(q)
    Pc = P; Pc(logical(eye(n))) = 0;
    qx = max(Pc, [], 2);
    if n == 1 || all(qx == 0), qx = max(P, [], 2); end
else
    qx = q(:) .* ones(n,1);
end
if nu == 0, qx = zeros(n,1); end
pm = max(P(:)); qm = max(qx)*(nu > 0);

Z = zeros(n*pm + nu*(qm+1)*(nu > 0), numel(rows));
for k = 1:pm
    Z((k-1)*n+1:k*n, :) = X(:, rows-k);
end
if nu > 0
    for k = 0:qm
        Z(n*pm + k*nu + (1:nu), :) = U(:, rows-k);
    end
end

cols = false(n, size(Z,1));
for i = 1:n
    for k = 1:pm
        cols(i, (k-1)*n + find(P(i,:) >= k)) = true;
    end
    if nu > 0
        cols(i, n*pm + (1:nu*(qx(i)+1))) = true;
    end
end
