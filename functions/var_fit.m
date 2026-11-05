% VAR_FIT - Least-squares fit of a vector autoregressive model, with
% optional exogenous inputs (VARX), to series with gaps.
%
%   X_t = A1*X_{t-1} + ... + Ap*X_{t-p} + B0*U_t + ... + Bq*U_{t-q} + e_t,
%   cov(e) = V
%
% Usage:
%   [A, V, info] = var_fit(X, p)
%   [A, V, info] = var_fit(X, p, rows, U, q)
%
% INPUTS:
%   X    - nVar x nSamples series, demeaned (no intercept is fitted). NaN
%          marks missing samples: a time point is used only if it and its
%          lags (of X and U) are all finite, so gaps and segment edges
%          never mix.
%   p    - model order (number of lags of X); or one order per equation
%          (nVar vector), or an nVar x nVar matrix P of the lags of
%          variable j in equation i (Hsiao, 1981; see VAR_DESIGN): A and B
%          are then zero beyond each equation's lags
%   rows - optional: time points to use (indices); default all valid ones.
%          VAR_ORDER passes the same points to every order so that the
%          information criteria compare like with like.
%   U    - optional nIn x nSamples exogenous inputs (e.g. respiration): they
%          drive X but are not modelled themselves (Baselli et al., 1997)
%   q    - lags of U (0..q), scalar or per equation (default: each
%          equation's largest cross-lag order)
%
% OUTPUTS:
%   A    - nVar x nVar*p coefficients [A1 ... Ap]
%   V    - nVar x nVar residual covariance (maximum likelihood: E*E'/nObs)
%   info - struct: B (nVar x nIn*(q+1), [B0 ... Bq]), nObs, rows, E
%          (residuals), rho (spectral radius of the companion matrix; < 1 =
%          stable), loglik (Gaussian, up to a constant)
%
% The regression is solved by QR (backslash), never through the normal
% equations. Granger causality between components of X given an
% exogenous U follows from A and V alone (VAR_SSGC): the input is known in
% both the full and the reduced models.
%
% Copyright (C) - Cedric Cannard, BrainBeats 2026

function [A, V, info] = var_fit(X, p, rows, U, q)

if nargin < 4, U = []; end
if nargin < 5, q = []; end
n = size(X,1);
nu = size(U,1);
[~, ~, P, qx] = var_design(X, U, [], p, q);
pm = max(P(:)); qm = max(qx);
if nargin < 3 || isempty(rows)
    rows = var_rows([X; U], max(pm, qm));
end
rows = rows(:)';
nObs = numel(rows);
[Z, cols] = var_design(X, U, rows, P, qx);
k = max(sum(cols, 2));
if nObs <= k
    error('var_fit: %d usable time points for %d coefficients per equation (order %d).', nObs, k, pm)
end

Y = X(:, rows);
if all(cols(:))
    C = (Z' \ Y')';       % QR least squares, all equations at once
else
    % equation-specific lags (Hsiao, 1981): each equation on its own
    % regressors, zero elsewhere
    C = zeros(n, size(Z,1));
    for i = 1:n
        C(i, cols(i,:)) = (Z(cols(i,:),:)' \ Y(i,:)')';
    end
end
A = C(:, 1:n*pm);
E = Y - C*Z;
V = (E*E') / nObs;
V = (V + V')/2;

info.B = C(:, n*pm+1:end);
info.q = qm;
info.P = P;
info.qx = qx;
info.order = max(P, [], 2)';
info.ncoef = sum(cols, 2)';
info.nObs = nObs;
info.rows = rows;
info.E = E;
info.rho = max(abs(eig(companion(A))));
R = chol(V);
info.loglik = -nObs/2 * 2*sum(log(diag(R))) - nObs*n/2;


function rows = var_rows(X, p)
% Time points whose value and p lags are all finite
ok = all(isfinite(X), 1);
N = numel(ok);
c = cumsum([0 ~ok]);                          % number of bad samples up to each point
t = p+1:N;
rows = t(c(t+1) - c(t-p) == 0);               % no bad sample in t-p..t


function Ac = companion(A)
n = size(A,1); p = size(A,2)/n;
Ac = [A; eye(n*(p-1)) zeros(n*(p-1), n)];
