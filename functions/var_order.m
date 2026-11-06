% VAR_ORDER - VAR (or VARX) model order by the Bayesian (BIC) and Akaike
% (AIC) information criteria: for the whole model, and per equation with
% separate own-lag and cross-lag orders.
%
% Usage:
%   [pBIC, pAIC, crit] = var_order(X, pmax)
%   [pBIC, pAIC, crit] = var_order(X, pmax, U, pown, maxcoef)
%
% INPUTS:
%   X       - nVar x nSamples series, demeaned; NaN = missing (see VAR_FIT)
%   pmax    - largest order tried (whole model; cross lags per equation)
%   U       - optional exogenous inputs, with lags 0..p at order p
%   pown    - largest own-lag order tried per equation (default pmax)
%   maxcoef - largest number of coefficients per equation (default Inf)
%
% OUTPUTS:
%   pBIC, pAIC - orders minimizing each criterion (whole model)
%   crit       - struct: p (1:pmax), bic, aic (whole model, per order);
%                P (nVar x nVar, BIC per equation: P(i,i) own lags, P(i,j)
%                cross lags) and qx (lags of U per equation, 0..pmax,
%                chosen separately: respiration can act with longer lags
%                than the cross terms). Series with different dynamics (EEG power
%                envelopes, the interpolated NN series) need different
%                orders; a smooth series needs long own lags but no long
%                cross lags (Hsiao, 1981).
%
% Every order is fitted on the same time points (those valid at the
% largest order), so the criteria compare models of the same data.
%
% Copyright (C) - Cedric Cannard, BrainBeats 2026

function [pBIC, pAIC, crit] = var_order(X, pmax, U, pown, maxcoef)

if nargin < 3, U = []; end
if nargin < 4 || isempty(pown), pown = pmax; end
if nargin < 5 || isempty(maxcoef), maxcoef = Inf; end
n = size(X,1); nu = size(U,1);
pl = max(pmax, pown);
[~, ~, info] = var_fit(X, pl, [], U, pl);
rows = info.rows;
nObs = info.nObs;

% whole model, same order for everything
bic = nan(1,pmax); aic = nan(1,pmax);
for p = 1:pmax
    [~, V] = var_fit(X, p, rows, U, p);
    ld = 2*sum(log(diag(chol(V))));
    k = n*(n*p + nu*(p+1));
    bic(p) = nObs*ld + k*log(nObs);
    aic(p) = nObs*ld + 2*k;
end
[~, pBIC] = min(bic);
[~, pAIC] = min(aic);

% per equation: own lags x cross lags x input lags (0..pmax)
Z = var_design(X, U, rows, pl, pl);
P = zeros(n); qx = zeros(n,1);
qlist = 0; if nu > 0, qlist = 0:pmax; end
for i = 1:n
    y = X(i, rows)';
    best = Inf;
    for po = 1:pown
        for pc = 1:pmax
            Pi = pc*ones(1,n); Pi(i) = po;
            mx = false(1, size(Z,1));
            for k = 1:pl, mx((k-1)*n + find(Pi >= k)) = true; end
            for pq = qlist
                m = mx;
                if nu > 0, m(n*pl + (1:nu*(pq+1))) = true; end
                kc = sum(m);
                if kc > maxcoef, continue; end
                Zi = Z(m,:)';
                e = y - Zi*(Zi \ y);
                b = nObs*log(e'*e/nObs) + kc*log(nObs);
                if b < best, best = b; P(i,:) = Pi; qx(i) = pq; end
            end
        end
    end
end
crit = struct('p', 1:pmax, 'bic', bic, 'aic', aic, 'P', P, 'qx', qx);
