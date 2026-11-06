% DATA_RANK - Rank of multichannel data: the number of eigenvalues of the
% channel covariance above a tolerance relative to the largest one.
%
% Tolerance: (n * eps('single'))^2 * the largest eigenvalue, n = number of
% channels (the singular-value rule n * eps * sigma_max, squared for
% variances; single precision because EEGLAB stores data in single). A
% direction removed by an average reference, an interpolation or a removed
% component keeps a variance of the order of eps('single')^2 relative to the
% channels', far below the tolerance, and real EEG directions stay far above
% it. Unlike the former absolute threshold (1e-7, in the data's units^2),
% it does not depend on the units (uV, V).
%
% Usage:
%   r = data_rank(X)
%
% Input:
%   X - data, channels x samples (or channels x samples x epochs)
%
% Output:
%   r - rank
%
% Copyright (C) - Cedric Cannard, 2026

function r = data_rank(X)
X = double(X(:,:));
X = X - mean(X, 2);
C = (X * X') / max(1, size(X,2) - 1);
e = eig((C + C') / 2);
r = sum(e > max(e) * (size(X,1) * eps('single'))^2);
end
