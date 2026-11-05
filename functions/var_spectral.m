% VAR_SPECTRAL - Frequency-domain coupling measures of a VAR model:
% coherence, partial coherence, directed coherence and generalized partial
% directed coherence (Faes, Erla & Nollo, 2012).
%
% Usage:
%   sp = var_spectral(A, V, nfreq)
%
% INPUTS:
%   A, V  - VAR coefficients [A1 ... Ap] (n x n*p) and residual covariance
%   nfreq - number of frequency intervals from 0 to Nyquist
%
% OUTPUTS (sp fields; n x n x nfreq+1; element (i,j,:) is from j to i;
% squared magnitudes, 0-1):
%   w    - frequencies (rad/sample, 0 to pi)
%   S    - spectral matrix (complex)
%   coh  - coherence: total, undirected coupling
%   pcoh - partial coherence: direct, undirected coupling (given all the
%          other variables of the model)
%   dc   - directed coherence: total (direct + indirect) influence of j on
%          i, as a share of the power of i (normalized by inflow)
%   gpdc - generalized partial directed coherence: direct influence of j
%          on i, as a share of what j sends (normalized by outflow, noise
%          weighted)
%
% DC and GPDC are ratios: their values change with the other variables in
% the model, so compare them only between models with the same variables.
% Granger causality (VAR_SSGC) has no such normalization. DC and GPDC use
% the diagonal of V (strictly causal model: instantaneous effects are not
% modelled).
%
% REFERENCE:
%   Faes, L., Erla, S. & Nollo, G. (2012). Measuring connectivity in linear
%     multivariate processes: definitions, interpretation, and practical
%     analysis. Comput Math Methods Med 2012, 140513.
%
% Copyright (C) - Cedric Cannard, BrainBeats 2026

function sp = var_spectral(A, V, nfreq)

n = size(A,1); p = size(A,2)/n;
w = pi*(0:nfreq)/nfreq; nf = numel(w);

% Abar(w) = I - sum_k A_k e^{-ikw}, all frequencies at once
E = exp(-1i*(1:p)'*w);                                   % p x nf
Abar = repmat(eye(n),1,1,nf) - reshape(reshape(A, n*n, p)*E, n, n, nf);
H = pagemldivide(Abar, repmat(eye(n),1,1,nf));           % transfer function (entries needed)
S = pagemtimes(pagemtimes(H, V), 'none', H, 'ctranspose');
Pinv = pagemtimes(pagemtimes(Abar, 'ctranspose', V \ eye(n), 'none'), Abar);   % inverse spectral matrix

dS = real(diag_pages(S)); dP = real(diag_pages(Pinv));     % n x nf
coh  = abs(S).^2    ./ (permute(dS,[1 3 2]) .* permute(dS,[3 1 2]));
pcoh = abs(Pinv).^2 ./ (permute(dP,[1 3 2]) .* permute(dP,[3 1 2]));

sig2 = diag(V)';                                           % 1 x n
num = abs(H).^2 .* sig2;                                   % |H_ij|^2 sigma_j^2
dc = num ./ sum(num, 2);                                   % normalized by row (inflow to i)
num = abs(Abar).^2 ./ sig2';                               % |Abar_ij|^2 / sigma_i^2
gpdc = num ./ sum(num, 1);                                 % normalized by column (outflow of j)

sp = struct('w', w, 'S', S, 'coh', coh, 'pcoh', pcoh, 'dc', dc, 'gpdc', gpdc);


function d = diag_pages(M)
[n, ~, nf] = size(M);
d = reshape(M(repmat(logical(eye(n)), 1, 1, nf)), n, nf);
