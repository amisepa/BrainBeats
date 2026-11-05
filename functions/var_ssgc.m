% VAR_SSGC - Granger causality of a VAR model, in the time and frequency
% domains, through its state-space form (Barnett & Seth, 2015).
%
% The VAR is written in innovations form; the model of any subset of its
% variables (without some of them) then follows exactly from a discrete
% algebraic Riccati equation instead of a second, truncated regression.
% This avoids the bias and the negative values of spectral Granger
% causality computed from two separately fitted models (Stokes & Purdon,
% 2017), and the spectral values average to the time-domain value.
%
% Usage:
%   [F, f, chk] = var_ssgc(A, V, tgt, src, cond, nfreq)
%
% INPUTS:
%   A, V  - VAR coefficients [A1 ... Ap] (n x n*p) and residual covariance
%           (see VAR_FIT)
%   tgt   - index of the target variable(s)
%   src   - index of the source variable(s)
%   cond  - indices of the variables conditioned on ([] = none). Variables
%           in none of tgt, src, cond are left out of the model: with
%           cond = [], F is the pairwise Granger causality of the [tgt src]
%           process.
%   nfreq - number of frequency intervals from 0 to Nyquist (f has nfreq+1
%           values, at pi*(0:nfreq)/nfreq rad/sample); 0 or omitted = time
%           domain only
%
% OUTPUTS:
%   F   - Granger causality from src to tgt given cond (time domain;
%         natural log of a ratio of prediction error variances; >= 0)
%   f   - spectral Granger causality, 1 x nfreq+1; its mean over the
%         frequencies (trapezoid) equals F
%   chk - struct of numerical checks: dare (largest relative Riccati
%         residual), flat (conditional case: largest relative deviation of
%         the spectrum of the reduced innovations from flat; should be ~0)
%
% Conditional spectral GC: the [tgt cond] part of the process is filtered
% by the inverse of its own (reduced) model, and the spectral GC from
% [cond src] to the filtered target is computed as in the unconditional
% case (Geweke, 1984; as in MVGC, Barnett & Seth 2014).
%
% REFERENCES:
%   Barnett, L. & Seth, A.K. (2015). Granger causality for state-space
%     models. Phys Rev E 91, 040101.
%   Geweke, J. (1984). Measures of conditional linear dependence and
%     feedback between time series. J Am Stat Assoc 79, 907-915.
%   Chu, E.K. et al. (2005). Structure-preserving algorithms for periodic
%     discrete-time algebraic Riccati equations (SDA). Int J Control 77.
%
% Copyright (C) - Cedric Cannard, BrainBeats 2026

function [F, f, chk] = var_ssgc(A, V, tgt, src, cond, nfreq)

if nargin < 5, cond = []; end
if nargin < 6, nfreq = 0; end
n = size(A,1); p = size(A,2)/n; m = n*p;
Ac = [A; eye(m-n) zeros(m-n, n)];      % companion (state transition)
C  = A;                                % observation: X_t = C z_t + e_t
K  = [eye(n); zeros(m-n, n)];          % state input: z_{t+1} = Ac z_t + K e_t

tgt = tgt(:)'; src = src(:)'; cond = cond(:)';
s = [tgt src cond];
nt = numel(tgt); ny = numel(src); nz = numel(cond);
it = 1:nt; iy = nt+(1:ny); iz = nt+ny+(1:nz);   % positions within s
chk.dare = 0;

% Innovations model of the [tgt src cond] process
Cs = C(s,:);
if numel(s) == n
    Ks = K(:,s); Vs = V(s,s);          % all variables: reorder only
else
    [Ks, Vs, r1] = ss_reduce(Ac, C, K, V, s);
    chk.dare = max(chk.dare, r1);
end

% Reduced model without the source
ir = [it iz];
[Kr, Vr, r2] = ss_reduce(Ac, Cs, Ks, Vs, ir);
chk.dare = max(chk.dare, r2);
F = logdet(Vr(it,it)) - logdet(Vs(it,it));
chk.flat = NaN;

f = [];
if nfreq == 0, return; end

w = pi*(0:nfreq)/nfreq; nf = numel(w);
ZA = repmat(-Ac, 1, 1, nf);            % e^{iw} I - Ac at each frequency
for k = 1:nf
    ZA(:,:,k) = ZA(:,:,k) + exp(1i*w(k))*eye(m);
end
Hs = eye(nt+ny+nz) + pagemtimes(Cs, pagemldivide(ZA, repmat(Ks,1,1,nf)));   % transfer function of s
f = zeros(1, nf);

if nt == 1
    % Scalar target (the usual case): all frequencies at once
    if nz == 0
        Ht = Hs(it,:,:);                                   % 1 x ns x nf
        Hty = Hs(it,iy,:);
        PS = Vs(iy,iy) - Vs(iy,it)*(Vs(it,it) \ Vs(it,iy));
        Stt = real(qform(Ht, Vs));
        f = log(Stt) - log(Stt - real(qform(Hty, PS)));
    else
        Cr = Cs(ir,:);
        Hr = eye(nt+nz) + pagemtimes(Cr, pagemldivide(ZA, repmat(Kr,1,1,nf)));
        G = pagemldivide(Hr, Hs(ir,:,:));
        izy = [iz iy];
        PS = Vs(izy,izy) - Vs(izy,it)*(Vs(it,it) \ Vs(it,izy));
        Vrt = Vr(it,it);
        if nargout > 2
            chk.flat = max(abs(real(qform(G(it,:,:), Vs)) - Vrt)) / Vrt;
        end
        f = log(Vrt) - log(Vrt - real(qform(G(it,izy,:), PS)));
    end
    f = reshape(f, 1, nf);
elseif nz == 0
    % Unconditional: f = log|S_tt| - log|S_tt - H_ty Sig_y|t H_ty*|
    PS = Vs(iy,iy) - Vs(iy,it)*(Vs(it,it) \ Vs(it,iy));
    for k = 1:nf
        Ht = Hs(it,:,k);
        Stt = Ht*Vs*Ht';
        Hty = Hs(it,iy,k);
        f(k) = logdet(Stt) - logdet(Stt - Hty*PS*Hty');
    end
else
    % Conditional: filter [tgt cond] by the inverse reduced model; the
    % filtered target is white with covariance Vr(it,it)
    Cr = Cs(ir,:);
    Hr = eye(nt+nz) + pagemtimes(Cr, pagemldivide(ZA, repmat(Kr,1,1,nf)));
    G = pagemldivide(Hr, Hs(ir,:,:));  % rows: filtered [tgt cond]; cols: innovations of s
    izy = [iz iy];
    PS = Vs(izy,izy) - Vs(izy,it)*(Vs(it,it) \ Vs(it,izy));
    Vrt = Vr(it,it); ldr = logdet(Vrt);
    flat = 0;
    for k = 1:nf
        Gt = G(it,:,k);
        flat = max(flat, norm(Gt*Vs*Gt' - Vrt, 1) / norm(Vrt, 1));
        Gzy = G(it,izy,k);
        f(k) = ldr - logdet(Vrt - Gzy*PS*Gzy');
    end
    chk.flat = flat;
end
f = real(f);


function [Kr, Vr, res] = ss_reduce(A, C, K, V, r)
% Innovations model of the sub-process r of x_t = C z_t + e_t,
% z_{t+1} = A z_t + K e_t (steady-state Kalman filter of x_r alone)
Q = K*V*K'; S = K*V(:,r); R = V(r,r); Cr = C(r,:);
P = dare_sda(A, Cr, Q, R, S);
Vr = Cr*P*Cr' + R;  Vr = (Vr + Vr')/2;
APC = A*P*Cr' + S;
Kr = APC / Vr;
Pn = A*P*A' + Q - APC*(Vr \ APC');                 % Riccati residual
res = norm(Pn - P, 1) / max(norm(P, 1), eps);


function P = dare_sda(A, C, Q, R, S)
% Solves P = A P A' + Q - (A P C' + S)(C P C' + R)^-1 (A P C' + S)' by the
% structure-preserving doubling algorithm (quadratic convergence), after
% removing the cross term S
m = size(A,1);
RC = R \ C;
Ab = A - S*RC;
Qb = Q - S*(R \ S');  Qb = (Qb + Qb')/2;
Ak = Ab'; G = C'*RC; G = (G + G')/2; H = Qb;
I = eye(m);
for iter = 1:60                         % error ~ rho^(2^iter): 60 is never reached
    W = I + G*H;
    WA = W \ Ak;
    Gn = G + Ak*(W \ G)*Ak';  Gn = (Gn + Gn')/2;
    Hn = H + Ak'*H*WA;        Hn = (Hn + Hn')/2;
    Ak = Ak*WA;
    dH = norm(Hn - H, 1);
    G = Gn; H = Hn;
    if dH <= 1e-12*max(norm(H, 1), eps) || norm(Ak, 1) < 1e-14, break; end
end
P = H;


function q = qform(H, M)
% h*M*h' for each page of a 1 x k x nf array h, as 1 x 1 x nf
q = sum(pagemtimes(H, M) .* conj(H), 2);


function ld = logdet(M)
% log-determinant of a Hermitian positive definite matrix (Cholesky)
M = (M + M')/2;
[R, flag] = chol(M);
if flag
    ld = real(log(det(M)));        % not positive definite: report, don't crash
else
    ld = 2*sum(log(real(diag(R))));
end
