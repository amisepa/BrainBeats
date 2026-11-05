function test_coupling_engine
% TEST_COUPLING_ENGINE - Checks of the VAR engine of the brain-heart
% coupling analysis (var_fit, var_order, var_ssgc, var_spectral, var_sim)
% on a VAR with known couplings (tier 1). No EEGLAB needed. Errors on the
% first failed check.
%
% Model (4 variables, order 2, as 2-Hz series): 1 = EEG band power (E),
% 2 = heart period (H), 3 = respiration (R), 4 = a second EEG band (E2).
%   R: 0.25-Hz oscillation; R -> H (respiratory sinus arrhythmia) and
%   R -> E (respiration-locked EEG power); E -> H (brain to heart, lag 2);
%   H -> E2 (heart to brain, lag 1). No H -> E: the pairwise H -> E
%   Granger causality is the respiratory common drive, the conditional one
%   is zero.
%   1. Riccati residuals ~0; the reduced innovations are white; spectral GC
%      averages to the time-domain GC (pairwise and conditional).
%   2. True GC: conditional H -> E = 0, pairwise H -> E > 0 (mediated);
%      E -> H, H -> E2, R -> H > 0. GPDC and DC are 0 where the model has
%      no direct / no path; DC rows and GPDC columns sum to 1.
%   3. State-space GC equals regression GC of a long (order 40) model
%      fitted to 400 000 simulated samples (independent check of the
%      Riccati reduction).
%   4. Estimation from 200 000 samples recovers the true GC; BIC finds
%      order 2; gaps (NaN) are skipped without bias.

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'functions'));
rng(11);

[A, V] = test_model();
n = 4; E = 1; H = 2; R = 3; E2 = 4;
nfreq = 256;

% 1-2. exact values from the true model
pairs = {H E [R E2]; H E []; E H [R E2]; E H []; R H [E E2]; H E2 [E R]; E2 H [E R]};
Ftrue = zeros(size(pairs,1),1);
for k = 1:size(pairs,1)
    [F, f, chk] = var_ssgc(A, V, pairs{k,2}, pairs{k,1}, pairs{k,3}, nfreq);
    Ftrue(k) = F;
    fm = trapz(f)/nfreq;   % mean over [0, pi]
    assert(chk.dare < 1e-10, 'Riccati residual %.2g (pair %d)', chk.dare, k)
    assert(isnan(chk.flat) || chk.flat < 1e-8, 'Reduced innovations not white: %.2g (pair %d)', chk.flat, k)
    assert(abs(fm - F) <= 1e-6*max(F,1e-3), 'Spectral GC mean %.6g vs time-domain %.6g (pair %d)', fm, F, k)
    assert(all(f > -1e-10), 'Negative spectral GC (pair %d)', k)
end
assert(Ftrue(1) < 1e-12, 'Conditional H->E should be 0, got %.3g', Ftrue(1))
assert(Ftrue(2) > 0.01, 'Pairwise H->E should be > 0 (respiratory common drive), got %.3g', Ftrue(2))
assert(all(Ftrue([3 5 6]) > 0.02), 'E->H, R->H, H->E2 should be > 0')
assert(Ftrue(7) < 1e-12, 'Conditional E2->H should be 0, got %.3g', Ftrue(7))
fprintf('1-2 OK: true GC  H->E|R,E2 %.1e, H->E pairwise %.3f, E->H|. %.3f, R->H|. %.3f, H->E2|. %.3f\n', Ftrue([1 2 3 5 6]));

sp = var_spectral(A, V, nfreq);
assert(max(abs(sum(sp.dc,2)-1),[],'all') < 1e-12 && max(abs(sum(sp.gpdc,1)-1),[],'all') < 1e-12, 'DC rows / GPDC columns do not sum to 1')
assert(max(sp.gpdc(E,H,:)) < 1e-20 && max(sp.dc(E,H,:)) < 1e-20, 'GPDC/DC H->E should be 0 (no path)')
assert(max(sp.gpdc(H,E2,:)) < 1e-20, 'GPDC E2->H should be 0 (no direct link)')
assert(max(sp.dc(H,E2,:)) < 1e-20, 'DC E2->H should be 0 (no path)')
assert(max(sp.dc(E2,E,:)) > 1e-4, 'DC E->E2 should be > 0 (indirect path E->H->E2)')
assert(all(sp.coh(:) >= -1e-12 & sp.coh(:) <= 1+1e-12) && all(sp.pcoh(:) <= 1+1e-12), 'Coherences outside [0 1]')
assert(max(sp.coh(E,H,:)) > 0.05 && max(sp.pcoh(E,H,:)) > 0.05, 'E-H coherence and partial coherence should be > 0')
fprintf('   DC/GPDC zeros where there is no path / no direct link; DC E->E2 (indirect) max %.3f\n', max(sp.dc(E2,E,:)));

% 3. Riccati reduction vs long regressions on simulated data
X = var_sim(A, V, 4e5, 1, 2000);
pl = 40;
Fols = zeros(4,1); idx = [1 2 3 6];
for k = 1:4
    pk = pairs(idx(k),:);
    [x, y, z] = deal(pk{2}, pk{1}, pk{3});
    Fols(k) = log(rss(X([x z],:), pl, 1) / rss(X([x y z],:), pl, 1));
end
d = abs(Fols - Ftrue(idx));
assert(all(d < 0.006), 'State-space GC vs order-%d regression: %s vs %s', pl, mat2str(Ftrue(idx)',3), mat2str(Fols',3))
fprintf('3 OK: state-space GC = order-%d regression GC within %.4f\n', pl, max(d));

% 4. estimation
X = var_sim(A, V, 2e5, 1, 2000);
pB = var_order(X, 8);
assert(pB == 2, 'BIC order %d (expected 2)', pB)
[Ah, Vh] = var_fit(X, 2);
Fh = zeros(size(Ftrue));
for k = 1:size(pairs,1)
    Fh(k) = var_ssgc(Ah, Vh, pairs{k,2}, pairs{k,1}, pairs{k,3});
end
assert(all(abs(Fh - Ftrue) < 0.004), 'Estimated GC %s vs true %s', mat2str(Fh',3), mat2str(Ftrue',3))
Xg = X; Xg(:, 50000:50100) = NaN; Xg(2, 120000:120003) = NaN;     % gaps
[Ag, ~, info] = var_fit(Xg, 2);
assert(info.nObs == size(X,2) - 2 - 103 - 6, 'Gap handling: %d time points used', info.nObs)
assert(max(abs(Ag(:) - Ah(:))) < 0.01, 'Gaps changed the fit')
fprintf('4 OK: BIC order %d; estimated GC within %.4f of the truth; gaps skipped (%d points)\n', pB, max(abs(Fh - Ftrue)), info.nObs);
disp('test_coupling_engine: all checks passed')


function [A, V] = test_model()
n = 4; p = 2;
A1 = zeros(n); A2 = zeros(n);
w = 2*pi*0.25/2; r = 0.9;                % respiration: 0.25 Hz at fs = 2 Hz
A1(3,3) = 2*r*cos(w); A2(3,3) = -r^2;
A1(2,2) = 0.5;  A1(2,3) = 0.4;  A2(2,1) = 0.3;    % H: self, R -> H, E -> H (lag 2)
A1(1,1) = 0.7;  A1(1,3) = 0.3;                    % E: self, R -> E
A1(4,4) = 0.6;  A1(4,2) = 0.25;                   % E2: self, H -> E2
A = [A1 A2];
V = eye(n);


function r = rss(X, p, target)
% residual variance of variable 'target' regressed on p lags of all rows
[n, N] = size(X);
Z = zeros(n*p, N-p);
for k = 1:p
    Z((k-1)*n+1:k*n, :) = X(:, p+1-k:N-k);
end
y = X(target, p+1:N);
b = Z' \ y';
r = mean((y - b'*Z).^2);
