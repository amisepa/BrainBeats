function test_coupling_stats(nSim)
% TEST_COUPLING_STATS - Calibration of the tests of COMPUTE_COUPLING on
% simulated 2-Hz series (5 min) with a known structure. No EEGLAB needed.
%
%   Null: respiration drives both the EEG band power and the heart period,
%   with no direct link between them. The conditional tests (bootstrap GC,
%   F-test, transfer entropy) must reject at ~5%; the pairwise F-test
%   rejects far more often (it sees the common drive), which is what
%   conditioning is for.
%   Effect: a direct EEG -> heart link (lag 2 s) is added; the conditional
%   tests must detect it, and not the reverse direction more than ~5%.
%
% nSim: datasets per condition (default 200). Errors if a rate is outside
% its binomial 99.9% range around 5%.

if nargin < 1, nSim = 200; end
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'functions'));
nw = 0;                                % parfor workers: only an open pool is used
if ~isempty(gcp('nocreate')), nw = gcp('nocreate').NumWorkers; end
fso = 2; T = 300*fso;

for cond = 1:2
    c = [0 0.25];
    [A, V] = model(c(cond));
    P = nan(nSim, 7);   % [boot E->H, boot H->E, F E->H, F H->E, pairF E->H, pairF H->E, TE E->H]
    parfor (s = 1:nSim, nw)
        rng(1000*cond + s);
        opt = struct('nboot', 100, 'nfreq', 64, 'maxlag', 5, 'seed', s);
        X = var_sim(A, V, T, 1, 500);
        S = struct('t', (0:T-1)/fso, 'eeg', reshape(X(1,:),1,1,T), 'heart', X(2,:), 'resp', X(3,:), 'resp_source', 'channel', ...
            'regions', struct('name','sim','chans',1,'labels',{{'sim'}}), 'bands', [8 13], ...
            'band_names', {{'alpha'}}, 'fs_out', fso);
        out = evalc_quiet(S, opt);
        t = out.table;
        g = @(src, tgt, mdl, meas, col) t.(col)(strcmp(t.source,src) & strcmp(t.target,tgt) & ...
            strcmp(t.model,mdl) & strcmp(t.freq,'all') & strcmp(t.measure,meas));
        P(s,:) = [g('alpha','heart','conditional','GC','p') g('heart','alpha','conditional','GC','p') ...
            g('alpha','heart','conditional','GC','p_F') g('heart','alpha','conditional','GC','p_F') ...
            g('alpha','heart','pairwise','GC','p_F') g('heart','alpha','pairwise','GC','p_F') ...
            g('alpha','heart','conditional','TE','p')];
    end
    rate = mean(P < 0.05, 1);
    fprintf('%s: rejection rate at .05 (n = %d)\n', ternary(cond==1,'NULL (respiration drives both)','EFFECT (alpha -> heart)'), nSim);
    fprintf('   cond. GC bootstrap  E->H %.3f  H->E %.3f\n', rate(1), rate(2));
    fprintf('   cond. GC F-test     E->H %.3f  H->E %.3f\n', rate(3), rate(4));
    fprintf('   pairwise F-test     E->H %.3f  H->E %.3f\n', rate(5), rate(6));
    fprintf('   cond. TE (copula)   E->H %.3f\n', rate(7));
    lim = 0.05 + [-1 1]*3.29*sqrt(0.05*0.95/nSim);
    if cond == 1
        chk = rate([1 2 3 4 7]);
        assert(all(chk >= lim(1) & chk <= lim(2)), 'Conditional tests not calibrated: %s (99.9%% range %.3f-%.3f)', mat2str(chk,3), lim)
        assert(min(rate(5:6)) > 0.2, 'Pairwise tests should see the common drive')
    else
        assert(all(rate([1 3 7]) > 0.8), 'Direct alpha -> heart link missed: %s', mat2str(rate([1 3 7]),3))
        assert(all(rate([2 4]) <= lim(2)), 'Reverse direction falsely detected: %s', mat2str(rate([2 4]),3))
    end
end
disp('test_coupling_stats: all checks passed')


function out = evalc_quiet(S, opt)
w = warning('off', 'BrainBeats:coupling');
out = compute_coupling(S, opt);
warning(w);


function [A, V] = model(cEH)
% 1 = alpha power, 2 = heart period, 3 = respiration (2 Hz)
n = 3; A1 = zeros(n); A2 = zeros(n); A4 = zeros(n);
w = 2*pi*0.25/2; r = 0.85;
A1(3,3) = 2*r*cos(w); A2(3,3) = -r^2;          % breathing at 0.25 Hz
A1(1,1) = 0.8; A1(1,3) = 0.35;                  % alpha: slow, respiration-modulated
A1(2,2) = 0.6; A1(2,3) = 0.5; A2(2,2) = -0.2;   % heart: RSA
A4(2,1) = cEH;                                  % alpha -> heart at lag 2 s
A = [A1 A2 zeros(n) A4];
V = eye(n);


function o = ternary(c, a, b)
if c, o = a; else, o = b; end
