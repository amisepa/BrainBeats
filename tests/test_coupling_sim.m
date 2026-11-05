function R = test_coupling_sim(nSim, nboot, pick)
% TEST_COUPLING_SIM - End-to-end check of the brain-heart coupling analysis
% (GET_COUPLING_SERIES + COMPUTE_COUPLING) on simulated EEG, ECG and
% respiration with a known coupling (SIM_BRAINHEART; tier 3). No EEGLAB
% needed.
%
% Scenarios (posterior region, 4 bands, 5 min each):
%   null      - respiration drives alpha power and heart rate; no direct
%               alpha-heart link; CFA in the EEG; respiration recorded;
%               heartbeat-evoked EEG removed ('evoked')
%   null_cfa  - same data, heartbeat-evoked EEG left in ('none')
%   null_noR  - no respiration in the model
%   null_edr  - respiration derived from the ECG (QRS amplitude)
%   b2h       - alpha -> heart (2 s)
%   h2b       - heart -> alpha (1 s)
%   null_noCFA, null_qrsK - diagnostics: no cardiac field in the EEG; a
%               cardiac field of constant amplitude (no respiratory
%               modulation of the QRS)
%   null_ssp, null_both - the cardiac field removed by SSP, or SSP + evoked
% pick: cell array of scenario names to run (default all)
% Prints the rejection rates (p < .05) of the conditional tests (bootstrap
% GC, F-test, transfer entropy) and of the pairwise F-test, for alpha <->
% heart and beta <-> heart (beta carries most of the QRS power), and
% returns them in R. Checks: conditional alpha tests near 5% in 'null', and
% the planted directions detected.

if nargin < 1, nSim = 20; end
if nargin < 2, nboot = 50; end
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'functions'));
if nargin < 3 || isempty(pick), pick = 'all'; end
scen = {
    'null',     struct('c_ah',0,  'c_ha',0),   'channel', 'evoked'
    'null_cfa', struct('c_ah',0,  'c_ha',0),   'channel', 'none'
    'null_noR', struct('c_ah',0,  'c_ha',0),   'none',    'evoked'
    'null_edr', struct('c_ah',0,  'c_ha',0),   'ecg',     'evoked'
    'b2h',      struct('c_ah',0.5,'c_ha',0),   'channel', 'evoked'
    'h2b',      struct('c_ah',0,  'c_ha',0.3), 'channel', 'evoked'
    'null_noCFA', struct('c_ah',0,'c_ha',0,'cfa_uv',0),     'channel', 'evoked'
    'null_qrsK',  struct('c_ah',0,'c_ha',0,'qrs_resp',0),   'channel', 'evoked'
    'null_ssp',   struct('c_ah',0,'c_ha',0),                'channel', 'ssp'
    'null_both',  struct('c_ah',0,'c_ha',0),                'channel', 'both'};
if ~strcmp(pick, 'all'), scen = scen(ismember(scen(:,1), pick), :); end
reg = {'posterior', {'Pz' 'P3' 'P4' 'O1' 'O2'}};
cols = {'aH_boot','Ha_boot','aH_F','Ha_F','aH_TE','Ha_TE','aH_pair','Ha_pair','bH_F','Hb_F','bH_pair','Hb_pair'};
R = struct;
for s = 1:size(scen,1)
    P = nan(nSim, numel(cols));
    for i = 1:nSim
        nm = double(scen{s,1});
        cfg = scen{s,2}; cfg.seed = 1000*sum(nm.*(1:numel(nm))) + i;   % own seeds per scenario
        sim = sim_brainheart(cfg);
        switch scen{s,3}
            case 'channel', resp = struct('source','channel','signal',sim.resp);
            case 'ecg'
                amp = arrayfun(@(b) sim.ecg(b) - min(sim.ecg(b:min(b+round(0.06*sim.fs), end))), sim.beats);
                resp = struct('source','ecg','t',(sim.beats-1)/sim.fs,'amp',amp);
            otherwise, resp = struct('source','none');
        end
        S = get_coupling_series(sim.data, sim.fs, sim.chanlocs, sim.beats, sim.NN_t, sim.NN, resp, ...
            struct('regions', {reg}, 'cfa', scen{s,4}));
        w = warning('off','BrainBeats:coupling');
        out = compute_coupling(S, struct('nboot', nboot, 'nfreq', 64, 'seed', i));
        warning(w);
        t = out.table;
        g = @(src,tgt,mdl,meas,col) t.(col)(strcmp(t.source,src) & strcmp(t.target,tgt) & ...
            strcmp(t.model,mdl) & strcmp(t.freq,'all') & strcmp(t.measure,meas));
        P(i,:) = [g('alpha','heart','conditional','GC','p') g('heart','alpha','conditional','GC','p') ...
            g('alpha','heart','conditional','GC','p_F') g('heart','alpha','conditional','GC','p_F') ...
            g('alpha','heart','conditional','TE','p') g('heart','alpha','conditional','TE','p') ...
            g('alpha','heart','pairwise','GC','p_F') g('heart','alpha','pairwise','GC','p_F') ...
            g('beta','heart','conditional','GC','p_F') g('heart','beta','conditional','GC','p_F') ...
            g('beta','heart','pairwise','GC','p_F') g('heart','beta','pairwise','GC','p_F')];
    end
    rate = mean(P < 0.05, 1);
    R.(scen{s,1}) = cell2struct(num2cell(rate), cols, 2);
    tmp = [cols; num2cell(rate)];
    fprintf('%-9s', scen{s,1}); fprintf(' %s %.2f', tmp{:}); fprintf('\n');
end
lim = 0.05 + 3.29*sqrt(0.05*0.95/nSim);
if ~all(isfield(R, {'null' 'b2h' 'h2b'})), return; end
assert(R.null.aH_F <= lim && R.null.Ha_F <= lim, 'Conditional alpha-heart F-tests exceed %.2f in the null scenario', lim)
assert(R.b2h.aH_F >= 0.8, 'alpha -> heart missed (F-test)')
assert(R.h2b.Ha_F >= 0.8, 'heart -> alpha missed (F-test)')
if nboot > 0
    assert(R.null.aH_boot <= lim && R.null.Ha_boot <= lim, 'Conditional alpha-heart bootstrap tests exceed %.2f in the null scenario', lim)
    assert(R.b2h.aH_boot >= 0.8 && R.h2b.Ha_boot >= 0.8, 'Planted direction missed (bootstrap)')
end
disp('test_coupling_sim: checks passed')
