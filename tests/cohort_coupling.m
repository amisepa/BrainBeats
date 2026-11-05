function R = cohort_coupling(cacheFiles, outFile, nboot)
% COHORT_COUPLING - Brain-heart coupling of a cohort under several analysis
% variants, from the caches of CACHE_COUPLING_DATA (no EEGLAB needed).
%
% Variants:
%   default    respiration from the ECG (R-S amplitude), cardiac field 'evoked'
%   ssp, both, nocfa   cardiac field removed by SSP, SSP + evoked, not removed
%   resp_none  no respiration in the model
%   resp_ppg   respiration from the PPG only (pulse foot and amplitude)
%   resp_all   ECG- and PPG-derived respiration together
%   vasc       default + PPG pulse arrival time and pulse amplitude modelled
%
% Usage:
%   R = cohort_coupling(cacheFiles, outFile, nboot)
%
% Inputs:
%   cacheFiles - cell array of *_coupling_cache.mat files
%   outFile    - .mat file for the results (saved after each recording;
%                recordings already in it are skipped)
%   nboot      - bootstrap samples (default 0: F-tests and transfer entropy)
%
% Output R: struct array, one per recording: name, and per variant the
%   coupling table, model diagnostics and respiration-heart coherence.
%
% Copyright (C) - Cedric Cannard, 2026

if nargin < 3, nboot = 0; end
V = struct('name', {'default' 'ssp' 'both' 'nocfa' 'resp_none' 'resp_ppg' 'resp_all' 'vasc'}, ...
    'resp', {{'ecg'} {'ecg'} {'ecg'} {'ecg'} {} {'ppg_riiv' 'ppg_riav'} {'ecg' 'ppg_riiv' 'ppg_riav'} {'ecg'}}, ...
    'cfa', {'evoked' 'ssp' 'both' 'none' 'evoked' 'evoked' 'evoked' 'evoked'}, ...
    'vascular', {{} {} {} {} {} {} {} {'pat' 'pamp'}});
R = struct('name', {});
if exist(outFile, 'file'), load(outFile, 'R'); end
done = {R.name};
w = warning('off', 'BrainBeats:coupling');
for k = 1:numel(cacheFiles)
    C = load(cacheFiles{k});
    if ismember(C.name, done), continue; end
    t0 = tic;
    r = struct('name', C.name);
    wEdr = round(0.06*C.fs);
    b = C.beats; s = double(C.ecg_sig);
    edr = struct('source','ecg', 't', (b-1)/C.fs, 'amp', arrayfun(@(x) s(x) - min(s(x:min(x+wEdr, numel(s)))), b));
    ppg = struct('onsets', C.ppg_onsets, 'raw', C.ppg, 'sig', C.ppg_sig, 'window', [0.04 0.5]);
    for v = 1:numel(V)
        resp = struct('source', 'none', 't', [], 'amp', []);
        for j = 1:numel(V(v).resp)
            resp(j) = struct('source', V(v).resp{j}, 't', [], 'amp', []);
            if strcmp(V(v).resp{j}, 'ecg'), resp(j).t = edr.t; resp(j).amp = edr.amp; end
        end
        try
            S = get_coupling_series(double(C.X), C.fs, C.chanlocs, C.beats, C.NN_t, C.NN, resp, ...
                struct('cfa', V(v).cfa, 'ppg', ppg, 'vascular', {V(v).vascular}));
            out = compute_coupling(S, struct('nboot', nboot, 'nfreq', 64));
            r.(V(v).name) = struct('table', out.table, 'models', out.models, ...
                'resp_heart_coh', out.resp_heart_coh, 'warnings', {out.warnings}, 'cfa', S.cfa);
        catch e
            r.(V(v).name) = struct('error', e.message);
            fprintf('%s %s: %s\n', C.name, V(v).name, e.message);
        end
    end
    if isempty(R), R = r; else, R(end+1) = r; end %#ok<AGROW>
    save(outFile, 'R', 'V');
    fprintf('%s: %d variants in %.0f s\n', C.name, numel(V), toc(t0));
end
warning(w);
