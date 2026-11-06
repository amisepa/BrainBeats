function R = cohort_coupling(cacheFiles, outFile, nboot, cfaAll)
% COHORT_COUPLING - Brain-heart coupling of a cohort under several analysis
% variants, from the caches of CACHE_COUPLING_DATA (no EEGLAB needed).
%
% Variants (series computed once per cardiac field option):
%   default    respiration from the ECG (R-S amplitude), cardiac field 'evoked'
%   resp_none  no respiration in the model
%   resp_ppg   respiration from the PPG only (pulse foot and amplitude)
%   resp_all   ECG- and PPG-derived respiration together
%   vasc       default + PPG pulse arrival time and pulse amplitude modelled
%   ssp, both, nocfa   cardiac field removed by SSP, SSP + evoked, not
%              removed (respiration from the ECG; only with cfaAll)
%
% Usage:
%   R = cohort_coupling(cacheFiles, outFile, nboot, cfaAll)
%
% Inputs:
%   cacheFiles - cell array of *_coupling_cache.mat files
%   outFile    - .mat file for the results (saved after each recording;
%                recordings already in it are skipped)
%   nboot      - bootstrap samples (default 0: F-tests and transfer entropy)
%   cfaAll     - also the cardiac field variants (default false)
%
% Output R: struct array, one per recording: name, the series (S, evoked),
%   and per variant the coupling table, model diagnostics and
%   respiration-heart coherence.
%
% Copyright (C) - Cedric Cannard, 2026

if nargin < 3, nboot = 0; end
if nargin < 4, cfaAll = false; end
cfas = {'evoked'}; if cfaAll, cfas = {'evoked' 'ssp' 'both' 'none'}; end
cfaName = struct('evoked','default', 'ssp','ssp', 'both','both', 'none','nocfa');
R = struct('name', {});
if exist(outFile, 'file'), load(outFile, 'R'); end
done = {R.name};
w = warning('off', 'BrainBeats:coupling');
copt = struct('nboot', nboot, 'nfreq', 64);
for k = 1:numel(cacheFiles)
    C = load(cacheFiles{k});
    if ismember(C.name, done), continue; end
    t0 = tic;
    r = struct('name', C.name, 'S', []);
    vn = {'default' 'resp_none' 'resp_ppg' 'resp_all' 'vasc'};
    if cfaAll, vn = [vn {'ssp' 'both' 'nocfa'}]; end
    for v = vn, r.(v{1}) = struct('error', 'not run'); end   % same fields for every recording
    wEdr = round(0.06*C.fs);
    b = C.beats; s = double(C.ecg_sig);
    resp = struct('source', {'ecg' 'ppg_riiv' 'ppg_riav'}, 'signal', [], ...
        't', {(b-1)/C.fs [] []}, 'amp', {arrayfun(@(x) s(x) - min(s(x:min(x+wEdr, numel(s)))), b) [] []});
    ppg = struct('onsets', C.ppg_onsets, 'raw', C.ppg, 'sig', C.ppg_sig, 'window', [0.04 0.5]);
    for c = 1:numel(cfas)
        try
            S = get_coupling_series(double(C.X), C.fs, C.chanlocs, C.beats, C.NN_t, C.NN, resp, ...
                struct('cfa', cfas{c}, 'ppg', ppg, 'vascular', {{'pat' 'pamp'}}));
            % respiration proxies by name (a poor PPG proxy is left out of S)
            ecg = {'ecg'}; ppgR = {'ppg_riiv' 'ppg_riav'};
            if strcmp(cfas{c}, 'evoked')
                r.S = S;
                variants = {'default', ecg, false; 'resp_none', {}, false; 'resp_ppg', ppgR, false; ...
                    'resp_all', [ecg ppgR], false; 'vasc', ecg, true};
            else
                variants = {cfaName.(cfas{c}), ecg, false};
            end
        catch e
            fprintf('%s %s: %s\n', C.name, cfas{c}, e.message);
            continue
        end
        for v = 1:size(variants,1)                        % a failing variant leaves the others
            try
                Sv = select(S, variants{v,2}, variants{v,3});
                out = compute_coupling(Sv, copt);
                r.(variants{v,1}) = struct('table', out.table, 'models', out.models, ...
                    'resp_heart_coh', out.resp_heart_coh, 'warnings', {out.warnings}, 'cfa', S.cfa);
            catch e
                fprintf('%s %s: %s\n', C.name, variants{v,1}, e.message);
                r.(variants{v,1}) = struct('error', e.message);
            end
        end
    end
    if isempty(R), R = r; else, R(end+1) = r; end %#ok<AGROW>
    save(outFile, 'R', '-v7.3');
    fprintf('%s: done in %.0f s\n', C.name, toc(t0));
end
warning(w);


function Sv = select(S, want, vasc)
% Respiration proxies named in want (those present); vascular series kept or not
names = strsplit(S.resp_source, '+');
rr = find(ismember(names, want));
if ~isempty(want) && isempty(rr), error('respiration %s not available', strjoin(want, '+')); end
Sv = S;
Sv.resp = S.resp(rr, :);
if isempty(rr), Sv.resp_source = 'none'; else, Sv.resp_source = strjoin(names(rr), '+'); end
if ~vasc, Sv.vasc = []; Sv.vasc_names = {}; end
