function T = summarize_cohort_coupling(resFile)
% SUMMARIZE_COHORT_COUPLING - Group summary of COHORT_COUPLING results.
%
% Prints, per analysis variant: the share of recordings with a significant
% conditional GC (F-test, p < .05; ~5% expected without coupling) for each
% band, direction and region, pooled over regions; the respiration checks
% (HF coherence of each respiration proxy with NN; respiration -> heart);
% the vascular links (PPG pulse arrival time and amplitude); and the model
% diagnostics (Ljung-Box failures, samples per coefficient).
%
% Usage:
%   T = summarize_cohort_coupling(resFile)
%
% Output T: table of the shares per variant, source, target, measure.
%
% Copyright (C) - Cedric Cannard, 2026

L = load(resFile, 'R');
R = L.R;
fn = fieldnames(R);
V = struct('name', fn(cellfun(@(f) isstruct(R(1).(f)) && isfield(R(1).(f), 'table'), fn)));
n = numel(R);
fprintf('%d recordings\n', n);
rows = {};
for v = 1:numel(V)
    vn = V(v).name;
    ok = arrayfun(@(r) isfield(r.(vn), 'table'), R);
    if ~any(ok), continue; end
    tabs = arrayfun(@(r) r.(vn).table, R(ok), 'UniformOutput', false);
    t1 = tabs{1};
    sel = strcmp(t1.measure,'GC') & strcmp(t1.freq,'all') & ismember(t1.model, {'conditional' 'exogenous'});
    keys = unique(strcat(t1.source(sel), '>', t1.target(sel), '|', t1.model(sel)), 'stable');
    fprintf('\n== %s (%d recordings) ==\n', vn, sum(ok));
    for k = 1:numel(keys)
        parts = split(keys{k}, {'>' '|'});
        hits = 0; tot = 0;
        for i = 1:numel(tabs)
            t = tabs{i};
            m = strcmp(t.measure,'GC') & strcmp(t.freq,'all') & strcmp(t.source, parts{1}) & ...
                strcmp(t.target, parts{2}) & strcmp(t.model, parts{3});
            p = t.p_F(m);
            hits = hits + sum(p < .05); tot = tot + sum(isfinite(p));
        end
        share = hits/max(tot,1);
        pb = 1 - binocdf(hits-1, tot, 0.05);           % above 5%?
        rows(end+1,:) = {vn, parts{1}, parts{2}, parts{3}, share, tot, pb}; %#ok<AGROW>
        flag = ''; if pb < 0.001, flag = '  ***'; elseif pb < 0.01, flag = '  **'; end
        fprintf('  %-10s -> %-8s %-11s %5.1f%% of %3d%s\n', parts{1}, parts{2}, parts{3}, 100*share, tot, flag);
    end
    % respiration and diagnostics
    coh = cell2mat(arrayfun(@(r) r.(vn).resp_heart_coh(:)', R(ok), 'UniformOutput', false)');
    if ~isempty(coh) && any(isfinite(coh(:)))
        fprintf('  respiration-NN HF coherence (median [IQR]) per proxy:');
        for j = 1:size(coh,2)
            q = prctile(coh(:,j), [25 50 75]);
            fprintf(' %.2f [%.2f %.2f]', q(2), q(1), q(3));
        end
        fprintf('\n');
    end
    lb = cell2mat(arrayfun(@(r) vertcat(r.(vn).models.ljungbox_p), R(ok), 'UniformOutput', false)');
    spc = cell2mat(arrayfun(@(r) [r.(vn).models.samples_per_coef], R(ok), 'UniformOutput', false));
    vars = R(find(ok,1)).(vn).models(1).variables;
    fprintf('  Ljung-Box p < .05: %s (%s); samples/coef median %.0f\n', ...
        mat2str(round(100*mean(lb < .05, 1)), 2), strjoin(vars, ' '), median(spc));
end
T = cell2table(rows, 'VariableNames', {'variant','source','target','model','share','n','p_binom'});
