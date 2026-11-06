% PLOT_HEP_ICS - HEP of the independent components (IC level of RUN_HEP,
% 'hep_level' 'ics' or 'both'): for each IC, its scalp map and its HEP
% against the surrogate heartbeat trains (gray: 95% of the surrogates, red:
% p < .05 with the statistic chosen in RUN_HEP, max statistic by default), with its
% ICLabel class. All ICs left after the cleaning are
% shown, in the order of the decomposition: no IC is selected (whether a
% heartbeat-locked brain response can be separated from the cardiac field
% artifact is not established; see RUN_HEP). A heartbeat-locked IC can be the
% cardiac field artifact itself: check its map and its time course.
%
% Usage:
%   plot_hep_ics(HEP.brainbeats.ics)
%   figs = plot_hep_ics(HEP.brainbeats.ics, perFig)
%
% Inputs:
%   I      - HEP.brainbeats.ics (see RUN_HEP): ic, maps, chanlocs, iclabel,
%            iclabel_classes, hep, hep_times, surrogate (optional)
%   perFig - ICs per figure (default 24)
%
% Output:
%   figs   - figure handles
%
% Copyright (C) - Cedric Cannard, 2026

function figs = plot_hep_ics(I, perFig)

if nargin < 2 || isempty(perFig), perFig = 24; end
figs = gobjects(0);
if isempty(I) || isempty(I.ic), return; end
nIC = numel(I.ic);
t = I.hep_times;
hasS = isfield(I,'surrogate') && ~isempty(I.surrogate);
if hasS
    % statistic of the red areas (max statistic by default; family-wise over all ICs
    % and latencies), as chosen in RUN_HEP
    stat = 'fdr';
    if isfield(I.surrogate,'stat') && isfield(I.surrogate.hep, ['p_' I.surrogate.stat]), stat = I.surrogate.stat; end
    sigAll = I.surrogate.hep.(['p_' stat]) < .05;
    statNames = struct('tfce','TFCE p < .05, family-wise over all ICs', 'tmax','max-statistic p < .05, family-wise over all ICs', 'fdr','FDR p < .05');
    nNull = I.surrogate.nSurr;
    if isfield(I.surrogate,'nSurrHep'), nNull = I.surrogate.nSurrHep; end
end
nCol = 4;                                   % ICs per row (map + HEP each)
for first = 1:perFig:nIC
    idx = first:min(first+perFig-1, nIC);
    nRow = ceil(numel(idx)/nCol);
    name = 'HEP of the independent components';
    if nIC > perFig, name = sprintf('%s (%d-%d)', name, I.ic(idx(1)), I.ic(idx(end))); end
    fig = figure('color','w','Name',name,'NumberTitle','Off','Position',[40 40 1400 70+175*nRow]);
    try icadefs; set(fig,'color',BACKCOLOR); catch; end
    tl = tiledlayout(fig, nRow, 3*nCol, 'TileSpacing','compact', 'Padding','compact');
    for k = idx
        % scalp map (1 tile), HEP (2 tiles)
        nexttile(tl);
        topoplot(I.maps(:,k), I.chanlocs, 'electrodes','off', 'style','map', 'gridscale',40);
        ax = nexttile(tl, [1 2]); hold(ax,'on')
        y = I.hep(k,:);
        if hasS
            S = I.surrogate.hep;
            yl = [min([S.null_lo(k,:) y]) max([S.null_hi(k,:) y])];
            yl = yl + [-1 1]*0.08*diff(yl);
            d = diff([false sigAll(k,:) false]);
            on = find(d == 1); off = find(d == -1) - 1;
            for j = 1:numel(on)
                patch(ax, [t(on(j)) t(off(j)) t(off(j)) t(on(j))], yl([1 1 2 2]), [1 .8 .8], 'EdgeColor','none');
            end
            fill(ax, [t fliplr(t)], [S.null_lo(k,:) fliplr(S.null_hi(k,:))], [.8 .8 .8], 'EdgeColor','none', 'FaceAlpha',.8);
            ylim(ax, yl)
        end
        plot(ax, t, y, 'k', 'LineWidth', 1.2);
        xline(ax, 0, 'k--');
        xlim(ax, t([1 end])); box(ax,'on')
        lab = '';
        if isfield(I,'iclabel') && ~isempty(I.iclabel)
            [pm, c] = max(I.iclabel(k,:));
            lab = sprintf(': %s %.0f%%', I.iclabel_classes{c}, 100*pm);
        end
        title(ax, sprintf('IC %d%s', I.ic(k), lab), 'FontSize', 9);
        if k == idx(1)
            ylabel(ax, 'HEP (IC activation)');
        end
        if k > idx(end) - nCol
            xlabel(ax, 'Latency (ms)');
        end
        set(ax, 'FontSize', 8);
    end
    if hasS
        title(tl, sprintf(['HEP of the independent components (scalp map, ICLabel class): gray = 95%% of %d ' ...
            'surrogate heartbeat trains, red = %s'], nNull, statNames.(stat)), ...
            'FontSize', 11, 'FontWeight', 'bold');
    else
        title(tl, 'HEP of the independent components (scalp map, ICLabel class)', 'FontSize', 11, 'FontWeight', 'bold');
    end
    colormap(fig, parula)
    finish_figure(fig)
    figs(end+1) = fig; %#ok<AGROW>
end
end
