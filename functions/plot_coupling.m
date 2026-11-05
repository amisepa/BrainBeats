% PLOT_COUPLING - Figures of the brain-heart coupling analysis.
%
% Usage:
%   plot_coupling(out, S)
%
% INPUTS:
%   out - output of COMPUTE_COUPLING
%   S   - series of GET_COUPLING_SERIES (optional: adds the series figure)
%
% Figures:
%   1. Conditional Granger causality, EEG -> heart and heart -> EEG
%      (regions x bands). * = significant (p < .05, bootstrap, or F-test
%      without bootstrap); o = significant pairwise only: the link runs
%      through the other variables (respiration, other bands).
%   2. Spectral conditional GC over the modulation frequencies (0-1 Hz;
%      LF and HF shaded), one line per region, dotted: 95th percentile of
%      the bootstrap null.
%   3. The series themselves (first region): band power envelopes, NN and
%      respiration, z-scored.
%
% Copyright (C) - Cedric Cannard, BrainBeats 2026

function plot_coupling(out, S)

t = out.table;
regs = unique(t.region, 'stable');
bands = unique(t.band(~strcmp(t.band,'')), 'stable');
nr = numel(regs); nb = numel(bands);
figs0 = findall(groot, 'Type', 'figure');

% 1. Conditional GC maps
figure('Color','w','Name','Brain-heart coupling: conditional Granger causality','NumberTitle','off');
dirs = {'EEG -> heart', 'heart -> EEG'};
G = nan(nr, nb, 2); sig = false(nr, nb, 2); med = false(nr, nb, 2);
for r = 1:nr
    for b = 1:nb
        for d = 1:2
            if d == 1, src = bands{b}; tgt = 'heart'; else, src = 'heart'; tgt = bands{b}; end
            sel = strcmp(t.region, regs{r}) & strcmp(t.source, src) & strcmp(t.target, tgt) & ...
                strcmp(t.measure, 'GC') & strcmp(t.freq, 'all');
            c = sel & strcmp(t.model, 'conditional'); pw = sel & strcmp(t.model, 'pairwise');
            if ~any(c), continue; end
            G(r,b,d) = t.value(c);
            pc = t.p(c); if isnan(pc), pc = t.p_F(c); end
            sig(r,b,d) = pc < 0.05;
            med(r,b,d) = ~sig(r,b,d) && t.p(pw) < 0.05;
        end
    end
end
cmax = max(G(:)); if ~(cmax > 0), cmax = 1; end
for d = 1:2
    subplot(1,2,d)
    imagesc(G(:,:,d), [0 cmax]); colormap(gca, parula); hold on
    [rr, bb] = find(sig(:,:,d)); text(bb, rr, '*', 'HorizontalAlignment','center', 'FontSize',18, 'Color','w', 'FontWeight','bold');
    [rr, bb] = find(med(:,:,d)); plot(bb, rr, 'wo', 'MarkerSize', 10, 'LineWidth', 1.5);
    set(gca, 'XTick', 1:nb, 'XTickLabel', bands, 'YTick', 1:nr, 'YTickLabel', regs, 'TickLength', [0 0]);
    title(dirs{d}); cb = colorbar; ylabel(cb, 'GC (log ratio)');
end
sgtitle({'Conditional Granger causality (given the other bands and respiration)', ...
    '* significant;  o pairwise only (mediated by the other variables)'}, 'FontSize', 10);

% 2. Spectral GC
sp = out.spectra; f = out.freqs;
figure('Color','w','Name','Brain-heart coupling: spectral Granger causality','NumberTitle','off');
col = lines(nr);
for d = 1:2
    for b = 1:nb
        subplot(2, nb, (d-1)*nb + b); hold on
        if d == 1, Y = sp.gc_eh; N95 = sp.gc_eh_null95; else, Y = sp.gc_he; N95 = sp.gc_he_null95; end
        yl = max([Y(:); N95(:); eps]);
        patch([out.lf fliplr(out.lf)], [0 0 yl yl]*1.05, [.9 .9 .9], 'EdgeColor','none', 'FaceAlpha', .6);
        patch([out.hf fliplr(out.hf)], [0 0 yl yl]*1.05, [.85 .9 1], 'EdgeColor','none', 'FaceAlpha', .6);
        for r = 1:nr
            plot(f, squeeze(Y(r,b,:)), 'Color', col(r,:), 'LineWidth', 1.5);
            if any(isfinite(N95(r,b,:)))
                plot(f, squeeze(N95(r,b,:)), ':', 'Color', col(r,:));
            end
        end
        xlim([0 min(1, f(end))]); ylim([0 yl*1.05]); box on
        if d == 1, title(bands{b}); end
        if b == 1, ylabel(sprintf('%s\nspectral GC', dirs{d})); end
        if d == 2, xlabel('Modulation frequency (Hz)'); end
    end
end
legend(findobj(gca,'Type','line','LineStyle','-'), flipud(regs(:)), 'Location', 'northeast', 'FontSize', 8);
sgtitle('Spectral conditional GC (grey: LF, blue: HF; dotted: bootstrap 95%)', 'FontSize', 10);

% 3. Series
if nargin > 1 && ~isempty(S)
    figure('Color','w','Name','Brain-heart coupling: series','NumberTitle','off');
    X = [reshape(S.eeg(1,:,:), size(S.eeg,2), []); S.heart];
    names = [S.band_names(:)' {'NN'}];
    if ~isempty(S.resp), X = [X; S.resp]; names{end+1} = sprintf('resp (%s)', S.resp_source); end
    Z = (X - mean(X,2,'omitnan')) ./ std(X,0,2,'omitnan');
    plot(S.t, Z' - 4*(0:size(Z,1)-1)); box on
    set(gca, 'YTick', fliplr(-4*(0:size(Z,1)-1)), 'YTickLabel', fliplr(names));
    xlabel('Time (s)'); xlim([S.t(1) S.t(end)]);
    title(sprintf('Series of region %s (z-scored, %g Hz)', S.regions(1).name, S.fs_out))
end

figs = findall(groot, 'Type', 'figure');
finish_figure(figs(~ismember(figs, figs0)))
