% REPORT_COUPLING - One-page report of the brain-heart coupling analysis:
% what was coupled, in which direction, directly or through respiration and
% the other EEG bands, in which frequency range, and how far to trust it.
%
% Usage:
%   lines = report_coupling(out, S)
%   lines = report_coupling(out, S, show)
%
% INPUTS:
%   out  - output of COMPUTE_COUPLING
%   S    - series of GET_COUPLING_SERIES
%   show - draw the figure (default true); the text is always returned
%
% OUTPUT:
%   lines - the interpretation, as a cell array of sentences (also printed)
%
% The figure shows, per region, the links between the EEG bands and the
% heart period as arrows: solid = direct (significant conditional Granger
% causality after FDR across all the recording's tests: not explained by
% respiration or the other bands); dashed grey = pairwise only (the
% association runs through respiration or the other bands); width = GC;
% colour = the frequency range carrying it (LF 0.04-0.15 Hz, HF 0.15-0.4
% Hz). Below: the strongest links, the quality checks, and the text.
%
% Copyright (C) - Cedric Cannard, BrainBeats 2026

function lines = report_coupling(out, S, show)

if nargin < 3, show = true; end
t = out.table;
if ~ismember('q', t.Properties.VariableNames)            % results from before the FDR column
    t.q = t.p_F;
end
regs = {S.regions.name};
bands = S.band_names(:)';
nr = numel(regs); nb = numel(bands);
pval = @(m) m.p_F;                                         % as the FDR (q)

% Link classification
L = struct('region', {}, 'band', {}, 'dir', {}, 'gc', {}, 'q', {}, 'p', {}, 'p_pair', {}, ...
    'lf', {}, 'hf', {}, 'class', {}, 'trgc_ok', {});
for r = 1:nr
    for b = 1:nb
        tr = t(strcmp(t.region, regs{r}) & strcmp(t.band, bands{b}) & strcmp(t.measure,'TRGC'), :);
        for d = 1:2
            if d == 1, src = bands{b}; tgt = 'heart'; else, src = 'heart'; tgt = bands{b}; end
            m = strcmp(t.region, regs{r}) & strcmp(t.source, src) & strcmp(t.target, tgt) & strcmp(t.measure,'GC');
            c = t(m & strcmp(t.model,'conditional') & strcmp(t.freq,'all'), :);
            if isempty(c), continue; end
            pw = t(m & strcmp(t.model,'pairwise') & strcmp(t.freq,'all'), :);
            lf = t.value(m & strcmp(t.model,'conditional') & strcmp(t.freq,'LF'));
            hf = t.value(m & strcmp(t.model,'conditional') & strcmp(t.freq,'HF'));
            p = pval(c); q = c.q;
            if q < 0.05, cls = 'direct';
            elseif ~isempty(pw) && pw.p_F < 0.05 && p >= 0.05, cls = 'mediated';
            elseif p < 0.05, cls = 'uncorrected';
            else, cls = 'none'; end
            ok = true;
            if ~isempty(tr) && strcmp(cls,'direct')
                ok = (d == 1 && tr.value > 0) || (d == 2 && tr.value < 0);
            end
            L(end+1) = struct('region', regs{r}, 'band', bands{b}, 'dir', d, 'gc', c.value, 'q', q, 'p', p, ...
                'p_pair', first(pw.p_F, NaN), 'lf', lf, 'hf', hf, 'class', cls, 'trgc_ok', ok); %#ok<AGROW>
        end
    end
end

% Text
lines = {};
nBeats = sum(isfinite(S.heart));
dur = nBeats*mean(S.heart,'omitnan')/1000;
lines{end+1} = sprintf('Data: %d heartbeats (%.1f min of clean EEG and heart signal), mean heart rate %.0f bpm, %d regions x %d bands.', ...
    nBeats, dur/60, 60000/mean(S.heart,'omitnan'), nr, nb);
nd = sum(strcmp({L.class}, 'direct')); nu = sum(strcmp({L.class}, 'uncorrected'));
lines{end+1} = sprintf(['Tests: %d EEG <-> heart links; %d significant after FDR (direct), %d only before correction ' ...
    '(about %.1f would be by chance at p < .05).'], numel(L), nd, nu, 0.05*numel(L));
% direct links, both directions of a pair in one sentence, strongest first
D = L(strcmp({L.class}, 'direct'));
o = [];
if ~isempty(D)
    pairs = unique(strcat({D.region}, '|', {D.band}), 'stable');
    pg = cellfun(@(p) max([D(strcmp(strcat({D.region}, '|', {D.band}), p)).gc]), pairs);
    [~, o] = sort(pg, 'descend');
end
for k = o(:)'
    d = D(strcmp(strcat({D.region}, '|', {D.band}), pairs{k}));
    rb = strrep(pairs{k}, '|', ' ');
    both = L(strcmp({L.region}, d(1).region) & strcmp({L.band}, d(1).band));
    if numel(d) == 2
        [~, i] = max([d.gc]); [~, j] = min([d.gc]);
        s = sprintf(['%s <-> heart: direct in both directions (not explained by respiration or the other bands); ' ...
            '%s stronger (GC %.3f vs %.3f), q <= %.3g, %s.'], rb, dir_word(d(i)), d(i).gc, d(j).gc, max([d.q]), ...
            freq_word(sum([d.lf]), sum([d.hf])));
        dom = d(i);
    else
        s = sprintf('%s: direct (not explained by respiration or the other bands), GC %.3f, q = %.3g, %s.', ...
            link_name(d), d.gc, d.q, freq_word(d.lf, d.hf));
        dom = d;
        if numel(both) == 2 && both(3-dom.dir).gc > dom.gc, dom = []; end   % not the net direction: no TRGC claim
    end
    if ~isempty(dom) && ~dom.trgc_ok
        s = [s ' The net direction does not reverse with time (time-reversed GC): it may come from noise or mixing.']; %#ok<AGROW>
    end
    lines{end+1} = s; %#ok<AGROW>
end
M = L(strcmp({L.class}, 'mediated'));
if ~isempty(M)
    lines{end+1} = ['Pairwise only (the association runs through respiration or the other bands): ' ...
        strjoin(arrayfun(@link_name, M, 'UniformOutput', false), ', ') '.'];
end
U = L(strcmp({L.class}, 'uncorrected'));
if ~isempty(U)
    lines{end+1} = ['Significant before correction only (treat as exploratory): ' ...
        strjoin(arrayfun(@link_name, U, 'UniformOutput', false), ', ') '.'];
end
if nd == 0 && isempty(M)
    lines{end+1} = 'No EEG <-> heart link stands out in this recording.';
end
hiBand = cellfun(@(b) S.bands(strcmp(bands, b), 1) >= 13, {D.band});   % bands starting at 13 Hz or above
if any(hiBand)
    lines{end+1} = ['Beta (and faster) power, at temporal and frontal sites above all, often carries scalp muscle ' ...
        'activity (EMG), which follows arousal and posture along with heart rate: check that the cleaning removed ' ...
        'muscle components before reading these links as cortical.'];
end

% Respiration
if isempty(S.resp)
    lines{end+1} = ['Respiration: none in the model. HF coupling may be the respiratory drive of both the heart ' ...
        'and the EEG: record respiration, or use the ECG- or PPG-derived respiration.'];
else
    coh = out.resp_heart_coh;
    qual = 'tracks breathing well'; if max(coh) < 0.5, qual = 'tracks breathing only partly'; end
    if max(coh) < 0.2, qual = 'does not seem to track breathing: HF results may still hold the respiratory drive'; end
    lines{end+1} = sprintf('Respiration (%s): coherence with the heart period in HF %s, so it %s.', ...
        strrep(S.resp_source, '+', ' + '), strjoin(compose('%.2f', coh), ' / '), qual);
    e = t(strcmp(t.model,'exogenous') & ~strcmp(t.target,'heart') & ismember(t.target, bands), :);
    hit = e(e.p_F < 0.05/height(e), :);
    if ~isempty(hit)
        lines{end+1} = ['Respiration modulates EEG power (Bonferroni): ' ...
            strjoin(strcat(hit.region, {' '}, hit.target), ', ') '.'];
    end
end
if isfield(S,'cfa') && isfield(S.cfa,'share_var') && ~isempty(S.cfa.share_var)
    lines{end+1} = sprintf(['Cardiac field: removed by ''%s'' (median %.1f%% of a channel''s variance, max %.1f%%) before ' ...
        'the band powers.'], S.cfa.method, 100*median(S.cfa.share_var), 100*max(S.cfa.share_var));
end

% Model checks
lb = vertcat(out.models.ljungbox_p);
bad = any(lb*size(lb,2) < 0.01, 2);
spc = [out.models.samples_per_coef];
if any(bad)
    lines{end+1} = sprintf('Model check: residuals not white in %s (Ljung-Box): those regions'' results are less reliable.', ...
        strjoin(regs(bad), ', '));
else
    rng_s = sprintf('%.0f-%.0f', min(spc), max(spc));
    if round(min(spc)) == round(max(spc)), rng_s = sprintf('%.0f', min(spc)); end
    lines{end+1} = sprintf('Model check: residuals white in every region; %s samples per coefficient.', rng_s);
end
if ~isempty(out.warnings)
    lines{end+1} = sprintf('%d warning(s): see out.warnings.', numel(out.warnings));
end
lines{end+1} = ['Read with care: a direction is temporal precedence (Granger), not proof of causation; an unrecorded ' ...
    'source (e.g., brainstem arousal) reaching the heart and the cortex at different delays makes one appear to lead. ' ...
    'GC values compare across recordings with the same settings; DC and GPDC only within the same model.'];
fprintf('\nBrain-heart coupling report\n');
fprintf('  - %s\n', lines{:});

if ~show, return; end

% Figure
f = figure('Color','w','Name','Brain-heart coupling report','NumberTitle','off','Position',[60 40 1250 940]);
annotation(f, 'textbox', [0.02 0.93 0.96 0.06], 'String', {'Brain-heart coupling report', lines{1}}, ...
    'EdgeColor','none', 'FontSize', 11, 'Interpreter','none', 'FontWeight','bold');
cLF = [0.15 0.45 0.80]; cHF = [0.90 0.45 0.10]; cBB = [0.35 0.25 0.55];
gmax = max([L.gc 0.01]);
if nr <= 6
    for r = 1:nr
        ax = axes(f, 'Position', [0.03 + (r-1)*0.95/nr 0.60 0.95/nr - 0.02 0.29]); hold(ax,'on'); axis(ax,'off');
        xb = linspace(1, nb, nb); xh = (1 + nb)/2;
        xlim(ax, [0.4 nb+0.6]); ylim(ax, [-0.35 1.35]);
        pos = getpixelposition(ax);                    % data -> pixel scale, for undistorted arrowheads
        sc = [pos(3)/diff(xlim(ax)) pos(4)/diff(ylim(ax))];
        for k = find(strcmp({L.region}, regs{r}))
            if strcmp(L(k).class, 'none') || strcmp(L(k).class, 'uncorrected'), continue; end
            b = find(strcmp(bands, L(k).band));
            off = 0.24*(2*(L(k).dir == 1) - 1);       % the two directions side by side
            if L(k).dir == 1, p1 = [xb(b) 1]; p2 = [xh 0]; else, p1 = [xh 0]; p2 = [xb(b) 1]; end
            if strcmp(L(k).class, 'direct')
                col = cBB; if L(k).lf > 1.5*L(k).hf, col = cLF; elseif L(k).hf > 1.5*L(k).lf, col = cHF; end
                draw_arrow(ax, p1, p2, off, 1 + 7*L(k).gc/gmax, col, '-', sc);
            else
                draw_arrow(ax, p1, p2, off, 1.2, [0.6 0.6 0.6], '--', sc);
            end
        end
        scatter(ax, xb, ones(1,nb), 900, [0.92 0.92 0.92], 'filled', 'MarkerEdgeColor', [0.4 0.4 0.4]);
        text(ax, xb, ones(1,nb), bands, 'HorizontalAlignment','center', 'FontSize', 8);
        scatter(ax, xh, 0, 1300, [1 0.85 0.85], 'filled', 'MarkerEdgeColor', [0.6 0.2 0.2]);
        text(ax, xh, 0, 'heart', 'HorizontalAlignment','center', 'FontSize', 9, 'FontWeight','bold');
        title(ax, regs{r}, 'FontSize', 10);
    end
    ax = axes(f, 'Position', [0.03 0.565 0.94 0.03]); axis(ax,'off'); hold(ax,'on');
    plot(ax, [0 0.04], [0.5 0.5], '-', 'Color', cLF, 'LineWidth', 3); text(ax, 0.045, 0.5, 'direct, LF', 'FontSize', 8);
    plot(ax, [0.17 0.21], [0.5 0.5], '-', 'Color', cHF, 'LineWidth', 3); text(ax, 0.215, 0.5, 'direct, HF', 'FontSize', 8);
    plot(ax, [0.34 0.38], [0.5 0.5], '-', 'Color', cBB, 'LineWidth', 3); text(ax, 0.385, 0.5, 'direct, both', 'FontSize', 8);
    plot(ax, [0.53 0.57], [0.5 0.5], '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.2); text(ax, 0.575, 0.5, 'pairwise only (mediated)', 'FontSize', 8);
    text(ax, 0.78, 0.5, 'arrowhead = direction; width = GC', 'FontSize', 8);
    xlim(ax, [0 1]); ylim(ax, [0 1]);
else
    plot_coupling(out);                                % many regions (channels): the maps instead
end

% Strongest links: GC in LF and HF
ax = axes(f, 'Position', [0.17 0.33 0.32 0.20]);
[~, o] = sort([L.gc], 'descend'); top = L(o(1:min(8, numel(L))));
if numel(top) < 2, top = L(1:min(2, numel(L))); end
bh = barh(ax, [[top.lf]' [top.hf]'], 'stacked'); bh(1).FaceColor = cLF; bh(2).FaceColor = cHF;
set(ax, 'YTick', 1:numel(top), 'YTickLabel', arrayfun(@link_name, top, 'UniformOutput', false), 'YDir','reverse', 'FontSize', 8);
for k = 1:numel(top)
    if strcmp(top(k).class,'direct'), text(ax, top(k).lf + top(k).hf, k, '  *', 'FontSize', 12, 'FontWeight','bold'); end
end
xlabel(ax, 'spectral GC: LF (blue) + HF (orange)');
title(ax, 'Strongest links (* = direct, after FDR)', 'FontSize', 9); box(ax, 'off');

% Checks
ax = axes(f, 'Position', [0.56 0.33 0.41 0.20]); axis(ax,'off');
chk = {};
if isempty(S.resp), chk(end+1,:) = {'Respiration', 'none', 'r'};
else
    v = max(out.resp_heart_coh); c = 'g'; if v < 0.5, c = 'y'; end; if v < 0.2, c = 'r'; end
    chk(end+1,:) = {'Respiration tracks breathing', sprintf('HF coherence %.2f', v), c};
end
chk(end+1,:) = {'Residuals white', sprintf('%d of %d regions', sum(~bad), nr), ifelse(any(bad), 'y', 'g')};
chk(end+1,:) = {'Samples per coefficient', sprintf('%.0f (>= 10)', min(spc)), ifelse(min(spc) < 10, 'r', 'g')};
if isfield(S.cfa,'share_var') && ~isempty(S.cfa.share_var)
    chk(end+1,:) = {'Cardiac field removed', sprintf('%s, median %.1f%%', S.cfa.method, 100*median(S.cfa.share_var)), 'g'};
else
    chk(end+1,:) = {'Cardiac field removed', 'no', 'y'};
end
chk(end+1,:) = {'Chance level', sprintf('%d tests: ~%.1f expected at p < .05', numel(L), 0.05*numel(L)), 'n'};
cmap = struct('g', [0.2 0.65 0.3], 'y', [0.95 0.7 0.1], 'r', [0.85 0.2 0.2], 'n', [0.6 0.6 0.6]);
for k = 1:size(chk,1)
    y = 1 - (k-0.5)/size(chk,1);
    text(ax, 0.02, y, char(9679), 'Color', cmap.(chk{k,3}), 'FontSize', 14);
    text(ax, 0.07, y, chk{k,1}, 'FontSize', 9, 'FontWeight','bold');
    text(ax, 0.55, y, chk{k,2}, 'FontSize', 9);
end
title(ax, 'Quality checks', 'FontSize', 9);

% Interpretation
annotation(f, 'textbox', [0.03 0.01 0.94 0.25], 'String', strcat([char(8226) ' '], lines(2:end)), ...
    'EdgeColor', 'none', 'FontSize', 9, 'Interpreter','none', 'VerticalAlignment','top');
finish_figure(f)


function s = link_name(l)
if l.dir == 1, s = sprintf('%s %s -> heart', l.region, l.band); else, s = sprintf('heart -> %s %s', l.region, l.band); end


function s = dir_word(l)
if l.dir == 1, s = sprintf('%s -> heart', l.band); else, s = sprintf('heart -> %s', l.band); end


function w = freq_word(lf, hf)
if lf > 1.5*hf, w = 'mostly in LF (0.04-0.15 Hz)';
elseif hf > 1.5*lf, w = 'mostly in HF (0.15-0.4 Hz)';
else, w = 'in both LF and HF'; end


function draw_arrow(ax, p1, p2, off, lw, col, ls, sc)
% Arrow from p1 to p2 (data units), drawn in pixel space so that the head
% keeps its shape (sc: pixels per data unit); off shifts it sideways so
% that the two directions of a pair sit side by side
P1 = p1.*sc; P2 = p2.*sc;                              % to pixels
d = P2 - P1; u = d/norm(d); n = [-u(2) u(1)];
if u(2) > 0, n = -n; end                               % same side reference for both directions
off = off*sc(1);
P1 = P1 + off*n + 26*u; P2 = P2 + off*n - 30*u;        % stop short of the nodes
h = 9 + lw; wdt = 4 + 0.6*lw;
base = P2 - h*u;
toD = @(P) P./sc;
a = toD(P1); b = toD(base);
plot(ax, [a(1) b(1)], [a(2) b(2)], ls, 'Color', col, 'LineWidth', lw);
tri = [P2; base + wdt*n; base - wdt*n];
tri = tri./sc;
patch(ax, tri(:,1), tri(:,2), col, 'EdgeColor', col);


function v = first(a, b)
if isempty(a) || isnan(a(1)), v = b; else, v = a(1); end
if isempty(v), v = NaN; end


function o = ifelse(c, a, b)
if c, o = a; else, o = b; end
