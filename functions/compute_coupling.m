% COMPUTE_COUPLING - Brain-heart coupling between EEG band power and the
% heart period (NN), with respiration as an exogenous input.
%
% For each region, one model holds the power envelopes of the EEG bands and
% the NN series (from GET_COUPLING_SERIES), with respiration as an
% exogenous input (VARX: it drives both but is not modelled itself; a
% smooth, quasi-periodic breathing signal is not an autoregressive process,
% and modelling it as one blew up the model order: Baselli et al., 1997;
% Porta et al., 2000). From that model:
%   - strength and direction: conditional Granger causality (GC) from each
%     EEG band to the heart and back, given the other bands and
%     respiration, in the time domain and in the LF (0.04-0.15 Hz) and HF
%     (0.15-0.4 Hz) bands of the fluctuations (state space, VAR_SSGC);
%   - direct vs mediated: the same GC from a model of the band and the heart
%     alone (pairwise). A link that is pairwise but not conditional runs
%     through respiration or the other bands;
%   - coherence (total, pairwise model) and partial coherence (direct, given
%     the other bands and respiration); directed coherence (total,
%     normalized by inflow) and generalized partial directed coherence
%     (direct, normalized by outflow) (VAR_SPECTRAL);
%   - a rank-based (Gaussian copula) conditional transfer entropy, robust
%     to non-Gaussian amplitudes and outliers (Ince et al., 2017);
%   - time-reversed GC, which flags a direction that may come from noise
%     or mixing rather than from causality (Haufe et al., 2013; Winkler et
%     al., 2016);
%   - respiration -> heart and respiration -> band (exogenous GC), and the
%     zero-lag residual partial correlations (coupling faster than one grid
%     sample, which no directed measure sees).
%
% Usage:
%   out = compute_coupling(S, opt)
%
% INPUTS:
%   S   - series from GET_COUPLING_SERIES
%   opt - struct (all optional):
%     .window - analysis window (s; default 300). Longer recordings are cut
%               into windows, each fitted and tested on its own (F-tests);
%               the summary model is fitted on all windows together.
%     .maxlag - longest cross lag (one variable on another) considered (s;
%               default 5)
%     .maxlag_own - longest own lag considered (s; default 10). Orders are
%               chosen by BIC per equation, own and cross lags separately
%               (VAR_ORDER), with >= 10 samples per coefficient.
%     .maxlag_white - own lags are then added, one at a time, to each
%               equation whose residuals fail Ljung-Box (p < .05), up to
%               this lag (s; default 30; 0 = BIC orders only). BIC left
%               the heart-period residuals autocorrelated in 40-58% of real
%               recordings (ds003838, HEP-NF; own lags predicting ~1.5% of
%               the residual variance); the added lags whiten them without
%               changing the EEG -> heart F-tests (null 4.6-5.3% at
%               p < .05, vs 4.7-5.8%).
%     .nfreq  - frequency intervals from 0 to Nyquist (default 128)
%     .nboot  - parametric bootstrap samples for the conditional GC tests
%               (time and LF/HF; default 200; 0 = F-tests only). The null
%               model is the fitted model with the tested link removed and
%               the target equation refitted, simulated with the recorded
%               respiration: the null keeps every other coupling and the
%               autocorrelations (Faes et al., 2010).
%     .seed   - random seed (default 0)
%
% OUTPUT (out fields):
%   table   - one row per region, band, direction and measure: region,
%             band, source, target, measure (GC, TE, TRGC, COH, PCOH, DC,
%             GPDC), model ('conditional', 'pairwise', 'undirected',
%             'exogenous'), freq ('all', 'LF', 'HF'), value, p (bootstrap
%             for conditional GC, F-test for pairwise and exogenous GC,
%             circular shifts for TE; NaN where no test), null_median (bias
%             reference: GC under the null model), p_F (F-test, with a
%             heteroskedasticity-robust covariance (HC3) and the degrees of
%             freedom reduced for the per-run detrending), q (FDR,
%             Benjamini-Hochberg, across the out.n_tests conditional
%             EEG <-> heart GC tests of the recording)
%   freqs   - frequency axis of the spectra (Hz)
%   spectra - nRegions x nBands x nFreqs arrays: gc_eh / gc_he (conditional,
%             EEG -> heart / heart -> EEG), gc_eh_pair / gc_he_pair, coh,
%             pcoh, dc_eh / dc_he, gpdc_eh / gpdc_he, and gc_eh_null95 /
%             gc_he_null95 (95th percentile of the bootstrap null)
%   models  - per region: variables, own orders (order_own_bic: before
%             the lags added for white residuals), cross orders (BIC), the
%             lag-order matrix P, nObs, samples per coefficient (largest
%             equation), spectral radius, Ljung-Box p of each residual
%             series, zero-lag partial correlations
%   windows - per window: start time, and conditional/pairwise GC with
%             F-test p (nRegions x nBands x 2 directions x nWindows)
%   resp_heart_coh - peak HF coherence between respiration and the heart
%             period (first region's series): does the respiration signal
%             track breathing? (warning below 0.2)
%   grid, fs - 'beat' or 'uniform', and samples per second (beat grid: at
%             the mean heart period; frequencies in Hz assume it)
%   warnings - cell array of text warnings
%
% GC values are natural-log ratios of prediction error variances (TE in
% bits is GC/(2 ln 2) for Gaussian data); they carry a small positive
% bias, given by null_median.
%
% REFERENCES:
%   Barnett, L. & Seth, A.K. (2015). Phys Rev E 91, 040101.
%   Baselli, G. et al. (1997). Med Biol Eng Comput 35, 488-496.
%   Faes, L., Porta, A. & Nollo, G. (2010). IEEE TBME 57, 1897-1906.
%   Faes, L., Erla, S. & Nollo, G. (2012). Comput Math Methods Med, 140513.
%   Ince, R.A.A. et al. (2017). Hum Brain Mapp 38, 1541-1573.
%   MacKinnon, J.G. & White, H. (1985). J Econometrics 29, 305-325.
%   Winkler, I. et al. (2016). IEEE Trans Signal Process 64, 3110-3121.
%
% Copyright (C) - Cedric Cannard, BrainBeats 2026

function out = compute_coupling(S, opt)

if nargin < 2, opt = struct; end
opt = set_default(opt, 'window', 300);
opt = set_default(opt, 'maxlag', 5);
opt = set_default(opt, 'maxlag_own', 10);
opt = set_default(opt, 'maxlag_white', 30);
opt = set_default(opt, 'nfreq', 128);
opt = set_default(opt, 'nboot', 200);
opt = set_default(opt, 'seed', 0);
fso = S.fs_out;
nr = numel(S.regions); nb = size(S.bands,1);
hasResp = ~isempty(S.resp);
lf = [0.04 0.15]; hf = [0.15 0.4];
freqs = (0:opt.nfreq)/opt.nfreq * fso/2;
iLF = freqs >= lf(1) & freqs < lf(2);
iHF = freqs >= hf(1) & freqs <= hf(2);
out.warnings = {};
out.resp_heart_coh = NaN;
if ~hasResp
    out.warnings{end+1} = ['No respiration: HF coupling may be the respiratory drive of both heart ' ...
        'rate and EEG; conditional GC conditions on the other EEG bands only.'];
end

% Windows (grid indices)
T = numel(S.t);
wlen = round(opt.window*fso);
if T < 1.5*wlen
    win = [1 T];
else
    st = 1:wlen:T-wlen+1;
    win = [st' st'+wlen-1];
    win(end,2) = T;
end
nw = size(win,1);

vasc = []; vascN = {};
if isfield(S,'vasc') && ~isempty(S.vasc), vasc = S.vasc; vascN = S.vasc_names; end
nvs = numel(vascN);
vnames = [S.band_names(:)' {'heart'} vascN];
nv = nb + 1 + nvs; ih = nb + 1; iv = ih + (1:nvs); nu = size(S.resp, 1);
respN = 'resp';

rows = {}; tabcols = {'region','band','source','target','measure','model','freq','value','p','null_median','p_F'};
sz = [nr nb opt.nfreq+1];
sp = struct('gc_eh',nan(sz),'gc_he',nan(sz),'gc_eh_pair',nan(sz),'gc_he_pair',nan(sz), ...
    'coh',nan(sz),'pcoh',nan(sz),'dc_eh',nan(sz),'dc_he',nan(sz),'gpdc_eh',nan(sz),'gpdc_he',nan(sz), ...
    'gc_eh_null95',nan(sz),'gc_he_null95',nan(sz));
wres = struct('start', S.t(win(:,1)), 'gc_cond', nan(nr,nb,2,nw), 'p_cond', nan(nr,nb,2,nw), ...
    'gc_pair', nan(nr,nb,2,nw), 'p_pair', nan(nr,nb,2,nw));
rng(opt.seed);

for r = 1:nr
    rname = S.regions(r).name;
    X = [reshape(S.eeg(r,:,:), nb, []); S.heart; vasc];
    if hasResp, X = [X; S.resp]; end %#ok<AGROW>

    % Detrend and standardize each window; windows separated by NaN
    Xw = nan(size(X,1), T + nw);
    pos = 0; wmask = zeros(1, T + nw);
    for w = 1:nw
        seg = detrend_std(X(:, win(w,1):win(w,2)));
        L = size(seg,2);
        Xw(:, pos+(1:L)) = seg; wmask(pos+(1:L)) = w;
        pos = pos + L + 1;                         % one NaN between windows
    end
    Xw = Xw(:, 1:pos-1); wmask = wmask(1:pos-1);
    Uw = Xw(nv+1:end, :); Xw = Xw(1:nv, :);       % exogenous respiration

    % Orders (BIC) per equation, own lags and cross lags separately: EEG
    % envelopes need 1-2 lags, the smooth (interpolated) NN series long own
    % lags but not long cross lags. Capped at 10 samples per coefficient.
    nValid = sum(all(isfinite([Xw; Uw]),1));
    maxcoef = floor(nValid/10);
    pmax = round(opt.maxlag*fso); pown = round(opt.maxlag_own*fso);
    if maxcoef < nv + 2*nu
        out.warnings{end+1} = sprintf('%s: too few valid samples (%d) for a %d-variable model.', rname, nValid, nv);
        continue
    end
    [~, ~, crit] = var_order(Xw, pmax, Uw, pown, maxcoef);
    qx = crit.qx; pwhite = round(opt.maxlag_white*fso);
    [P, A, V, info, lbp] = fit_white(Xw, Uw, crit.P, qx, pwhite, maxcoef);
    lbmin = min(lbp);
    exoN = {}; if nu > 0, exoN = strsplit(S.resp_source, '+'); end
    mdl = struct('name', rname, 'variables', {vnames}, 'exogenous', {exoN}, ...
        'order_own', diag(P)', 'order_own_bic', diag(crit.P)', 'order_cross', max(P - diag(diag(P)), [], 2)', 'order_resp', qx', 'P', P, 'nObs', info.nObs, ...
        'samples_per_coef', info.nObs/max(info.ncoef), 'rho', info.rho, ...
        'ljungbox_minp', lbmin, 'ljungbox_p', lbp', 'zerolag_pcorr', partial_corr(V));
    if any(info.ncoef >= maxcoef - nv)
        out.warnings{end+1} = sprintf('%s: model order limited by the number of samples (10 per coefficient).', rname);
    end
    if info.rho >= 1
        out.warnings{end+1} = sprintf('%s: the fitted model is unstable (spectral radius %.3f).', rname, info.rho);
    end
    if lbmin*nv < 0.01                                % Bonferroni over the variables
        out.warnings{end+1} = sprintf(['%s: residuals are not white (Ljung-Box p = %.2g, %s) even with %g s of ' ...
            'own lags: the model may miss slow structure.'], rname, lbmin, vnames{find(lbp == lbmin, 1)}, opt.maxlag_white);
    end

    spec = var_spectral(A, V, opt.nfreq);
    [Arev, Vrev] = var_fit(fliplr(Xw), P, [], fliplr(Uw), qx);    % reversed time, for TRGC
    Xc = copnorm_rows(Xw); Uc = copnorm_rows(Uw);                 % for the transfer entropy

    for b = 1:nb
        bname = S.band_names{b};
        others = setdiff(1:nv, [b ih]);

        % pairwise model: the band and the heart alone (no respiration)
        X2 = Xw([b ih], :);
        mc2 = floor(sum(all(isfinite(X2),1))/10);
        [~, ~, c2] = var_order(X2, pmax, [], pown, mc2);
        [P2, A2, V2, info2] = fit_white(X2, [], c2.P, [], pwhite, mc2);
        spec2 = var_spectral(A2, V2, opt.nfreq);

        for dir = 1:2
            if dir == 1, src = b; tgt = ih; srcN = bname; tgtN = 'heart'; s2 = 1; t2 = 2;
            else,        src = ih; tgt = b; srcN = 'heart'; tgtN = bname; s2 = 2; t2 = 1; end

            [Fc, fc] = var_ssgc(A, V, tgt, src, others, opt.nfreq);
            [Fp, fp] = var_ssgc(A2, V2, t2, s2, [], opt.nfreq);
            pFc = ftest(Xw, Uw, info.rows, P, qx, tgt, src);
            pFp = ftest(X2, [], info2.rows, P2, [], t2, s2);

            [pb, nullmed, null95] = deal(nan(1,3), nan(1,3), nan(1,opt.nfreq+1));
            if opt.nboot > 0
                [pb, nullmed, null95] = boot_gc(Xw, Uw, info, A, tgt, src, others, ...
                    [Fc mean(fc(iLF)) mean(fc(iHF))], opt.nboot, opt.nfreq, iLF, iHF);
            end
            [te, pte, temed] = gcmi_te(Xc, Uc, info.rows, P, qx, tgt, src, max(opt.nboot, 100), round(20*fso));

            rows(end+1,:) = {rname, bname, srcN, tgtN, 'GC', 'conditional', 'all', Fc, pb(1), nullmed(1), pFc}; %#ok<AGROW>
            rows(end+1,:) = {rname, bname, srcN, tgtN, 'GC', 'conditional', 'LF', mean(fc(iLF)), pb(2), nullmed(2), NaN}; %#ok<AGROW>
            rows(end+1,:) = {rname, bname, srcN, tgtN, 'GC', 'conditional', 'HF', mean(fc(iHF)), pb(3), nullmed(3), NaN}; %#ok<AGROW>
            rows(end+1,:) = {rname, bname, srcN, tgtN, 'GC', 'pairwise', 'all', Fp, pFp, NaN, pFp}; %#ok<AGROW>
            rows(end+1,:) = {rname, bname, srcN, tgtN, 'GC', 'pairwise', 'LF', mean(fp(iLF)), NaN, NaN, NaN}; %#ok<AGROW>
            rows(end+1,:) = {rname, bname, srcN, tgtN, 'GC', 'pairwise', 'HF', mean(fp(iHF)), NaN, NaN, NaN}; %#ok<AGROW>
            rows(end+1,:) = {rname, bname, srcN, tgtN, 'TE', 'conditional', 'all', te, pte, temed, NaN}; %#ok<AGROW>
            for fb = {'LF' iLF; 'HF' iHF}'
                rows(end+1,:) = {rname, bname, srcN, tgtN, 'DC', 'conditional', fb{1}, mean(spec.dc(tgt,src,fb{2})), NaN, NaN, NaN}; %#ok<AGROW>
                rows(end+1,:) = {rname, bname, srcN, tgtN, 'GPDC', 'conditional', fb{1}, mean(spec.gpdc(tgt,src,fb{2})), NaN, NaN, NaN}; %#ok<AGROW>
            end

            if dir == 1
                sp.gc_eh(r,b,:) = fc; sp.gc_eh_pair(r,b,:) = fp; sp.gc_eh_null95(r,b,:) = null95;
                sp.dc_eh(r,b,:) = spec.dc(tgt,src,:); sp.gpdc_eh(r,b,:) = spec.gpdc(tgt,src,:);
                Feh = Fc; sigEH = min(pb(1), pFc) < 0.05;
            else
                sp.gc_he(r,b,:) = fc; sp.gc_he_pair(r,b,:) = fp; sp.gc_he_null95(r,b,:) = null95;
                sp.dc_he(r,b,:) = spec.dc(tgt,src,:); sp.gpdc_he(r,b,:) = spec.gpdc(tgt,src,:);
                Fhe = Fc; sigHE = min(pb(1), pFc) < 0.05;
            end
        end

        % time-reversed GC: net flow EEG -> heart minus the same on reversed
        % time; flagged when a significant link's direction does not reverse
        net = Feh - Fhe;
        netRev = var_ssgc(Arev, Vrev, ih, b, others) - var_ssgc(Arev, Vrev, b, ih, others);
        rows(end+1,:) = {rname, bname, bname, 'heart', 'TRGC', 'conditional', 'all', net - netRev, NaN, NaN, NaN}; %#ok<AGROW>
        if ((net > 0 && sigEH) || (net < 0 && sigHE)) && sign(net) == sign(netRev)
            out.warnings{end+1} = sprintf(['%s %s: the dominant direction does not reverse with time ' ...
                '(time-reversed GC): it may come from noise or mixing rather than causality.'], rname, bname);
        end

        % undirected: total coherence (pairwise model), partial coherence (full model)
        sp.coh(r,b,:) = spec2.coh(1,2,:); sp.pcoh(r,b,:) = spec.pcoh(b,ih,:);
        for fb = {'LF' iLF; 'HF' iHF}'
            rows(end+1,:) = {rname, bname, bname, 'heart', 'COH', 'undirected', fb{1}, mean(spec2.coh(1,2,fb{2})), NaN, NaN, NaN}; %#ok<AGROW>
            rows(end+1,:) = {rname, bname, bname, 'heart', 'PCOH', 'undirected', fb{1}, mean(spec.pcoh(b,ih,fb{2})), NaN, NaN, NaN}; %#ok<AGROW>
        end

        % respiration -> band
        if hasResp
            [g, pr] = exog_gc(Xw, Uw, info.rows, P, qx, b);
            rows(end+1,:) = {rname, bname, 'resp', bname, 'GC', 'exogenous', 'all', g, pr, NaN, pr}; %#ok<AGROW>
        end

        % per window: F-tests only
        if nw > 1
            for w = 1:nw
                Xk = Xw; Xk(:, wmask ~= w) = NaN;
                [Ak, Vk, ik] = var_fit(Xk, P, [], Uw, qx);
                X2k = Xk([b ih], :);
                [A2k, V2k, i2k] = var_fit(X2k, P2);
                for dir = 1:2
                    if dir == 1, src = b; tgt = ih; s2 = 1; t2 = 2; else, src = ih; tgt = b; s2 = 2; t2 = 1; end
                    wres.gc_cond(r,b,dir,w) = var_ssgc(Ak, Vk, tgt, src, others);
                    wres.gc_pair(r,b,dir,w) = var_ssgc(A2k, V2k, t2, s2, []);
                    wres.p_cond(r,b,dir,w) = ftest(Xk, Uw, ik.rows, P, qx, tgt, src);
                    wres.p_pair(r,b,dir,w) = ftest(X2k, [], i2k.rows, P2, [], t2, s2);
                end
            end
        end
    end

    % respiration -> heart (respiratory sinus arrhythmia), and whether the
    % respiration signal tracks breathing: its coherence with the heart
    % period in HF (peak). Granger causality alone is no check: the heart's
    % own past already predicts a regular respiratory oscillation.
    if hasResp
        [g, pr] = exog_gc(Xw, Uw, info.rows, P, qx, ih);
        rows(end+1,:) = {rname, '', respN, 'heart', 'GC', 'exogenous', 'all', g, pr, NaN, pr}; %#ok<AGROW>
        if r == 1                                     % each proxy on its own
            cohR = nan(1, nu);
            for u = 1:nu
                Xr = [Xw(ih,:); Uw(u,:)];
                [~, ~, c3] = var_order(Xr, pmax, [], pown, floor(sum(all(isfinite(Xr),1))/10));
                [A3, V3] = var_fit(Xr, c3.P);
                sp3 = var_spectral(A3, V3, opt.nfreq);
                cohR(u) = max(sp3.coh(1,2,iHF));
            end
            out.resp_heart_coh = cohR;
            if max(cohR) < 0.2
                out.warnings{end+1} = sprintf(['Respiration (%s) shares little with the heart period in HF ' ...
                    '(peak coherence %.2f): it may not track breathing, and conditioning on it may not remove ' ...
                    'the respiratory drive.'], S.resp_source, max(cohR));
            end
        end
    end

    % Vascular series (PPG): conditional GC with each band and with the
    % heart period (F-tests; e.g. pulse arrival time -> heart period is
    % baroreflex-like, EEG -> pulse amplitude a vasomotor path)
    for vi = 1:nvs
        k = iv(vi);
        pairs = [num2cell(1:nb) {ih}];
        for jj = 1:numel(pairs)
            j = pairs{jj};
            for dir = 1:2
                if dir == 1, src = j; tgt = k; else, src = k; tgt = j; end
                oth = setdiff(1:nv, [src tgt]);
                [Fv, fv] = var_ssgc(A, V, tgt, src, oth, opt.nfreq);
                pF = ftest(Xw, Uw, info.rows, P, qx, tgt, src);
                bnd = ''; if j <= nb, bnd = vnames{j}; end
                rows(end+1,:) = {rname, bnd, vnames{src}, vnames{tgt}, 'GC', 'conditional', 'all', Fv, pF, NaN, pF}; %#ok<AGROW>
                rows(end+1,:) = {rname, bnd, vnames{src}, vnames{tgt}, 'GC', 'conditional', 'LF', mean(fv(iLF)), NaN, NaN, NaN}; %#ok<AGROW>
                rows(end+1,:) = {rname, bnd, vnames{src}, vnames{tgt}, 'GC', 'conditional', 'HF', mean(fv(iHF)), NaN, NaN, NaN}; %#ok<AGROW>
            end
        end
        if hasResp
            [g, pr] = exog_gc(Xw, Uw, info.rows, P, qx, k);
            rows(end+1,:) = {rname, '', respN, vnames{k}, 'GC', 'exogenous', 'all', g, pr, NaN, pr}; %#ok<AGROW>
        end
    end
    out.models(r) = mdl;
end

if isempty(rows)
    error('compute_coupling: no region could be modelled: %s', strjoin(out.warnings, ' '));
end
out.table = cell2table(rows, 'VariableNames', tabcols);
% FDR (Benjamini-Hochberg) across the conditional EEG <-> heart GC tests of
% the recording (all regions, bands, directions), on the F-test p-values:
% the bootstrap's resolution, 1/(nboot+1), is too coarse for ~32 tests
% (the smallest q would need p < .05/32)
t = out.table;
sel = strcmp(t.measure,'GC') & strcmp(t.model,'conditional') & strcmp(t.freq,'all') & ...
    (strcmp(t.source,'heart') | strcmp(t.target,'heart')) & ~ismember(t.source, vascN) & ~ismember(t.target, vascN);
t.q = nan(height(t), 1);
t.q(sel) = fdr_bh(t.p_F(sel));
out.table = t;
out.n_tests = sum(sel);
if isfield(S, 'grid'), out.grid = S.grid; else, out.grid = 'uniform'; end
out.fs = fso;
out.freqs = freqs;
out.spectra = sp;
out.windows = wres;
out.lf = lf; out.hf = hf;
out.opt = opt;
for k = 1:numel(out.warnings), warning('BrainBeats:coupling', '%s', out.warnings{k}); end


function [pb, med, null95] = boot_gc(X, U, info, A, tgt, src, others, obs, nboot, nfreq, iLF, iHF)
% Parametric bootstrap under "no direct src -> tgt": the target equation is
% refitted without the source lags; data are simulated from that model,
% driven by the recorded respiration, with the observed gaps, and refitted
n = size(X,1); nu = size(U,1);
pm = max(info.P(:));
rows = info.rows;
[Z, cols] = var_design(X, U, rows, info.P, info.qx);
m = cols(tgt,:) & ~src_cols(n, pm, size(Z,1), src);
C = [A info.B];
C(tgt,:) = 0; C(tgt, m) = (Z(m,:)' \ X(tgt, rows)')';
E = X(:, rows) - C*Z;
V0 = (E*E')/numel(rows);
A0 = C(:, 1:n*pm); B0 = C(:, n*pm+1:end);
mask = isfinite(X);
null = zeros(nboot, 3); nullf = zeros(nboot, nfreq+1);
if nu > 0
    Xs = var_sim(A0, V0, size(X,2), nboot, [], B0, U);
else
    Xs = var_sim(A0, V0, size(X,2), nboot);
end
for k = 1:nboot
    x = Xs(:,:,k); x(~mask) = NaN;
    [Ab, Vb] = var_fit(x, info.P, rows, U, info.qx);
    [F, f] = var_ssgc(Ab, Vb, tgt, src, others, nfreq);
    null(k,:) = [F mean(f(iLF)) mean(f(iHF))];
    nullf(k,:) = f;
end
pb = (1 + sum(null >= obs, 1)) / (nboot + 1);
med = median(null, 1);
null95 = prctile_cols(nullf, 95);


function [P, A, V, info, lbp] = fit_white(X, U, P, qx, pcap, maxcoef)
% Fit the model, then add one own lag at a time to each equation whose
% residuals fail Ljung-Box (p < .05), until they pass, reach pcap lags or
% the coefficient cap, or too few time points remain
[A, V, info] = var_fit(X, P, [], U, qx);
[~, lbp] = ljungbox(info.E, info.rows, max(10, 2*max(P(:))), max(P,[],2));
while true
    grow = find(lbp(:) < .05 & diag(P) < pcap & info.ncoef(:) < maxcoef);
    if isempty(grow), break; end
    P1 = P;
    d = sub2ind(size(P), grow, grow);
    P1(d) = P1(d) + 1;
    try
        [A1, V1, info1] = var_fit(X, P1, [], U, qx);
    catch
        break
    end
    P = P1; A = A1; V = V1; info = info1;
    [~, lbp] = ljungbox(info.E, info.rows, max(10, 2*max(P(:))), max(P,[],2));
end


function pval = ftest(X, U, rows, P, qx, tgt, src)
% Granger F-test: the target equation (its own lags, VAR_DESIGN) vs the
% same without the source lags
n = size(X,1);
[Z, cols] = var_design(X, U, rows, P, qx);
m = cols(tgt,:);
drop = src_cols(n, max(P(:)), size(Z,1), src);
pval = nested_f(X(tgt, rows)', Z(m,:)', drop(m), nruns([X; U], rows));


function [g, pval] = exog_gc(X, U, rows, P, qx, tgt)
% Respiration -> target: log ratio of the target's residual variances
% without and with respiration, and its F-test
n = size(X,1);
[Z, cols] = var_design(X, U, rows, P, qx);
m = cols(tgt,:);
drop = false(1, size(Z,1)); drop(n*max(P(:))+1:end) = true;
[pval, RSSu, RSSr] = nested_f(X(tgt, rows)', Z(m,:)', drop(m), nruns([X; U], rows));
g = log(RSSr/RSSu);


function nr = nruns(X, rows)
% Runs of complete samples holding fitted time points: DETREND_STD removed
% the mean and slope of each
ok = all(isfinite(X), 1);
id = cumsum(diff([false ok]) == 1);
nr = numel(unique(id(rows)));


function s = src_cols(n, pm, nz, src)
% Regressors holding the lags of the source variable(s)
s = false(1, nz);
for k = 1:pm, s((k-1)*n + src) = true; end


function [pval, RSSu, RSSr] = nested_f(y, Z, drop, nr)
% Wald F-test of the dropped regressors with a heteroskedasticity-robust
% covariance (HC3: MacKinnon & White, 1985) and the residual degrees of
% freedom reduced by the mean and slope removed from each of the nr runs.
% Real EEG power and heart period residuals are heteroskedastic: with the
% classical F the null datasets (segment and subject swaps) rejected 5.7-7.3%
% at p < .05, with this test 5.0-6.5% (simulation 5.4 -> 5.3%), power unchanged.
[Q, R] = qr(Z, 0);
b = R \ (Q'*y);
ru = y - Z*b;
Zr = Z(:, ~drop);
rr = y - Zr*(Zr \ y);
RSSu = ru'*ru; RSSr = rr'*rr;
h = sum(Q.^2, 2);                                  % leverage
G = R \ (Q' .* (ru ./ max(1 - h, eps))');
C = G*G';                                          % HC3 covariance of b
d1 = sum(drop); d0 = numel(y) - size(Z,2); d2 = max(d0 - 2*nr, 1);
F = b(drop)' * (C(drop,drop) \ b(drop)) / d1 * d2/d0;
pval = betainc(d2/(d2 + d1*max(F,0)), d2/2, d1/2);   % upper tail of F(d1,d2)


function [te, pte, med] = gcmi_te(Xc, Uc, rows, P, qx, tgt, src, nshift, minshift)
% Gaussian-copula conditional transfer entropy (bits) src -> tgt given the
% other variables and respiration: Gaussian CMI of the copula-normalized
% series (Ince et al., 2017), i.e. half the log ratio of the target's
% residual variances without and with the source's past, with the
% regressors of the target equation. Null: the source's lag vectors
% circularly shifted (>= minshift samples) against the rest. The restricted
% regression is factored once; each shift only re-residualizes the source.
n = size(Xc,1);
[Z, cols] = var_design(Xc, Uc, rows, P, qx);
m = cols(tgt,:);
isSrc = src_cols(n, max(P(:)), size(Z,1), src) & m;
y = Xc(tgt, rows)';
[Q, ~] = qr(Z(m & ~isSrc,:)', 0);
ry = y - Q*(Q'*y);
RSSr = ry'*ry;
W = Z(isSrc,:)';
te = te_fw(W, Q, ry, RSSr);
nObs = numel(rows);
null = zeros(nshift,1);
for k = 1:nshift
    s = minshift + randi(nObs - 2*minshift);
    null(k) = te_fw(circshift(W, s, 1), Q, ry, RSSr);
end
pte = (1 + sum(null >= te)) / (nshift + 1);
med = median(null);


function te = te_fw(W, Q, ry, RSSr)
% Frisch-Waugh: residual of ry on the source lags residualized on Q
Wr = W - Q*(Q'*W);
e = ry - Wr*(Wr \ ry);
te = 0.5*log2(RSSr/(e'*e));


function Y = detrend_std(X)
% Remove the mean and linear trend of each run of complete samples (runs
% are separated by gaps: removed artifacts, separate blocks), then scale
% each row to unit variance. Per run: a level shift between blocks (e.g.
% a change of state that moves EEG power and heart period together) would
% otherwise look like lagged prediction across the whole recording.
Y = nan(size(X));
ok = all(isfinite(X), 1);
d = diff([false ok false]);
runs = [find(d == 1); find(d == -1) - 1]';
for k = 1:size(runs,1)
    idx = runs(k,1):runs(k,2);
    if numel(idx) < 10, continue; end                  % too short to fit lags
    t = (1:numel(idx))';
    D = [ones(numel(idx),1) t - mean(t)];
    Y(:, idx) = X(:, idx) - (D*(D \ X(:, idx)'))';
end
s = std(Y, 0, 2, 'omitnan');
Y = Y ./ s;


function Y = copnorm_rows(X)
% Copula normalization: ranks -> standard normal (NaN kept)
Y = nan(size(X));
for i = 1:size(X,1)
    ok = find(isfinite(X(i,:)));
    [~, o] = sort(X(i,ok));
    rk = zeros(1,numel(ok)); rk(o) = 1:numel(ok);
    Y(i,ok) = -sqrt(2)*erfcinv(2*rk/(numel(ok)+1));
end


function [lb, pv] = ljungbox(E, rows, h, p)
% Smallest Ljung-Box p-value across the residual series (lags 1..h, only
% between time points that were both fitted), and the p-value of each
n = size(E,1);
pos = nan(1, max(rows)); pos(rows) = 1:numel(rows);
pv = zeros(n,1);
for i = 1:n
    e = E(i,:) - mean(E(i,:));
    Q = 0; N = numel(e);
    for k = 1:h
        a = rows(rows - k > 0);
        a = a(isfinite(pos(a - k)));
        r = sum(e(pos(a)).*e(pos(a-k))) / sum(e.^2);
        Q = Q + r^2/(N - k);
    end
    Q = N*(N+2)*Q;
    df = max(h - p(min(i, numel(p))), 1);
    pv(i) = gammainc(Q/2, df/2, 'upper');            % chi2 upper tail
end
lb = min(pv);


function pc = partial_corr(V)
P = V \ eye(size(V));
d = sqrt(diag(P));
pc = -P ./ (d*d');
pc(logical(eye(size(V)))) = 1;


function q = fdr_bh(p)
% Benjamini-Hochberg adjusted p-values (q-values)
p = p(:); m = numel(p);
[ps, o] = sort(p);
qs = ps .* m ./ (1:m)';
qs = min(1, flipud(cummin(flipud(qs))));
q = nan(m,1); q(o) = qs;


function q = prctile_cols(X, pc)
% Percentile per column (no Statistics Toolbox)
Xs = sort(X, 1);
n = size(X,1);
k = min(max(round(pc/100*n), 1), n);
q = Xs(k, :);


function s = set_default(s, f, v)
if ~isfield(s, f) || isempty(s.(f)), s.(f) = v; end
