% GET_COUPLING_SERIES - Slow brain, heart and respiration series for the
% brain-heart coupling analysis, on a common uniform time grid (2 Hz by
% default).
%
% Heart rate variability lives below 0.4 Hz (LF 0.04-0.15 Hz, HF 0.15-0.4
% Hz), so the coupling is measured between slow fluctuations: the power
% envelopes of the EEG bands (per region), the NN interval series and
% respiration, all sampled at a few Hz (Zanetti et al., 2019; Schulz et
% al., 2019; Pernice et al., 2022).
%
% Usage:
%   S = get_coupling_series(X, fs, chanlocs, beats, NN_t, NN, resp, opt)
%
% INPUTS:
%   X        - nChan x nPts EEG (uV) on its original time line; NaN where
%              samples were removed (e.g. by ASR). Each run of valid
%              samples is processed separately.
%   fs       - sampling rate (Hz)
%   chanlocs - EEGLAB channel locations (for the 'auto' regions)
%   beats    - sample indices of the heartbeats (real beats only), for the
%              removal of the heartbeat-evoked EEG ('cfa')
%   NN_t, NN - NN intervals (s) and their times (s, beat closing each one)
%   resp     - struct: .source ('channel', 'ecg' or 'none'); .signal
%              (1 x nPts at fs) for 'channel'; .t and .amp (beat times, s,
%              and QRS amplitudes) for 'ecg' (ECG-derived respiration)
%   opt      - struct (all optional):
%     .grid       - 'beat' (default): one sample per heartbeat, as in
%                   cardiovascular variability analysis (Baselli et al.,
%                   1997; Faes et al., 2014): the NN interval itself, the
%                   mean band power and respiration over that interval. No
%                   interpolation: a spline-interpolated NN series on a
%                   uniform grid faster than the heart rate has structure
%                   between the beats that autoregressive models cannot
%                   whiten (in simulation: Ljung-Box failures in 88% of the
%                   heart equations at 2 Hz, and inflated tests). Frequencies
%                   are then in Hz at the mean heart period.
%                   'uniform': a grid at .fs_out
%     .fs_out     - rate of the uniform grid (Hz; default 2)
%     .bands      - nBands x 2 band limits (Hz; default [1 4; 4 8; 8 13; 13 30])
%     .band_names - names (default {'delta' 'theta' 'alpha' 'beta'})
%     .regions    - 'auto' (frontal, central, temporal, posterior; default),
%                   'all', 'channels', or {name, {labels}; ...}
%     .cfa        - cardiac field artifact, before the band powers:
%                   'evoked' (default): the heartbeat-locked average EEG
%                   (cardiac field and heartbeat-evoked potential, by
%                   least-squares deconvolution so that overlapping windows
%                   are handled) is subtracted from each channel; 'ssp': the
%                   first .ssp_n spatial components of that average are
%                   projected out at all times (removes the field whatever
%                   its beat-to-beat amplitude, and the brain activity with
%                   the same topographies); 'both'; 'none'
%     .cfa_window - [start end] (s) around each beat (default [-0.3 0.7])
%     .ssp_n      - SSP components (default 2)
%     .min_seg    - shortest run of valid samples kept (s; default 30)
%     .bridge     - gaps up to this length (s; default 1) are bridged for
%                   the filters only (linear interpolation); the beats or
%                   grid windows touching them are left out
%
% OUTPUT (S fields):
%   t          - grid times (s), 1 x T ('beat': the beat closing each
%                interval)
%   grid       - 'beat' or 'uniform'; fs_out - samples per second (beat
%                grid: 1 / mean NN)
%   eeg        - nRegions x nBands x T, log10 band power (uV^2), NaN where
%                not available
%   heart      - 1 x T NN interval (ms), cubic interpolation of the beats
%   resp       - 1 x T respiration (arbitrary units), empty if none
%   regions    - struct array: name, chans (indices), labels
%   bands, band_names, fs_out, resp_source
%   cfa        - struct: evoked template variance as a share of the EEG
%                variance per channel (share_var), window
%
% Why remove the heartbeat-evoked EEG: each heartbeat adds the cardiac
% field (QRS, T wave) to the EEG, so the power of any window grows with the
% number of beats it contains. Left in, it would make EEG power follow the
% heart rate (a spurious heart -> brain coupling at lag 0-1 beat). Its
% average waveform is removed here as for the induced HRSP of the HEP
% analysis; its beat-to-beat amplitude changes (respiration) remain, which
% is one reason to keep respiration in the model.
%
% Copyright (C) - Cedric Cannard, BrainBeats 2026

function S = get_coupling_series(X, fs, chanlocs, beats, NN_t, NN, resp, opt)

if nargin < 8, opt = struct; end
opt = set_default(opt, 'fs_out', 2);
opt = set_default(opt, 'bands', [1 4; 4 8; 8 13; 13 30]);
opt = set_default(opt, 'band_names', {'delta' 'theta' 'alpha' 'beta'});
opt = set_default(opt, 'regions', 'auto');
opt = set_default(opt, 'cfa', 'evoked');
opt = set_default(opt, 'cfa_window', [-0.3 0.7]);
opt = set_default(opt, 'ssp_n', 2);
opt = set_default(opt, 'vascular', {});
opt = set_default(opt, 'resp_phase', false);
opt = set_default(opt, 'bridge', 1);
opt = set_default(opt, 'min_seg', 30);
opt = set_default(opt, 'decim', 'bin');
opt = set_default(opt, 'win', 1/opt.fs_out);
opt = set_default(opt, 'nn_interp', 'spline');
opt = set_default(opt, 'grid', 'beat');
[nChan, nPts] = size(X);
nb = size(opt.bands,1);
isBeat = strcmpi(opt.grid, 'beat');
if isBeat
    fso = 1; step = 1;                      % runs need no grid alignment
else
    fso = opt.fs_out;
    if abs(fs/fso - round(fs/fso)) > 1e-9
        error('get_coupling_series: the EEG sampling rate (%g Hz) must be a multiple of the grid rate (%g Hz).', fs, fso)
    end
    step = round(fs/fso);                   % EEG samples per grid point
end

% Regions
regions = make_regions(opt.regions, chanlocs, nChan);
nr = numel(regions);

% Short gaps (e.g. ASR bursts) bridged for the filters only: linear
% interpolation across gaps of <= .bridge s; every interval or grid
% window touching one is set to NaN below. Otherwise many short removals
% would leave no run long enough for the filters and the model.
ok = all(isfinite(X), 1);
bridged = false(1, nPts);
d = diff([true ok true]);
g0 = find(d == -1); g1 = find(d == 1) - 1;          % gap start / end samples
for k = 1:numel(g0)
    if g0(k) > 1 && g1(k) < nPts && g1(k) - g0(k) + 1 <= opt.bridge*fs
        a = g0(k) - 1; b = g1(k) + 1;
        X(:, g0(k):g1(k)) = X(:, a) + (X(:, b) - X(:, a)) .* ((g0(k):g1(k)) - a)/(b - a);
        bridged(g0(k):g1(k)) = true;
    end
end
cb = [0 cumsum(bridged)];

% Runs of valid samples, starting on the grid
ok = all(isfinite(X), 1);
d = diff([false ok false]);
seg = [find(d == 1); find(d == -1) - 1]';
for k = 1:size(seg,1)
    seg(k,1) = seg(k,1) + mod(step - mod(seg(k,1)-1, step), step);   % first sample on the grid
end
seg(seg(:,2) - seg(:,1) + 1 < opt.min_seg*fs, :) = [];
if isempty(seg)
    error('get_coupling_series: no run of clean EEG longer than %g s.', opt.min_seg)
end

% Filters: broadband (lowest to highest band edge), then each band
fbb = bandpass_fir(fs, [min(opt.bands(:,1)) max(opt.bands(:,2))]);
fband = cell(1,nb);
for b = 1:nb, fband{b} = bandpass_fir(fs, opt.bands(b,:)); end
% Grid samples lost at each run edge: filter transients (+ the anti-alias
% filter of resample, 10 lobes each side)
edgeS = ceil(max(cellfun(@numel, [fband {fbb}]))/2);    % EEG samples
edge = ceil(edgeS/fs * fso) + strcmpi(opt.decim,'resample')*10;
Lwin = round(opt.win*fs);

% 1. Broadband filter per run (time x channels)
Y = nan(nPts, nChan);
for k = 1:size(seg,1)
    idx = seg(k,1):seg(k,2);
    Y(idx,:) = filtfilt(fbb, 1, double(X(:,idx))');
end

% 2. Cardiac field artifact. 'ssp': its first spatial components (from the
% heartbeat-locked average) projected out at all times, whatever the
% beat-to-beat amplitude; 'evoked': the heartbeat-locked average waveform
% subtracted (least-squares deconvolution, so overlapping windows are
% handled); 'both': SSP, then the evoked waveform of what remains
S.cfa = struct('method', opt.cfa, 'window', opt.cfa_window, 'share_var', [], 'ssp_patterns', [], 'ssp_explained', []);
if ~strcmpi(opt.cfa, 'none') && ~isempty(beats)
    lag = round(opt.cfa_window(1)*fs):round(opt.cfa_window(2)*fs);
    L = numel(lag);
    beats = beats(:);
    r = beats + lag;                                  % nBeats x L sample indices
    c = repmat(1:L, numel(beats), 1);
    keep = r >= 1 & r <= nPts;
    D = sparse(r(keep), c(keep), 1, nPts, L);
    valid = all(isfinite(Y), 2) & ~bridged(:);
    Dv = D(valid,:);
    G = full(Dv'*Dv);
    G = G + 1e-10*trace(G)/L*eye(L);                  % only guards a degenerate (strictly periodic) beat train
    v0 = var(Y(valid,:));
    if any(strcmpi(opt.cfa, {'ssp' 'both'}))
        T = G \ full(Dv'*Y(valid,:));                 % L x nChan heartbeat-locked average
        [~, sv, W] = svd(T - mean(T,1), 'econ');      % spatial components
        ns = min(opt.ssp_n, nChan-1);
        W = W(:, 1:ns);
        S.cfa.ssp_patterns = W;
        S.cfa.ssp_explained = diag(sv(1:ns,1:ns)).^2 / sum(diag(sv).^2);
        Y = Y - (Y*W)*W';                             % projection, all times
    end
    if any(strcmpi(opt.cfa, {'evoked' 'both'}))
        T = G \ full(Dv'*Y(valid,:));                 % L x nChan template (Cholesky through backslash)
        Y = Y - D*T;
    end
    S.cfa.share_var = 1 - var(Y(valid,:)) ./ v0;      % share of each channel's variance removed
end

% 3. Band power envelopes, averaged per region, to the grid ('beat': mean
% power over each RR interval, the interval that ends with the beat)
NN_t = NN_t(:)'; NN = NN(:)';
if isBeat
    T = numel(NN);
    S.t = NN_t;
    bEnd = round(NN_t*fs);                        % last EEG sample of each interval (sample 1 = t 0)
    bStart = bEnd - round(NN*fs) + 1;
else
    T = floor((nPts-1)/step) + 1;
    S.t = (0:T-1)/fso;
end
S.eeg = nan(nr, nb, T);
for k = 1:size(seg,1)
    idx = seg(k,1):seg(k,2);
    g = (seg(k,1)-1)/step + 1;                        % grid index of the run's first sample
    if isBeat                                     % intervals inside the run, away from filter edges and bridged gaps
        inRun = find(bStart >= seg(k,1) + edgeS & bEnd <= seg(k,2) - edgeS);
        inRun = inRun(cb(bEnd(inRun)+1) - cb(bStart(inRun)) == 0);
    end
    for b = 1:nb
        P = abs(hilbert(filtfilt(fband{b}, 1, Y(idx,:)))).^2;    % time x channels
        for ir = 1:nr
            pr = mean(P(:, regions(ir).chans), 2);
            if isBeat
                cs = [0; cumsum(pr)];
                lo = bStart(inRun)' - seg(k,1) + 1; hi = bEnd(inRun)' - seg(k,1) + 1;
                S.eeg(ir, b, inRun) = reshape(log10((cs(hi+1) - cs(lo)) ./ (hi - lo + 1)), 1, 1, []);
                continue
            end
            if strcmpi(opt.decim, 'bin')
                z = log10(to_grid(pr, step, Lwin, 'bin'));    % mean power in each window
            else
                z = to_grid(log10(pr), step, Lwin, 'resample');
            end
            z([1:edge end-edge+1:end]) = NaN;
            c = seg(k,1) + (0:numel(z)-1)*step;           % grid samples: NaN where the window touches a bridged gap
            lo = max(c - floor(Lwin/2), 1); hi = min(c + ceil(Lwin/2) - 1, nPts);
            z(cb(hi+1) - cb(lo) > 0) = NaN;
            S.eeg(ir, b, g:g+numel(z)-1) = z;
        end
    end
end

% 4. Heart: NN interval (ms); 'beat': the intervals themselves (implausible
% ones, e.g. unfilled gaps, set to NaN); uniform grid: interpolated, with
% nothing across gaps > 3 s
if isBeat
    S.heart = NN*1000;
    S.heart(NN > 2 | NN < 0.25) = NaN;
else
    S.heart = interp1(NN_t, NN*1000, S.t, opt.nn_interp, NaN);
    gap = find(diff(NN_t) > 3);
    for k = gap(:)'
        S.heart(S.t > NN_t(k) & S.t < NN_t(k+1)) = NaN;
    end
end

% 5. PPG pulse of each interval: the first pulse onset .ppg.window (s) after
% the beat that opens the interval (ECG beats: [0.04 0.5]; PPG beats: the
% onset itself). Per interval: pulse arrival time from that beat (ms; with
% ECG beats), pulse amplitude (foot to peak) and foot value (baseline
% intensity), from the PPG low-passed at 8 Hz (no high-pass, which would
% remove the respiratory baseline variation)
ppgF = [];
if isfield(opt,'ppg') && isstruct(opt.ppg) && isfield(opt.ppg,'onsets') && ~isempty(opt.ppg.onsets)
    ppgF = ppg_features(opt.ppg, fs, NN_t - NN);
    if ~isBeat
        for f = {'pat' 'amp' 'foot'}
            ppgF.(f{1}) = interp1(NN_t, ppgF.(f{1}), S.t, 'pchip', NaN);
        end
    end
end

% 6. Respiration: one or several proxies (struct array), all exogenous
% inputs of the model. 'channel': a respiration belt or airflow; 'ecg':
% the R-S amplitude of each QRS; 'ppg_riiv' / 'ppg_riav': the PPG pulse
% foot (intensity) / amplitude variation
S.resp = []; names = {};
if isstruct(resp) && isfield(resp,'source')
    for ir = 1:numel(resp)
        z = [];
        switch lower(resp(ir).source)
            case 'channel'
                z = nan(1,T); zc = nan(1,T); zs = nan(1,T);
                flp = lowpass_fir(fs, 1);
                for k = 1:size(seg,1)
                    idx = seg(k,1):seg(k,2);
                    rs = filtfilt(flp, 1, double(resp(ir).signal(idx))');
                    if isBeat                     % mean over each RR interval
                        inRun = find(bStart >= seg(k,1) + edgeS & bEnd <= seg(k,2) - edgeS);
                        cs = [0; cumsum(rs)];
                        lo = bStart(inRun)' - seg(k,1) + 1; hi = bEnd(inRun)' - seg(k,1) + 1;
                        z(inRun) = (cs(hi+1) - cs(lo)) ./ (hi - lo + 1);
                        if opt.resp_phase         % respiratory phase at the beat closing the interval
                            ph = angle(hilbert(rs - movmean(rs, round(8*fs))));
                            zc(inRun) = cos(ph(hi)); zs(inRun) = sin(ph(hi));
                        end
                        continue
                    end
                    g = (seg(k,1)-1)/step + 1;
                    zz = to_grid(rs, step, Lwin, opt.decim);
                    zz([1:edge end-edge+1:end]) = NaN;
                    z(g:g+numel(zz)-1) = zz;
                end
            case 'ecg'
                z = interp1(resp(ir).t, resp(ir).amp, S.t, 'pchip', NaN);   % 'beat': the QRS amplitudes themselves
            case 'ppg_riiv'
                if ~isempty(ppgF), z = ppgF.foot; else, warning('No PPG: respiration ''ppg_riiv'' ignored.'); end
            case 'ppg_riav'
                if ~isempty(ppgF), z = ppgF.amp; else, warning('No PPG: respiration ''ppg_riav'' ignored.'); end
        end
        if ~isempty(z)
            S.resp(end+1,:) = z;
            names{end+1} = lower(resp(ir).source); %#ok<AGROW>
        end
        if strcmpi(resp(ir).source, 'channel') && opt.resp_phase && isBeat
            S.resp(end+1:end+2,:) = [zc; zs];
            names(end+1:end+2) = {'channel_cos' 'channel_sin'};
        end
    end
end
if isempty(names), S.resp_source = 'none'; else, S.resp_source = strjoin(names, '+'); end

% 7. Vascular series from the PPG, modelled like the heart period (not as
% inputs): pulse arrival time ('pat', ECG + PPG) and pulse amplitude ('pamp')
S.vasc = []; S.vasc_names = {};
for v = cellstr(opt.vascular)
    if isempty(v{1}), continue; end
    if isempty(ppgF), warning('No PPG: vascular series ignored.'); break; end
    switch lower(v{1})
        case 'pat',  S.vasc(end+1,:) = ppgF.pat;
        case 'pamp', S.vasc(end+1,:) = ppgF.amp;
        otherwise, error('Unknown vascular series ''%s'' (''pat'' or ''pamp'').', v{1})
    end
    S.vasc_names{end+1} = lower(v{1});
end

S.regions = regions;
S.bands = opt.bands;
S.band_names = opt.band_names;
S.grid = lower(opt.grid);
if isBeat
    S.fs_out = 1/mean(NN(NN >= 0.25 & NN <= 2));  % beats per second: frequency axis in Hz at the mean heart period
else
    S.fs_out = fso;
end
S.cfa_method = opt.cfa;


function F = ppg_features(ppg, fs, tStart)
% Pulse of each interval (first onset in tStart + ppg.window), its arrival
% time (ms), amplitude and foot value, from the PPG low-passed at 8 Hz
% (sign set to match ppg.sig, the polarity-corrected signal of GET_RR)
if ~isfield(ppg,'window') || isempty(ppg.window), ppg.window = [0.04 0.5]; end
raw = double(ppg.raw(:)');
ok = isfinite(raw);
x = raw; x(~ok) = 0;
x = filtfilt(lowpass_fir(fs, 8), 1, x); x(~ok) = NaN;
if isfield(ppg,'sig') && ~isempty(ppg.sig)
    r = corr(x(ok)' - mean(x(ok)), double(ppg.sig(ok))' - mean(ppg.sig(ok)));
    if r < 0, x = -x; end
end
to = (double(ppg.onsets(:)') - 1)/fs;               % onset times (s)
K = numel(tStart);
[F.pat, F.amp, F.foot] = deal(nan(1,K));
w = round(0.45*fs);
for k = 1:K
    j = find(to >= tStart(k) + ppg.window(1), 1);
    if isempty(j) || to(j) > tStart(k) + ppg.window(2), continue; end
    o = ppg.onsets(j);
    if o + w > numel(x) || any(isnan(x(o:o+w))), continue; end
    F.pat(k) = (to(j) - tStart(k))*1000;
    F.foot(k) = x(o);
    F.amp(k) = max(x(o:o+w)) - x(o);
end


function z = to_grid(x, step, L, method)
% Values of the column x (starting on the grid) at the grid points: mean
% over a window of L samples centered on each point ('bin'), or
% anti-aliased decimation ('resample')
n = numel(x);
nz = floor((n-1)/step) + 1;
if strcmpi(method, 'bin')
    c = 1 + (0:nz-1)'*step;                   % grid samples
    lo = c - floor(L/2); hi = lo + L - 1;
    ok = lo >= 1 & hi <= n;
    cs = [0; cumsum(x)];
    z = nan(nz, 1);
    z(ok) = (cs(hi(ok)+1) - cs(lo(ok))) / L;
else
    mu = mean(x);
    z = resample(x - mu, 1, step) + mu;
    z = z(1:nz);
end


function b = bandpass_fir(fs, band)
% Zero-phase use (filtfilt) of a Hamming-windowed FIR; transition 25% of
% the low edge (0.5-2 Hz)
tw = max(0.5, min(0.25*band(1), 2));
ord = 2*ceil(3.3/tw*fs/2);
b = fir1(ord, band/(fs/2));


function b = lowpass_fir(fs, fc)
ord = 2*ceil(3.3/(0.5*fc)*fs/2);
b = fir1(ord, fc/(fs/2));


function regions = make_regions(spec, chanlocs, nChan)
labels = {chanlocs.labels};
if iscell(spec)
    regions = struct('name', {}, 'chans', {}, 'labels', {});
    for k = 1:size(spec,1)
        idx = find(ismember(lower(labels), lower(spec{k,2})));
        if isempty(idx), warning('Region %s: none of its channels found.', spec{k,1}); continue; end
        regions(end+1) = struct('name', spec{k,1}, 'chans', idx, 'labels', {labels(idx)}); %#ok<AGROW>
    end
elseif strcmpi(spec, 'channels')
    regions = struct('name', labels, 'chans', num2cell(1:nChan), 'labels', num2cell(labels));
elseif strcmpi(spec, 'all')
    regions = struct('name', 'all', 'chans', 1:nChan, 'labels', {labels});
else  % 'auto'
    th = [chanlocs.theta]; rd = [chanlocs.radius];
    if numel(th) ~= nChan || numel(rd) ~= nChan
        warning('No channel locations: one region with all channels.')
        regions = struct('name', 'all', 'chans', 1:nChan, 'labels', {labels});
        return
    end
    x = rd.*cosd(th); y = rd.*sind(th);             % x toward the nose, y to the right
    grp = repmat({'central'}, 1, nChan);
    grp(abs(y) >= 0.35 & abs(x) < 0.15) = {'temporal'};
    grp(x >= 0.15) = {'frontal'};
    grp(x <= -0.15) = {'posterior'};
    names = {'frontal' 'central' 'temporal' 'posterior'};
    regions = struct('name', {}, 'chans', {}, 'labels', {});
    for k = 1:numel(names)
        idx = find(strcmp(grp, names{k}));
        if isempty(idx), continue; end
        regions(end+1) = struct('name', names{k}, 'chans', idx, 'labels', {labels(idx)}); %#ok<AGROW>
    end
end


function s = set_default(s, f, v)
if ~isfield(s, f) || isempty(s.(f)), s.(f) = v; end
