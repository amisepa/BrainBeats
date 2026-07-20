function [MFE, scales, CI] = compute_MFE(data, varargin)
% compute_MFE  Multiscale Fuzzy Entropy (Costa-style coarse-graining) across channels.
%
%   [MFE, scales, CI] = compute_MFE(data, 'm', 2, 'r', 0.15, 'tau', 1, ...
%                               'n', 2, 'coarsing','mean', 'num_scales', 30, ...
%                               'zNorm', 0, 'IncludeScale1', false, ...
%                               'Mode','local', 'Kernel','exponential', ...
%                               'MinSamplesPerBin', 4, 'StableMinBins', 100, ...
%                               'BlockSize', 2000, 'pdistMaxGB', 2.0, ...
%                               'Parallel', true, 'Progress', true)
%
%   ALGORITHM (classic, default):
%     1. Each channel is z-scored ONCE over the whole recording.
%     2. At each temporal scale s, the series is coarse-grained into
%        non-overlapping windows of length s (mean by default).
%     3. Fuzzy entropy of the coarse-grained series is computed with a
%        tolerance r that is FIXED across scales (r is a fraction of the
%        original SD). The coarse-grained series is NOT re-normalised, so the
%        natural decrease of variance across scales is preserved. This is done
%        by calling compute_FuzzEn with 'Normalize', false.
%
%   Outputs:
%     MFE    : [n_channels x n_scales] fuzzy entropy per channel and scale
%     scales : vector of scale factors actually used
%     CI     : [n_channels x 1] complexity index = sum of MFE across scales
%
%   GUI-tunable parameters (same names/behaviour as before):
%     'm'          embedding dimension               (default 2)
%     'r'          tolerance, fraction of SD          (default 0.15)
%     'tau'        time delay for embedding           (default 1)
%     'n'          fuzzy exponent                      (default 2)
%     'coarsing'   coarse-graining statistic          (default 'mean')
%                  {'mean','median','std','var'}
%     'num_scales' requested number of scales         (default 30)
%
%   Command-line / advanced parameters:
%     'zNorm'         varying-tolerance rescaling, integer 0-4 (default 0 = OFF)
%                     0  fixed tolerance across scales (classic)
%                     1  r*std  of each coarse-grained series
%                     2  r*var  of each coarse-grained series
%                     3  r*mad  (mean absolute deviation)
%                     4  r*mad  (median absolute deviation)
%     'IncludeScale1' include scale 1 (= plain FuzzEn)  (default false)
%     'Mode'          FuzzEn detrending, 'local' | 'global' (default 'local')
%     'Kernel'        'exponential' | 'gaussian'          (default 'exponential')
%     'MinSamplesPerBin','StableMinBins','BlockSize','pdistMaxGB',
%     'Parallel','Progress' as before.
%
%   NOTE ON DEFAULTS: the default coarse-graining is now 'mean' (classic
%   multiscale entropy). 'std'/'var' coarse-graining yield generalized
%   (variance-based) MFE and remain available as options.

% ---------------- Parse inputs ----------------
p = inputParser;
p.addRequired('data', @(x) (isstruct(x) && isfield(x,'data')) || (isnumeric(x) && ndims(x)==2));
p.addParameter('m', 2,                 @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('r', 0.15,              @(x) isnumeric(x) && isscalar(x) && x>0 && x<2);
p.addParameter('tau', 1,               @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('n', 2,                 @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('coarsing','mean',      @(s) any(strcmpi(s,{'median','mean','std','sd','standard deviation','var','variance'})));
p.addParameter('num_scales', 30,       @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('zNorm', 0,             @(x) isnumeric(x) && isscalar(x) && ismember(x,0:4));
p.addParameter('IncludeScale1', false, @(x) islogical(x) && isscalar(x));
p.addParameter('Mode', 'local',        @(s) any(strcmpi(s,{'local','global'})));
p.addParameter('Kernel', 'exponential',@(s) any(strcmpi(s,{'exponential','gaussian'})));
p.addParameter('MinSamplesPerBin', 4,  @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('StableMinBins', 100,   @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('BlockSize', 2000,      @(x) isnumeric(x) && isscalar(x) && x>=100);
p.addParameter('pdistMaxGB', 2.0,      @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('Parallel', true,       @(x) islogical(x) && isscalar(x));
p.addParameter('Progress', true,       @(x) islogical(x) && isscalar(x));
p.parse(data, varargin{:});

m             = p.Results.m;
r             = p.Results.r;
tau           = p.Results.tau;
n_exp         = p.Results.n;
coarseType    = p.Results.coarsing;
nScales_req   = p.Results.num_scales;
zNorm         = p.Results.zNorm;
incScale1     = p.Results.IncludeScale1;
modeType      = p.Results.Mode;
kernelType    = p.Results.Kernel;
minBinsAll    = p.Results.MinSamplesPerBin;
minBinsStable = p.Results.StableMinBins;
blockSize     = p.Results.BlockSize;
pdistMaxGB    = p.Results.pdistMaxGB;
parallelMode  = p.Results.Parallel;
showProgress  = p.Results.Progress;

% ---------------- Get numeric data ----------------
if isstruct(data)
    X = double(data.data);
else
    X = double(data);
end
if size(X,1) > size(X,2)
    X = X.';
end
[nch, nSamp] = size(X);

% ---------------- Cap S so every scale has >= minBinsAll coarse bins -----
S    = max(1, floor(nScales_req));
maxS = floor(nSamp / max(1, minBinsAll));
if S > maxS
    warning('compute_MFE:ReducingScales', ...
        'Reducing num_scales from %d to %d to keep >=%d coarse bins at every scale.', ...
        S, maxS, minBinsAll);
    S = maxS;
end
if S < 1
    S = 1;
end

if incScale1
    scales = 1:S;
else
    scales = 2:S;      % skip scale 1 (plain FuzzEn), preserved default
end
nScales = numel(scales);

% ---------------- Fill missing & z-score per channel (ONCE) --------------
for ch = 1:nch
    xi = X(ch,:);
    if any(~isfinite(xi))
        try
            xi = fillmissing(xi, 'linear', 'EndValues', 'nearest');
        catch
            isn = ~isfinite(xi);
            if any(isn)
                idx = find(~isn, 1, 'first');
                if ~isempty(idx), xi(1:idx-1) = xi(idx); end
                idx = find(~isn, 1, 'last');
                if ~isempty(idx), xi(idx+1:end) = xi(idx); end
                prev = [xi(1), xi(1:end-1)];
                next = [xi(2:end), xi(end)];
                bad  = isn & isfinite(prev) & isfinite(next);
                xi(bad) = 0.5 * (prev(bad) + next(bad));
            end
        end
        X(ch,:) = xi;
    end
end

Xz = X;
for c = 1:nch
    x  = X(c,:);
    mu = mean(x, 'omitnan');
    sd = std(x, 0, 'omitnan');
    if ~isfinite(sd) || sd == 0
        Xz(c,:) = 0;
    else
        Xz(c,:) = (x - mu) ./ sd;
    end
end
% After z-scoring SD = 1, so the fixed absolute tolerance is simply r.

% ---------------- Outputs & progress header ------------------------------
MFE = nan(nch, nScales);

coarseLabel = coarse_label(coarseType);

if showProgress
    parOn = parallelMode && ~isempty(ver('parallel'));
    fprintf('MFE: %d ch | m=%g, tau=%g, r=%g, n=%g | coarse=%s | zNorm=%d | scales=%d..%d | parallel=%s\n', ...
        nch, m, tau, r, n_exp, coarseLabel, zNorm, scales(1), scales(end), string(parOn));
end

useWB = showProgress && ~parallelMode && usejava('desktop');
hWB = [];
if useWB
    try
        hWB = waitbar(0, 'Computing Multiscale Fuzzy Entropy...', 'Name', 'compute_MFE');
    catch
        hWB = [];
    end
end

useDQ = parallelMode && ~isempty(ver('parallel')) && showProgress;
if useDQ
    dq = parallel.pool.DataQueue;
    nDone = 0;
    afterEach(dq, @notifyProgress);
end

% ---------------- Compute per channel ------------------------------------
if parallelMode && ~isempty(ver('parallel'))
    parfor ch = 1:nch
        MFE(ch,:) = mfe_one_channel(Xz(ch,:), m, r, n_exp, tau, coarseType, zNorm, ...
            modeType, kernelType, scales, minBinsAll, minBinsStable, ...
            blockSize, pdistMaxGB, false, ch, nch);
        if useDQ
            send(dq, 1);
        end
    end
else
    for ch = 1:nch
        MFE(ch,:) = mfe_one_channel(Xz(ch,:), m, r, n_exp, tau, coarseType, zNorm, ...
            modeType, kernelType, scales, minBinsAll, minBinsStable, ...
            blockSize, pdistMaxGB, showProgress, ch, nch);

        if ~isempty(hWB) && isvalid(hWB)
            try
                waitbar(ch/nch, hWB, sprintf('Computing MFE... (%d/%d)', ch, nch));
            catch
            end
        end
    end
    if ~isempty(hWB) && isvalid(hWB)
        try
            close(hWB);
        catch
        end
    end
end

% ---------------- Complexity index (area under the MFE curve) ------------
CI = sum(MFE, 2, 'omitnan');

    function notifyProgress(~)
        nDone = nDone + 1;
        step = max(1, round(0.05 * nch));
        if nDone == 1 || nDone == nch || mod(nDone, step) == 0
            fprintf('  progress: ch %d/%d\n', nDone, nch);
        end
    end
end

% ========================================================================
function v = mfe_one_channel(sig, m, r, n_exp, tau, coarseType, zNorm, ...
    modeType, kernelType, scales, minBinsAll, minBinsStable, ...
    blockSize, pdistMaxGB, showProgress, ch, nch)

v = nan(1, numel(scales));

for si = 1:numel(scales)
    s = scales(si);

    if s == 1
        cg = sig;
    else
        L  = floor(numel(sig)/s) * s;
        Y  = reshape(sig(1:L), s, []);   % s samples per non-overlapping window
        cg = coarsegrain(Y, coarseType); % 1 x nBins
    end

    nBins     = numel(cg);
    minNeeded = max([minBinsAll, m+1, minBinsStable]);
    if nBins < minNeeded
        if showProgress
            fprintf('  [drop] ch %d: scale %d nBins=%d (<%d)\n', ch, s, nBins, minNeeded);
        end
        continue
    end

    % Tolerance: FIXED across scales (zNorm=0), else rescaled per scale.
    if zNorm == 0
        r_s = r;                          % r is a fraction of the (unit) original SD
    else
        r_s = r * zNorm_scale(cg, zNorm);
    end

    % Fast FuzzEn engine. Normalize=false so the coarse-grained series is used
    % as-is and r_s is applied as an absolute tolerance (no per-scale z-score).
    v(si) = compute_FuzzEn(cg, 'm', m, 'n', n_exp, 'tau', tau, 'r', r_s, ...
        'Mode', modeType, 'Kernel', kernelType, 'Normalize', false, ...
        'BlockSize', blockSize, 'pdistMaxGB', pdistMaxGB, ...
        'Parallel', false, 'Progress', false);
end

if showProgress
    fprintf('  ch %3d/%3d: done\n', ch, nch);
end
end

% ========================================================================
function cg = coarsegrain(Y, coarseType)
% Coarse-grain down dim 1 (over the s samples of each window) -> 1 x nBins.
cl = lower(strtrim(coarseType));
switch cl
    case {'mean'}
        cg = mean(Y, 1);
    case {'median'}
        cg = median(Y, 1);
    case {'std','sd','standard deviation'}
        cg = std(Y, 0, 1);
    case {'var','variance'}
        cg = var(Y, 0, 1);
    otherwise
        cg = mean(Y, 1);
end
end

% ========================================================================
function f = zNorm_scale(x, zNorm)
% Per-scale tolerance rescaling factor.
switch zNorm
    case 1, f = std(x, 1);      % population standard deviation
    case 2, f = var(x, 1);      % population variance
    case 3, f = mad(x, 0);      % mean absolute deviation
    case 4, f = mad(x, 1);      % median absolute deviation
    otherwise, f = 1;
end
end

% ========================================================================
function lbl = coarse_label(coarseType)
cl = lower(strtrim(coarseType));
if any(strcmp(cl, {'sd','std','standard deviation'}))
    lbl = 'STD';
elseif any(strcmp(cl, {'var','variance'}))
    lbl = 'VAR';
elseif strcmp(cl, 'mean')
    lbl = 'MEAN';
elseif strcmp(cl, 'median')
    lbl = 'MEDIAN';
else
    lbl = upper(coarseType);
end
end