function [RCMFE, scales] = compute_RCMFE(data, varargin)
% compute_RCMFE  Refined Composite Multiscale Fuzzy Entropy (RCMFE).
%
%   [RCMFE, scales] = compute_RCMFE(data, 'm', 2, 'r', 0.15, 'tau', 1, ...
%                                   'n', 2, 'coarsing', 'std', ...
%                                   'Mode', 'local', 'num_scales', 15)
%
% Inputs
%   data : EEGLAB EEG struct with .data OR numeric [n_ch x n_samp]
%
% Name-value parameters
%   'm'               : embedding dimension (default = 2)
%   'r'               : similarity bound (default = 0.15)
%   'tau'             : embedding delay for FuzEn (default = 1)
%   'n'               : fuzzy exponent (default = 2)
%   'coarsing'        : 'mean' | 'median' | 'trimmed mean' | 'std' | 'var'
%                       (default = 'std')
%   'Mode'            : 'local' | 'global' (default = 'local')
%   'num_scales'      : requested number of scales (default = 15)
%   'MinSamplesPerBin': minimum number of coarse bins required (default = 4)
%   'Parallel'        : true/false (default = true)
%   'Progress'        : true/false (default = true)
%
% Outputs
%   RCMFE  : [n_channels x numel(scales)] refined composite multiscale fuzzy entropy
%   scales : scales retained by the selected coarse-graining
%
% Notes
%   • RCMFE averages phi_m and phi_m+1 across shifted coarse-grained series,
%     then computes the entropy as log(phi_m_bar / phi_m1_bar).
%   • This differs from CMFE, which averages Fuzzy Entropy values directly.
%   • For 'std' and 'var', scale 1 is omitted because spread of a single-point
%     segment is not meaningful in the multiscale formulation.
%   • Use 'Mode','local' for closest agreement with the original local-centered
%     FuzEn-style multiscale literature.

p = inputParser;
p.addRequired('data', @(x) (isstruct(x) && isfield(x,'data')) || (isnumeric(x) && ndims(x)==2));
p.addParameter('m', 2,                 @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('r', 0.15,              @(x) isnumeric(x) && isscalar(x) && x>0 && x<2);
p.addParameter('tau', 1,               @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('n', 2,                 @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('coarsing','std',       @(s) any(strcmpi(s,{'median','mean','trimmed mean','trimmed','tmean','trim20','std','sd','standard deviation','var','variance'})));
p.addParameter('Mode','local',         @(s) any(strcmpi(s,{'local','global'})));
p.addParameter('num_scales', 15,       @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('MinSamplesPerBin', 4,  @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('Parallel', true,       @(x) islogical(x) && isscalar(x));
p.addParameter('Progress', true,       @(x) islogical(x) && isscalar(x));
p.parse(data, varargin{:});

m            = p.Results.m;
r            = p.Results.r;
tau          = p.Results.tau;
n_exp        = p.Results.n;
coarseType   = p.Results.coarsing;
modeType     = lower(p.Results.Mode);
nScales_req  = p.Results.num_scales;
minBinsHard  = p.Results.MinSamplesPerBin;
parallelMode = p.Results.Parallel;
showProg     = p.Results.Progress;

if isstruct(data), X = double(data.data); else, X = double(data); end
if size(X,1) > size(X,2), X = X.'; end
[nch, nSamp] = size(X);

S    = max(1, floor(nScales_req));
maxS = floor(nSamp / max(1, minBinsHard));
if S > maxS
    warning('compute_RCMFE:ReducingScales', ...
        'Reducing num_scales from %d to %d to keep >=%d coarse bins at every scale.', ...
        S, maxS, minBinsHard);
    S = maxS;
end
if S < 1, S = 1; end

isSpread = ismember(lower(strtrim(coarseType)), {'sd','std','standard deviation','var','variance'});
if isSpread
    scales = 2:S;
else
    scales = 1:S;
end

Xz = zscore_channels_local(X);
RCMFE = nan(nch, numel(scales));

if showProg
    state = ternary(parallelMode && ~isempty(ver('parallel')), 'on', 'off');
    fprintf('RCMFE: %d ch | m=%g, tau=%g, r=%g, n=%g | coarse=%s | mode=%s | scales=%d:%d | parallel=%s\n', ...
        nch, m, tau, r, n_exp, upperLabel(coarseType), upper(modeType), scales(1), scales(end), state);
end

useWB = ~parallelMode && usejava('desktop') && showProg;
hWB = [];
if useWB
    try, hWB = waitbar(0, 'Computing RCMFE...', 'Name', 'compute_RCMFE'); catch, hWB = []; end
end

useDQ = parallelMode && ~isempty(ver('parallel')) && showProg;
if useDQ
    dq = parallel.pool.DataQueue;
    nDone = 0;
    afterEach(dq, @notifyProgress);
end

if parallelMode && ~isempty(ver('parallel'))
    parfor ch = 1:nch
        sig = Xz(ch,:);
        v = nan(1, numel(scales));

        for ii = 1:numel(scales)
            s = scales(ii);

            phi_m_all  = nan(1, s);
            phi_m1_all = nan(1, s);

            for off = 1:s
                xoff = sig(off:end);
                Loff = floor(numel(xoff) / s) * s;
                nBins_off = Loff / s;
                if nBins_off < max(minBinsHard, m+1), continue; end

                Y  = reshape(xoff(1:Loff), s, []);
                cg = coarsegrain(Y, coarseType);

                [~, pm, pm1] = fuzz_engine_raw(cg, m, r, n_exp, tau, ...
                    'exponential', false, 2000, 2.0, modeType);
                phi_m_all(off)  = double(pm);
                phi_m1_all(off) = double(pm1);
            end

            good = isfinite(phi_m_all) & isfinite(phi_m1_all) & (phi_m_all > 0) & (phi_m1_all > 0);
            if any(good)
                phi_m_bar  = mean(phi_m_all(good));
                phi_m1_bar = mean(phi_m1_all(good));
                v(ii) = log(phi_m_bar / phi_m1_bar);
            end
        end

        RCMFE(ch,:) = v;
        if useDQ, send(dq,1); end
    end
else
    for ch = 1:nch
        sig = Xz(ch,:);
        v = nan(1, numel(scales));

        for ii = 1:numel(scales)
            s = scales(ii);

            phi_m_all  = nan(1, s);
            phi_m1_all = nan(1, s);

            for off = 1:s
                xoff = sig(off:end);
                Loff = floor(numel(xoff) / s) * s;
                nBins_off = Loff / s;
                if nBins_off < max(minBinsHard, m+1), continue; end

                Y  = reshape(xoff(1:Loff), s, []);
                cg = coarsegrain(Y, coarseType);

                [~, pm, pm1] = fuzz_engine_raw(cg, m, r, n_exp, tau, ...
                    'exponential', false, 2000, 2.0, modeType);
                phi_m_all(off)  = double(pm);
                phi_m1_all(off) = double(pm1);
            end

            good = isfinite(phi_m_all) & isfinite(phi_m1_all) & (phi_m_all > 0) & (phi_m1_all > 0);
            if any(good)
                phi_m_bar  = mean(phi_m_all(good));
                phi_m1_bar = mean(phi_m1_all(good));
                v(ii) = log(phi_m_bar / phi_m1_bar);
            end
        end

        RCMFE(ch,:) = v;
        if ~isempty(hWB) && isvalid(hWB)
            try, waitbar(ch/nch, hWB, sprintf('Computing RCMFE... (%d/%d)', ch, nch)); catch, end
        end
    end
    if ~isempty(hWB) && isvalid(hWB), try, close(hWB); catch, end, end
end

    function notifyProgress(~)
        nDone = nDone + 1;
        step = max(1, round(0.05*nch));
        if nDone == 1 || nDone == nch || mod(nDone, step) == 0
            fprintf('  progress: ch %d/%d\n', nDone, nch);
        end
    end
end

% -------------------------------------------------------------------------
function Xz = zscore_channels_local(X)
Xz = X;
for c = 1:size(X,1)
    x = X(c,:);
    mu = mean(x, 'omitnan');
    sd = std(x, 0, 'omitnan');
    if ~isfinite(sd) || sd == 0
        Xz(c,:) = 0;
    else
        Xz(c,:) = (x - mu) ./ sd;
    end
end
end

% -------------------------------------------------------------------------
function cg = coarsegrain(Y, ct)
switch lower(strtrim(ct))
    case 'mean'
        cg = mean(Y, 1, 'omitnan');
    case 'median'
        cg = median(Y, 1, 'omitnan');
    case {'trimmed mean','trimmed','tmean','trim20'}
        pct = 20;
        if exist('trimmean','file') == 2 && all(isfinite(Y(:)))
            cg = trimmean(Y, pct, 1);
        else
            nSeg = size(Y,2); cg = nan(1,nSeg);
            kfrac = pct/200;
            for j = 1:nSeg
                col = Y(:,j); col = col(isfinite(col));
                if isempty(col), cg(j)=NaN; continue; end
                k = floor(kfrac*numel(col));
                if 2*k >= numel(col)
                    cg(j) = NaN;
                else
                    col = sort(col);
                    cg(j) = mean(col(k+1:end-k));
                end
            end
        end
    case {'sd','std','standard deviation'}
        cg = std(Y, 0, 1, 'omitnan');
    case {'variance','var'}
        cg = var(Y, 0, 1, 'omitnan');
    otherwise
        error('compute_RCMFE:BadCoarsing', 'Unknown coarsing "%s".', ct);
end
cg = cg(:).';
end

% -------------------------------------------------------------------------
function s = upperLabel(coarseType)
cl = lower(strtrim(coarseType));
if any(strcmp(cl,{'sd','std','standard deviation'}))
    s = 'STD';
elseif any(strcmp(cl,{'var','variance'}))
    s = 'VAR';
elseif strcmp(cl,'mean')
    s = 'MEAN';
elseif strcmp(cl,'median')
    s = 'MEDIAN';
elseif any(strcmp(cl,{'trimmed mean','trimmed','tmean','trim20'}))
    s = 'TRIM20';
else
    s = upper(coarseType);
end
end

% -------------------------------------------------------------------------
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end