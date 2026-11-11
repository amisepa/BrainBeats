function [R, T] = validate_ppg_polarity(dsPath)
% VALIDATE_PPG_POLARITY - PPG polarity of GET_RR: the skewness of the slope
% (0.5-8 Hz band, as GET_RR computes it) vs quantile skewnesses of the same
% slope, on clean PPGs and with injected artifacts.
%
% Statistics (v: slope of the band-passed PPG; qN: N-th percentile of v):
%   skewness  - skewness(v), GET_RR until 2026-11-10
%   q1        - (q99 + q01 - 2 q50) / (q99 - q01)
%   q5        - (q95 + q05 - 2 q50) / (q95 - q05)
% Positive: the fast edge rises (blood volume, systolic peak up). Each is
% odd in the signal (s(-x) = -s(x)), so with artifacts of either sign the
% accuracy does not depend on the input polarity.
%
% Ground truth, independent of the edge asymmetry: from the ECG, the
% R-locked average of the band-passed PPG over one median RR interval; the
% pulse foot (minimum of a blood-volume PPG) comes before the systolic
% peak (maximum), so the input is inverted when its maximum comes first.
% The truth is kept when both halves of the recording agree. The Muse 2
% optics (tests/data/test_optics.csv, 8 channels, no ECG) are photodiode
% light intensity, inverted by construction.
%
% Recordings: the ds003838 PPGs (1000 Hz, with ECG), the sample dataset
% (sample_data/dataset.set, 250 Hz, PPG and ECG), the 8 Muse channels
% (64 Hz). Conditions, on the PPG at 250 Hz (ds003838 resampled after
% removing the mean; Muse at 64 Hz), 5 random repetitions each:
%   clean, and at the native rate for ds003838 ('native'), resampled with
%   zero-padded edges ('zero-pad edges', the trap of 2026-11-09)
%   steps     - 1 or 3 permanent baseline steps of 10 or 30 x the robust SD
%               of the band-passed PPG, random sign and position
%   motion    - 1, 5 or 15 bursts of 2 s of 0.5-4 Hz noise, 5 or 15 x that SD
%   dropouts  - 5 or 20 runs of 1-5 samples at zero (raw values)
%   clipping  - raw values clipped at their 90th (top) or 10th (bottom)
%               percentile
%
% Usage:
%   [R, T] = validate_ppg_polarity(dsPath)    % folder of sub-*/ecg/*.set
%
% Outputs:
%   R - struct array, one element per recording: name, fs, truth (+1 blood
%       volume, -1 inverted, NaN unclear), foot and peak latency of the
%       R-locked average (ms after R), stat (conditions x repetitions x
%       statistics)
%   T - table: condition, then the share of correct decisions per statistic
%       (recordings with a truth x repetitions)
%
% Results (2026-11-10; 73 recordings, 70 with a truth, all inverted; the
% halves disagree on sub-052, 056, 080, PPGs with scattered pulses; 89 s):
%   condition          n    skewness    q1      q5     (correct)
%   clean             70      98.6%   100%    100%    (skewness wrong on sub-084)
%   native (1000 Hz)  61       100%   100%    100%
%   zero-pad edges    61       8.2%   100%    100%
%   step 1 x10       350      97.4%   100%    100%
%   step 3 x10       350      93.1%   100%    100%
%   step 3 x30       350      52.9%   100%    100%
%   motion 1 x5      350      97.1%   100%    100%
%   motion 5 x5      350      92.3%   100%    100%
%   motion 5 x15     350      52.0%  88.3%    100%
%   motion 15 x15    350      54.0%  84.9%    100%
%   dropouts 5       350      67.4%   100%    100%
%   dropouts 20      350      61.4%  97.1%    100%
%   clip top 10%      70      98.6%   100%    100%
%   clip bottom 10%   70      98.6%   100%    100%
% Smallest q5 margin (truth x q5): clean 0.46, artifacts 0.30 (20 dropouts),
% 0.10 (15 motion bursts: 10% of the samples, beyond the 5% tails). On the
% three unclear PPGs |q5| < 0.03 in every condition: GET_RR warns below 0.2
% (never on a clean PPG; on 4 of 350 runs with 15 motion bursts). Since
% 2026-11-10 GET_RR uses q5: right on all 61 ds003838 PPGs with a truth at
% 1000 and 250 Hz, warning on sub-052, 056, 080 only.
%
% Copyright (C) - Cedric Cannard, 2026

bb = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(bb, 'functions'));
if ~exist('pop_loadset', 'file'), eeglab; close; end
statNames = {'skewness' 'q1' 'q5'};
cond = {'clean' 'native' 'zero-pad edges' 'step 1 x10' 'step 3 x10' 'step 3 x30' ...
    'motion 1 x5' 'motion 5 x5' 'motion 5 x15' 'motion 15 x15' 'dropouts 5' 'dropouts 20' 'clip top 10%' 'clip bottom 10%'};
nRep = 5;

% Recordings: raw PPG, its rate, and the ECG (empty for Muse)
recs = struct('name', {}, 'ppg', {}, 'ecg', {}, 'fs', {});
d = dir(fullfile(dsPath, 'sub-*'));
for k = 1:numel(d)
    f = fullfile(d(k).folder, d(k).name, 'ecg', [d(k).name '_task-rest_ecg.set']);
    if ~exist(f, 'file'), continue; end
    [~, E] = evalc('pop_loadset(f)');
    recs(end+1) = load_rec(E, d(k).name); %#ok<AGROW>
end
[~, E] = evalc('pop_loadset(''filename'', ''dataset.set'', ''filepath'', fullfile(bb, ''sample_data''))');
recs(end+1) = load_rec(E, 'sample dataset');
M = readmatrix(fullfile(bb, 'tests', 'data', 'test_optics.csv'));
fsM = 1/median(diff(M(:,1)));
for c = 1:8
    recs(end+1) = struct('name', sprintf('Muse optics ch%d', c), 'ppg', M(:,c+1)', 'ecg', [], 'fs', fsM); %#ok<AGROW>
end

R = struct('name', {}, 'fs', {}, 'truth', {}, 'foot', {}, 'peak', {}, 'stat', {});
for k = 1:numel(recs)
    x = recs(k).ppg;  fs0 = recs(k).fs;
    if isempty(recs(k).ecg)
        tr = -1;  lat = [NaN NaN];                     % Muse: light intensity
    else
        [tr, lat] = truth_from_ecg(x, recs(k).ecg, fs0);
    end
    if fs0 > 300                                       % ds003838: to 250 Hz
        [p, q] = rat(250/fs0, 1e-6);
        x250 = resample(x - mean(x), p, q) + mean(x);
        fs = fs0*p/q;
    else
        x250 = x;  fs = fs0;
    end
    S = nan(numel(cond), nRep, 3);
    S(1,1,:) = stats(x250, fs);
    if fs0 > 300
        S(2,1,:) = stats(x, fs0);
        S(3,1,:) = stats(resample(x, p, q), fs);
    end
    sdp = robust_sd(x250, fs);
    for c = 4:numel(cond)
        for rep = 1:nRep
            rng(100000*c + 1000*rep + k);
            S(c,rep,:) = stats(inject(x250, fs, cond{c}, sdp), fs);
            if startsWith(cond{c}, 'clip'), break; end     % deterministic
        end
    end
    R(end+1) = struct('name', recs(k).name, 'fs', fs0, 'truth', tr, 'foot', lat(1), 'peak', lat(2), 'stat', S); %#ok<AGROW>
    fprintf('%-18s truth %+d (foot %4.0f, peak %4.0f ms) clean %s\n', recs(k).name, tr, lat(1), lat(2), mat2str(squeeze(S(1,1,:))', 3));
end

% Summary: share of correct decisions per condition and statistic
ok = ~isnan([R.truth]);
fprintf('\n%d recordings, %d with a truth (%s unclear)\n', numel(R), nnz(ok), strjoin({R(~ok).name}, ', '));
T = table('Size', [numel(cond) 5], 'VariableTypes', {'string' 'double' 'double' 'double' 'double'}, ...
    'VariableNames', [{'condition' 'n'} statNames]);
fprintf('%-16s %6s %9s %9s %9s\n', 'condition', 'n', statNames{:});
for c = 1:numel(cond)
    C = []; tr = [];
    for k = find(ok)
        s = squeeze(R(k).stat(c,:,:));
        s = s(all(isfinite(s), 2), :);
        C = [C; sign(s)]; tr = [tr; repmat(R(k).truth, size(s,1), 1)]; %#ok<AGROW>
    end
    acc = mean(C == tr, 1);
    T(c,:) = {cond{c}, numel(tr), acc(1), acc(2), acc(3)};
    fprintf('%-16s %6d %8.1f%% %8.1f%% %8.1f%%\n', cond{c}, numel(tr), 100*acc);
end
for c = 1:3
    fprintf('\nWrong on "%s":\n', cond{c});
    for k = find(ok)
        s = squeeze(R(k).stat(c,1,:))';
        if any(isfinite(s)) && any(sign(s) ~= R(k).truth)
            fprintf('  %-18s truth %+d: %s\n', R(k).name, R(k).truth, mat2str(s, 3));
        end
    end
end
end

function r = load_rec(E, name)
iP = find(strcmpi({E.chanlocs.labels}, 'PPG'), 1);
iE = find(strcmpi({E.chanlocs.labels}, 'ECG'), 1);
x = double(E.data(iP,:));  e = double(E.data(iE,:));
ok = isfinite(x) & isfinite(e);
r = struct('name', name, 'ppg', x(ok), 'ecg', e(ok), 'fs', E.srate);
end

function s = stats(x, fs)
% polarity statistics of the slope, as GET_RR filters it
[b, a] = butter(2, [0.5 min(8, 0.45*fs)]/(fs/2));
v = diff(filtfilt(b, a, x - mean(x)));
q = quantile(v, [0.01 0.05 0.5 0.95 0.99]);
s = [skewness(v), (q(5) + q(1) - 2*q(3))/(q(5) - q(1)), (q(4) + q(2) - 2*q(3))/(q(4) - q(2))];
end

function [pol, lat] = truth_from_ecg(ppg, ecg, fs)
% polarity from the order of foot and systolic peak after the R-peak; NaN
% when the two halves of the recording disagree
tms = (0:numel(ecg)-1)/fs*1000; %#ok<NASGU> used in evalc
[~, ~, ~, pk] = evalc('get_RR(ecg, tms, struct(''heart_signal'', ''ecg'', ''fs'', fs))');
[b, a] = butter(2, [0.5 min(8, 0.45*fs)]/(fs/2));
y = filtfilt(b, a, ppg(:) - mean(ppg));
L = round(median(diff(pk)));
pk = pk(:);  pk = pk(pk + L - 1 <= numel(y));
h = round(numel(pk)/2);
P = zeros(1, 3);
for part = 1:3
    sel = {1:numel(pk), 1:h, h+1:numel(pk)};
    avg = mean(y(pk(sel{part}) + (0:L-1)), 1);
    [~, iMin] = min(avg);  [~, iMax] = max(avg);
    P(part) = 1 - 2*(iMax < iMin);
    if part == 1, lat = ([iMin iMax] - 1)/fs*1000; end
end
pol = P(1);
if numel(pk) < 20 || any(P(2:3) ~= P(1)), pol = NaN; end
if pol < 0, lat = lat([2 1]); end      % foot: the minimum once flipped
end

function sd = robust_sd(x, fs)
[b, a] = butter(2, [0.5 min(8, 0.45*fs)]/(fs/2));
y = filtfilt(b, a, x - mean(x));
sd = 1.4826*median(abs(y - median(y)));
end

function x = inject(x, fs, name, sdp)
% artifacts in raw units at random positions
N = numel(x);
tok = sscanf(regexprep(name, '[^0-9 ]', ' '), '%d')';
switch strtok(name)
    case 'step'
        for k = 1:tok(1)
            i0 = randi([round(0.05*N) round(0.95*N)]);
            x(i0:end) = x(i0:end) + sign(randn)*tok(2)*sdp;
        end
    case 'motion'
        [b, a] = butter(2, [0.5 4]/(fs/2));
        L = round(2*fs);  w = hann(L)';
        for k = 1:tok(1)
            n = filtfilt(b, a, randn(1, L + 2*L));
            n = n(L+1:2*L);  n = n/std(n);
            i0 = randi(N - L);
            x(i0:i0+L-1) = x(i0:i0+L-1) + tok(2)*sdp*w.*n;
        end
    case 'dropouts'
        for k = 1:tok(1)
            len = randi(5);  i0 = randi(N - len);
            x(i0:i0+len-1) = 0;
        end
    case 'clip'
        if contains(name, 'top'), x = min(x, prctile(x, 90)); else, x = max(x, prctile(x, 10)); end
end
end
