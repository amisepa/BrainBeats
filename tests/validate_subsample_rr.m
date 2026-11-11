function R = validate_subsample_rr(part, dataPath)
% VALIDATE_SUBSAMPLE_RR - Sub-sample beat times of GET_RR: are the RR
% intervals of a low-rate ECG or PPG closer to the truth with each beat
% refined to a fraction of a sample (vertex of the parabola through the
% peak and its two neighbors, what GET_RR returns since 2026-11-09) than on
% the sample grid?
%
% Parts:
%   'sim' - synthetic ECG with known beat times (20 recordings of 300 s) at
%           128, 250, 500 and 1000 Hz
%   'ecg' - ds003838 ECG (1000 Hz) resampled to 250 and 128 Hz; reference:
%           the beats detected at the native rate, refined with a spline
%   'ppg' - the same for the ds003838 PPG (pulse onsets, 'valleys')
%   'mit' - MIT-BIH Arrhythmia Database, MLII lead (360 Hz), resampled to
%           180 and 128 Hz; the paced records (102, 104, 107, 217) left out
% The detected beats are matched one-to-one to the reference within 50 ms,
% and intervals are compared over consecutive matched beats. Three beat
% times: 'grid' (the peak sample), 'get_RR' (the times GET_RR returns),
% 'spline' (maximum, or minimum for PPG onsets, of a cubic spline through
% +/-5 samples, on a 1/100-sample grid). The mean is removed before
% resampling: RESAMPLE pads with zeros, and the edge transients of a raw
% PPG (offset ~21 SD) flipped GET_RR's polarity on 59 of 64 recordings.
%
% Usage:
%   R = validate_subsample_rr('sim');
%   R = validate_subsample_rr('ecg', ds003838Path);  % folder of sub-*/ecg
%   R = validate_subsample_rr('ppg', ds003838Path);
%   R = validate_subsample_rr('mit', mitdbPath);     % .hea, .dat, RECORDS
%
% Output:
%   R - struct array, one element per recording and rate: rec, fs, fs0,
%       nbeat, nmatch, rr_err (ms RMS; grid, get_RR, spline), rmssd (ms) and
%       pnn50 (%) of the three, then of the reference on the same intervals.
%       A summary per rate is printed (recordings with at least 90% of their
%       beats matched).
%
% Results (2026-11-09; median over recordings of the interval error, ms
% RMS, sample grid -> GET_RR; RMSSD bias: median, ms; pNN50: largest error
% of a recording, points):
%   sim  (20)  128 Hz  3.23 -> 0.29   RMSSD +0.39 -> 0   pNN50 3.4 -> 0.9
%              250 Hz  1.65 -> 0.18         +0.12 -> 0         1.7 -> 0.9
%              500 Hz  0.83 -> 0.12         +0.03 -> 0         1.4 -> 0.3
%             1000 Hz  0.42 -> 0.09         +0.01 -> 0         0.9 -> 0.3
%   ecg  (64)  250 Hz  1.65 -> 0.04         +0.13 -> 0         3.0 -> 0.4
%              128 Hz  3.19 -> 0.19         +0.35 -> 0         3.3 -> 0.8
%   mit  (44)  180 Hz  2.29 -> 0.11         +0.09 -> 0         0.9 -> 0.1
%              128 Hz  3.22 -> 0.24         +0.18 -> 0         1.1 -> 0.3
%   ppg  (63)  250 Hz  1.87 -> 0.90         +0.13 -> 0         2.3 -> 0.7
%              128 Hz  3.30 -> 0.64         +0.40 -> 0         2.8 -> 0.4
% On the grid an interval is off by (1/fs)/sqrt(6) RMS (1.63 ms at 250 Hz).
% The spline is closer on real ECGs (ecg 0.02 and 0.11 ms, mit 0.04 and
% 0.12) with the same RMSSD and pNN50, and equal on the synthetic ECG and
% the PPG. The worst recordings (3-9 ms) are the same for all three beat
% times: detection, not timing. PPG: sub-080, with scattered pulses, is
% left out (its polarity is unclear, and GET_RR warns). With the skewness
% of the slope for the polarity (until 2026-11-10), 3 recordings changed
% polarity between rates and 61 were compared, with the same medians
% (tests/validate_ppg_polarity.m).
%
% Runtime: sim ~0.5 min, mit ~2 min, ecg and ppg ~4 min each.
%
% Copyright (C) - Cedric Cannard, 2026

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'functions'));
R = struct('rec', {}, 'fs', {}, 'fs0', {}, 'nbeat', {}, 'nmatch', {}, ...
    'rr_err', {}, 'rmssd', {}, 'pnn50', {});
switch part
    case 'sim'
        for rep = 1:20
            rng(1000 + rep);
            r = beats(300);
            for fs = [128 250 500 1000]
                rng(2000 + rep);
                ecg = ecg_from_beats(r, fs, round((r(end) + 1)*fs));
                T = detect(ecg, fs, 'ecg');
                [ia, ib] = match(T{1}, r, 0.05);
                R(end+1) = compare(sprintf('sim%02d', rep), fs, fs, T, ia, ib, r); %#ok<AGROW>
            end
        end
    case {'ecg', 'ppg'}
        if ~exist('pop_loadset', 'file'), error('EEGLAB is needed to read ds003838'); end
        d = dir(fullfile(dataPath, 'sub-*'));
        for k = 1:numel(d)
            f = fullfile(d(k).folder, d(k).name, 'ecg', [d(k).name '_task-rest_ecg.set']);
            if ~exist(f, 'file'), continue; end
            [~, E] = evalc('pop_loadset(f)');
            c = find(strcmpi({E.chanlocs.labels}, part), 1);
            if isempty(c), continue; end
            x = double(E.data(c,:));
            R = [R compare_rates(d(k).name, x(isfinite(x)), E.srate, [250 128], part)]; %#ok<AGROW>
        end
    case 'mit'
        recs = strsplit(strtrim(fileread(fullfile(dataPath, 'RECORDS'))));
        recs = setdiff(recs, {'102' '104' '107' '217'});
        for k = 1:numel(recs)
            [x, fs0] = read_mlii(fullfile(dataPath, recs{k}));
            R = [R compare_rates(recs{k}, x, fs0, [180 128], 'ecg')]; %#ok<AGROW>
        end
    otherwise
        error('part must be ''sim'', ''ecg'', ''ppg'' or ''mit''');
end
summary(R, part);
end

function R = compare_rates(name, x, fs0, rates, kind)
% reference: the beats at the native rate, spline-refined
T0 = detect(x, fs0, kind);
ref = T0{3};
R = struct([]);
for fs = rates
    [p, q] = rat(fs/fs0, 1e-6);
    y = resample(x - mean(x), p, q) + mean(x);   % resample pads with zeros
    T = detect(y, fs0*p/q, kind);
    [ia, ib] = match(T{1}, ref, 0.05);
    R = [R compare(name, fs, fs0, T, ia, ib, ref)]; %#ok<AGROW>
end
end

function s = compare(rec, fs, fs0, T, ia, ib, ref)
% interval errors (ms RMS) of the three beat times over consecutive matched
% beats; RMSSD and pNN50 of each, and of the reference on the same intervals
s = struct('rec', rec, 'fs', fs, 'fs0', fs0, 'nbeat', numel(T{1}), 'nmatch', numel(ia), ...
    'rr_err', nan(1,3), 'rmssd', nan(1,4), 'pnn50', nan(1,4));
ok = diff(ia) == 1 & diff(ib) == 1;
if ~any(ok), return; end
rrRef = diff(ref(ib));  rrRef = rrRef(ok);
for m = 1:3
    rr = diff(T{m}(ia));  rr = rr(ok);
    s.rr_err(m) = 1000*sqrt(mean((rr - rrRef).^2));
    [s.rmssd(m), s.pnn50(m)] = hrv(rr);
end
[s.rmssd(4), s.pnn50(4)] = hrv(rrRef);
end

function [rmssd, pnn50] = hrv(rr)
% over successive intervals only (a gap in the matched beats breaks a pair
% only once, at the same place for every beat time)
dd = diff(rr);
rmssd = 1000*sqrt(mean(dd.^2));
pnn50 = 100*mean(abs(dd) > 0.05);
end

function T = detect(x, fs, kind)
% beat times (s, column): peak sample, GET_RR's times, spline extremum
tms = (0:numel(x)-1)/fs*1000; %#ok<NASGU> used in evalc
[~, RR, RR_t, pk, sig] = evalc('get_RR(x, tms, struct(''heart_signal'', kind, ''fs'', fs))');
tb = [RR_t(1) - RR(1); RR_t(:)];            % first beat restored (get_RR drops it)
pk = [round(tb(1)*fs) + 1; pk(:)];
keep = pk > 5 & pk < numel(sig) - 5;        % room for the spline
pk = pk(keep);  tb = tb(keep);
t0 = (pk - 1)/fs;
sgn = 1 - 2*strcmp(kind, 'ppg');            % R-peaks are maxima, PPG onsets minima
u = -1:0.01:1;
t2 = t0;
for k = 1:numel(pk)
    v = sgn * spline(-5:5, sig(pk(k) + (-5:5)), u);
    [~, i] = max(v);
    t2(k) = t0(k) + u(i)/fs;
end
T = {t0, tb, t2};
end

function [ia, ib] = match(a, b, tol)
% one-to-one nearest matches of the times a in b within tol (s)
ia = zeros(0,1);  ib = zeros(0,1);
j = 1;
for i = 1:numel(a)
    while j < numel(b) && abs(b(j+1) - a(i)) < abs(b(j) - a(i)), j = j + 1; end
    if abs(b(j) - a(i)) <= tol && (isempty(ib) || ib(end) ~= j)
        ia(end+1,1) = i;  ib(end+1,1) = j; %#ok<AGROW>
    end
end
end

function summary(R, part)
fprintf('\n== %s\n', part);
for fs = unique([R.fs], 'stable')
    x = R([R.fs] == fs);
    g = [x.nmatch] >= 0.9*[x.nbeat] & [x.nbeat] > 0;
    g = g & arrayfun(@(s) all(isfinite(s.rr_err)), x);
    E = cat(1, x(g).rr_err);  M = cat(1, x(g).rmssd);  P = cat(1, x(g).pnn50);
    dm = M(:,1:3) - M(:,4);  dp = P(:,1:3) - P(:,4);
    fprintf('%4g Hz (from %s Hz): %d recordings, %d left out (< 90%% of the beats matched: %s)\n', fs, ...
        mat2str(unique([x.fs0])), nnz(g), nnz(~g), strjoin({x(~g).rec}, ', '));
    fprintf('  interval error, ms RMS, median [range]: grid %.2f [%.2f-%.2f], get_RR %.2f [%.2f-%.2f], spline %.2f [%.2f-%.2f]\n', ...
        median(E(:,1)), min(E(:,1)), max(E(:,1)), median(E(:,2)), min(E(:,2)), max(E(:,2)), median(E(:,3)), min(E(:,3)), max(E(:,3)));
    fprintf('  RMSSD bias, ms, median / largest:      grid %+.2f / %.2f, get_RR %+.2f / %.2f, spline %+.2f / %.2f (RMSSD median %.1f)\n', ...
        median(dm(:,1)), max(abs(dm(:,1))), median(dm(:,2)), max(abs(dm(:,2))), median(dm(:,3)), max(abs(dm(:,3))), median(M(:,4)));
    fprintf('  pNN50 error, points, median / largest: grid %+.2f / %.2f, get_RR %+.2f / %.2f, spline %+.2f / %.2f (pNN50 median %.1f)\n', ...
        median(dp(:,1)), max(abs(dp(:,1))), median(dp(:,2)), max(abs(dp(:,2))), median(dp(:,3)), max(abs(dp(:,3))), median(P(:,4)));
end
end

function r = beats(dur)
% beat times (s): 70 bpm, LF (0.1 Hz, 30 ms) and HF (0.25 Hz, 40 ms)
% modulation, 2 ms of noise
m = 60/70;  r = 1;
while r(end) < dur + 1
    r(end+1,1) = r(end) + m + 0.030*sin(2*pi*0.1*r(end)) + 0.040*sin(2*pi*0.25*r(end) + 1) + 0.002*randn; %#ok<AGROW>
end
end

function ecg = ecg_from_beats(r, fs, N)
% P, Q, R, S, T waves (Gaussians, R at the beat time), respiration, noise
w = [0.15 -0.10 1.0 -0.25 0.30]*1000;  mu = [-0.18 -0.025 0 0.025 0.25];  sd = [0.025 0.008 0.010 0.010 0.040];
tt = (0:N-1)/fs;  ecg = zeros(1, N);
for k = 1:numel(r)
    i0 = max(1, round((r(k)-0.4)*fs));  i1 = min(N, round((r(k)+0.6)*fs));
    seg = tt(i0:i1) - r(k);
    for j = 1:5
        ecg(i0:i1) = ecg(i0:i1) + w(j)*exp(-(seg - mu(j)).^2/(2*sd(j)^2));
    end
end
ecg = ecg + 80*sin(2*pi*0.25*tt) + 20*randn(1, N) + 500;
end

function [x, fs] = read_mlii(base)
% MLII lead (or the first signal) of a MIT-BIH record, format 212, in mV
h = strsplit(strtrim(fileread([base '.hea'])), newline);
h0 = sscanf(h{1}, '%*s %d %f %d');
nsig = h0(1);  fs = h0(2);
lines = strtrim(h(2:nsig+1));
c = find(endsWith(lines, 'MLII'), 1);
if isempty(c), c = 1; end
f = strsplit(lines{c});
gain = str2double(regexp(f{3}, '^[\d.]+', 'match', 'once'));  zero = str2double(f{5});
fid = fopen([base '.dat'], 'r');
B = fread(fid, [3 inf], 'uint8=>double');
fclose(fid);
s = [B(1,:) + 256*mod(B(2,:), 16); B(3,:) + 256*floor(B(2,:) / 16)];
s(s > 2047) = s(s > 2047) - 4096;
x = (s(c,:) - zero) / gain;
end
