function compare_hrv_windows(files, winLengths)
% COMPARE_HRV_WINDOWS - HRV band powers with LF and HF from the same 120-s
% windows (GET_HRV_FEATURES since v1.7) vs the v1.6 windows (5 cycles of each
% band's lowest frequency: 125 s for LF, 34 s for HF), on recordings and on
% synthetic NN series of known modulation. Stops with an error if a check
% fails.
%
% Usage:
%   compare_hrv_windows                          % the two sample datasets
%   compare_hrv_windows(files)                   % datasets with an 'ECG' channel
%   compare_hrv_windows(files, [60 90 180])      % also other LF/HF windows
%
% Inputs:
%   files      - EEGLAB datasets with an 'ECG' channel (default:
%                sample_data/dataset.set and sample_data/sub-047_rest.set)
%   winLengths - LF/HF window lengths (s) to compare (default: none):
%                copies of GET_HRV_FEATURES with another window, in tempdir
%
% The v1.6 function is read from git (commit cb963bd) and written to
% tempdir as GET_HRV_FEATURES_V16. Each recording: beats from the ECG
% (brainbeats_process, HRV time domain only), then both functions on the
% same NN series. Synthetic series: mean NN 60/HR s, sinusoids of known
% amplitude A (band power A^2/2) and white noise of SD s (band power
% 2*s^2*(60/HR)*bandwidth), 100 series per case, fresh seeds per case.
%
% Results (2026-11-05; 6 recordings and winLengths [60 90 180]: 151 s)
% Synthetic, bias of the mean (SD across series), v1.6 -> current:
%   0.1 + 0.25 Hz (30 and 20 ms), noise 15 ms, 70 bpm, 300 s (expected LF
%   492, HF 296 ms^2): Lomb-Scargle LF +0.5% (8.9%) -> -0.1% (7.8%), HF
%   +2.6% (8.8%) -> +0.8% (8.6%); welch LF +1.0% (10.9%) -> +0.1% (8.9%),
%   HF -2.3% (9.9%) -> -2.3% (9.7%). At 186 s and at 90 bpm: all within 2%.
%   Peaks near the LF/HF border (0.13, 0.17 Hz): HF SD 15.0 -> 9.5%
%   (Lomb-Scargle), 33.4 -> 10.3% (welch).
%   1/f-like (9 sinusoids, 0.01-0.33 Hz, 240 s): Lomb-Scargle LF +10.9 ->
%   +10.2%, HF +14.2 -> +4.9%; welch LF +1.0 -> +1.5%, HF -1.9 -> -2.2%.
%   Normalized LF/HF 1.02 -> 1.79 (ms^2: 1.67 -> 1.80).
%   Other windows, LF on the 1/f-like series (Lomb-Scargle / welch): 60 s
%   +24/+28%, 90 s +13/+7%, 120 s +10/+2%, 180 s +8/0%; with a VLF
%   sinusoid (0.025 Hz, 50 ms): 60 s +11/+33%, 90 s +11/+2%, 120 s +7/-1%,
%   180 s +5/0%. On 186 s, 180 s (one window) had the largest welch SD
%   (LF 13.5%, HF 15.6% vs 12.1%, 13.7% with 120 s).
%   Consistency: normalized LF/HF = ms^2 LF/HF (1.92); 'hrv_norm' LF =
%   LF/(LF+HF). 150 s: LF and HF; 90 s: HF only (whole recording, LF
%   skipped with a warning); 50 s: NaN.
% Recordings (ECG; beats from the pipeline of that day), LF, HF (ms^2),
% LF/HF, LF/(LF+HF), v1.6 -> current:
%                  Lomb-Scargle                          welch
%   dataset.set    1071 -> 1071, 1896 -> 1525, .57 -> .70  1442 -> 1313, 2268 -> 1850, .64 -> .71
%   sub-047        297 -> 335, 579 -> 559, .51 -> .60      296 -> 282, 538 -> 545, .55 -> .52
%   sub-062        201 -> 189, 291 -> 285, .69 -> .66      167 -> 151, 293 -> 306, .57 -> .49
%   sub-033        7105 -> 7188, 3174 -> 2800, 2.24 -> 2.57  5844 -> 6158, 2295 -> 2767, 2.55 -> 2.23
%   sub-040        249 -> 216, 157 -> 143, 1.59 -> 1.51    179 -> 197, 142 -> 138, 1.26 -> 1.42
%   sub-070        760 -> 649, 1380 -> 1290, .55 -> .50    506 -> 513, 1344 -> 1322, .38 -> .39
%   LF/(LF+HF), Lomb-Scargle: .361 -> .412, .339 -> .375, .408 -> .398,
%   .691 -> .720, .614 -> .602, .355 -> .335. The v1.6 HF was higher on all
%   6 (2-24%, Lomb-Scargle). Normalized Lomb-Scargle LF/HF (ms^2 in
%   brackets): .48 -> .69 (.70), .48 -> .59 (.60), .54 -> .66 (.66), 2.03 ->
%   2.60 (2.57), .83 -> 1.50 (1.51), .44 -> .45 (.50).
%   LF windows: 2, 3, 2, 1, 1, 2 (v1.6) -> 3, 4, 3, 2, 3, 3; a fixed 90-s
%   step would have left out 3-35% of the recordings.
%
% Copyright (C) - Cedric Cannard, 2026

% Paths: EEGLAB, then this repository instead of any other BrainBeats copy
if ~exist('pop_loadset','file'), eeglab; close; end
bb = fileparts(fileparts(mfilename('fullpath')));
p = strsplit(path, pathsep);
old = p(contains(p,'BrainBeats','IgnoreCase',true) & ~startsWith(p,bb));
if ~isempty(old), rmpath(old{:}); end
addpath(bb, fullfile(bb,'functions'));
if nargin < 1 || isempty(files)
    files = {fullfile(bb,'sample_data','dataset.set'), fullfile(bb,'sample_data','sub-047_rest.set')};
end
if nargin < 2, winLengths = []; end
files = cellstr(files);

% Reference: the v1.6 function from git, and copies of the current one with
% other LF/HF windows
d = fullfile(tempdir, 'brainbeats_hrvwin');
if ~exist(d,'dir'), mkdir(d); end
[status, src] = system(sprintf('git -C "%s" show cb963bd:functions/get_hrv_features.m', bb));
assert(status == 0, 'git show failed: %s', src)
writefun(d, 'get_hrv_features_v16', strrep(src, 'function [HRV, params] = get_hrv_features(', ...
    'function [HRV, params] = get_hrv_features_v16('));
cur = fileread(fullfile(bb,'functions','get_hrv_features.m'));
key = 'lfhfWin = 120;';
assert(numel(strfind(cur, key)) == 1, 'LF/HF window line not found')
for w = winLengths(:)'
    v = strrep(cur, key, sprintf('lfhfWin = %g;', w));
    writefun(d, sprintf('get_hrv_features_w%g', w), strrep(v, 'function [HRV, params] = get_hrv_features(', ...
        sprintf('function [HRV, params] = get_hrv_features_w%g(', w)));
end
addpath(d);
funs = [{'get_hrv_features_v16' 'get_hrv_features'} arrayfun(@(w) sprintf('get_hrv_features_w%g', w), winLengths(:)', 'UniformOutput', false)];
P = struct('hrv_time',false,'hrv_frequency',true,'hrv_nonlinear',false);

%% 1. Synthetic: known band powers
% case: name, duration (s), heart rate (bpm), sinusoids [f (Hz) A (ms)],
% noise SD (ms), checked (current function within 10% of the expected
% powers on average; not with a large VLF component, which leaks into LF
% with Lomb-Scargle: no taper)
cases = {'0.1+0.25 Hz'         300 70 [.10 30; .25 20] 15 true; ...
         '0.13+0.17 Hz'        300 70 [.13 30; .17 20] 15 true; ...
         '+VLF 0.025 Hz 50 ms' 300 70 [.025 50; .10 30; .25 20] 15 false; ...
         '0.1+0.25 Hz'         186 70 [.10 30; .25 20] 15 true; ...
         '0.1+0.25 Hz'         300 90 [.10 30; .25 20] 15 true; ...
         '1/f-like, 9 sinusoids' 240 75 [.01 40; .02 30; .03 25; .06 15; .09 20; .12 12; .2 10; .27 15; .33 8] 10 false};
nSer = 100;
for c = 1:size(cases,1)
    [name, T, hr, comps, sd, checked] = cases{c,:};
    inLF = comps(:,1) >= .04 & comps(:,1) < .15;
    inHF = comps(:,1) >= .15 & comps(:,1) <= .40;
    expLF = sum(comps(inLF,2).^2/2) + 2*sd^2*(60/hr)*.11;
    expHF = sum(comps(inHF,2).^2/2) + 2*sd^2*(60/hr)*.25;
    fprintf('\nSYNTH %s, %g s, %g bpm: expected LF %.0f, HF %.0f ms^2, LF/HF %.2f\n', name, T, hr, expLF, expHF, expLF/expHF);
    for spec = {'LombScargle' 'welch' 'LombScargle_norm'}
        P.hrv_spec = spec{1};
        R = nan(nSer, numel(funs), 2);
        for s = 1:nSer
            [NN, t] = sim_nn(T, hr, comps, sd, 1e5*c + s);
            for f = 1:numel(funs)
                H = run_quiet(funs{f}, NN, t, P);
                R(s,f,:) = [H.frequency.lf H.frequency.hf];
            end
        end
        for f = 1:numel(funs)
            lf = R(:,f,1); hf = R(:,f,2);
            if strcmp(spec{1}, 'LombScargle_norm')
                fprintf('SYNTH   %-26s %-17s LF %.3f  HF %.3f  LF/HF %.2f  LF/(LF+HF) %.3f\n', funs{f}, spec{1}, ...
                    mean(lf), mean(hf), mean(lf)/mean(hf), mean(lf./(lf+hf)));
            else
                fprintf('SYNTH   %-26s %-17s LF %4.0f (%+5.1f%%, SD %4.1f%%)  HF %4.0f (%+5.1f%%, SD %4.1f%%)  LF/HF %.2f  LF/(LF+HF) %.3f\n', ...
                    funs{f}, spec{1}, mean(lf), 100*(mean(lf)/expLF - 1), 100*std(lf)/expLF, ...
                    mean(hf), 100*(mean(hf)/expHF - 1), 100*std(hf)/expHF, mean(lf)/mean(hf), mean(lf./(lf+hf)));
            end
        end
        % Check (current function): band powers within 10% of the expected
        % ones, on average
        if checked && ~strcmp(spec{1}, 'LombScargle_norm')
            bias = [mean(R(:,2,1))/expLF mean(R(:,2,2))/expHF] - 1;
            assert(all(abs(bias) < .10), '%s, %s: LF and HF %+.1f%% and %+.1f%% of the expected powers', ...
                name, spec{1}, 100*bias(1), 100*bias(2))
        end
    end
end

%% 2. Consistency: normalized units and 'hrv_norm' (current function)
[NN, t] = sim_nn(300, 70, [.10 30; .25 20], 15, 7e5);
P.hrv_spec = 'LombScargle';
H = run_quiet('get_hrv_features', NN, t, P);
Pn = P; Pn.hrv_spec = 'LombScargle_norm';
Hn = run_quiet('get_hrv_features', NN, t, Pn);
Pr = P; Pr.hrv_norm = true;
Hr = run_quiet('get_hrv_features', NN, t, Pr);
fprintf(['\nCONSIST LF/HF %.2f in ms^2, %.2f normalized; hrv_norm LF %.4f = LF/(LF+HF) %.4f, ' ...
    'HF %.4f = HF/(LF+HF) %.4f\n'], H.frequency.lfhf, Hn.frequency.lfhf, Hr.frequency.lf, ...
    H.frequency.lf/(H.frequency.lf + H.frequency.hf), Hr.frequency.hf, H.frequency.hf/(H.frequency.lf + H.frequency.hf));
assert(abs(Hn.frequency.lfhf/H.frequency.lfhf - 1) < .05, 'Normalized LF/HF differs from the ms^2 one')
assert(abs(Hr.frequency.lf - H.frequency.lf/(H.frequency.lf + H.frequency.hf)) < 1e-3, 'hrv_norm is not LF/(LF+HF)')
assert(Hr.frequency.lfhf == H.frequency.lfhf, 'hrv_norm changed LF/HF')
assert(all(diff(H.frequency.pwr_freqs) > 0) && H.frequency.pwr_freqs(1) >= .04 && H.frequency.pwr_freqs(end) <= .40, ...
    'Exported spectrum: not one 0.04-0.40 Hz spectrum')

%% 3. Short recordings: 150 s (LF and HF), 90 s (HF from the whole recording, LF skipped), 50 s (none)
for T = [150 90 50]
    [NN, t] = sim_nn(T, 70, [.10 30; .25 20], 15, 8e5 + T);
    [H, txt] = run_quiet('get_hrv_features', NN, t, P);
    F = H.frequency;
    fprintf('SHORT %g s: LF %s, HF %s, LF/HF %s\n', T, num2str(getf(F,'lf')), num2str(getf(F,'hf')), num2str(getf(F,'lfhf')));
    switch T
        case 150
            assert(isfield(F,'lf') && isfield(F,'hf') && isfield(F,'lfhf') && ~isnan(F.lf) && ~isnan(F.hf), '150 s: LF or HF missing')
        case 90
            assert(~isfield(F,'lf') && isfield(F,'hf') && ~isnan(F.hf) && ~isfield(F,'lfhf'), '90 s: LF present or HF missing')
            assert(contains(txt, 'too short for estimating LF') && contains(txt, 'HF power from the whole recording'), '90 s: no warning')
        case 50
            assert(all(isnan([F.ulf F.vlf F.lf F.hf])) && isempty(F.pwr), '50 s: band powers not NaN')
    end
end

%% 4. Recordings: v1.6 vs current on the same NN series
for k = 1:numel(files)
    [fp, fn, ext] = fileparts(files{k});
    EEG = pop_loadset('filename',[fn ext],'filepath',fp);
    EEG = pop_select(EEG,'channel',{'ECG'});
    EEG = brainbeats_process(EEG,'analysis','features','eeg','off','heart_signal','ECG','heart_channels',{'ECG'}, ...
        'hrv_features',{'time'},'vis_cleaning',0,'vis_outputs',0,'gong',0,'save',0);
    NN = EEG.brainbeats.preprocessings.NN(:); t = EEG.brainbeats.preprocessings.NN_times(:);
    span = t(end) - t(1) + NN(1);
    fprintf('\nREC %s: %.0f s, %d NN, %.0f bpm\n', fn, span, numel(NN), 60/mean(NN));
    for spec = {'LombScargle' 'welch' 'LombScargle_norm'}
        P.hrv_spec = spec{1};
        for f = 1:numel(funs)
            H = run_quiet(funs{f}, NN, t, P);
            F = H.frequency;
            fprintf('REC   %-26s %-17s LF %8.4g  HF %8.4g  LF/HF %5.2f  LF/(LF+HF) %.3f\n', funs{f}, spec{1}, ...
                F.lf, F.hf, F.lfhf, F.lf/(F.lf + F.hf));
        end
    end
end
fprintf('\ncompare_hrv_windows: all checks passed\n');
end

% -------------------------------------------------------------------------
function [H, txt] = run_quiet(fun, NN, t, P)
% Run an HRV function with its printed output captured
[txt, H] = evalc(sprintf('%s(NN, t, P)', fun));
end

function v = getf(S, f)
if isfield(S, f), v = S.(f); else, v = 'absent'; end
end

function writefun(d, name, src)
fid = fopen(fullfile(d, [name '.m']), 'w');
fwrite(fid, src);
fclose(fid);
clear(name);
end

function [NN, t] = sim_nn(T, hr, comps, sd, seed)
% NN series (s) and beat times (s) of a T-s recording: mean 60/hr s plus
% sinusoids comps = [f (Hz) A (ms)] with random phases, taken at the start
% of each interval, plus white noise of SD sd (ms)
rng(seed);
m = 60/hr;
phi = 2*pi*rand(size(comps,1), 1);
NN = zeros(ceil(1.5*T/m) + 10, 1); t = NN;
tk = 0.8;   % first R-peak (s)
k = 0;
while true
    v = m + sum(comps(:,2)/1000 .* sin(2*pi*comps(:,1)*tk + phi)) + sd/1000*randn;
    if tk + v > T, break; end
    k = k + 1;
    tk = tk + v;
    NN(k) = v; t(k) = tk;
end
NN = NN(1:k); t = t(1:k);
end
