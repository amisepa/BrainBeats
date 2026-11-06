function [R, U] = compare_ectopic_rules(dataPath, nRep)
% COMPARE_ECTOPIC_RULES - Ectopic-beat step of CLEAN_RR on real RR series
% with injected artifacts: no ectopic step ('none'), the fixed rule ('ratio':
% an interval 20% shorter than the local median followed by one 20% longer)
% and the ectopic beats of the beat classification of Lipponen & Tarvainen
% (2019), each moved to the middle of its neighbors: all of them
% ('adaptive_0': 'ectopic_min_shift' 0), or only those the move shifts by
% at least 10% of their interval ('adaptive', 'ectopic_min_shift' 0.1).
%
% The R-peaks GET_RR finds in each recording's ECG are the truth. Each
% artifact type is injected in its own copies of the series, 6 per copy at
% random beats >= 15 beats apart (nRep copies per recording and type, each
% with its own seed); the rules run on the same copies:
%   premature .70/.75/.80 - a beat moved earlier, to .70-.80 x its interval,
%                           the next beat unchanged (full compensatory pause)
%   missed                - a beat removed
%   extra                 - a false beat at a random 20-80% of an interval
% A premature or extra beat is detected when CLEAN_RR removes it, a missed
% beat when it inserts a beat in the gap. False positives (FP): changes away
% from the artifacts (outside the two true intervals around each): true
% beats removed and beats inserted; new FP: those the same rule does not
% make on the series as recorded (caused by the injected artifacts). RMSSD
% and SDNN errors: after cleaning vs the series as recorded.
% The series as recorded are cleaned too: on clean recordings every change
% is a false positive. sub-038 is not injected: it has real ectopic beats
% (5 abnormal beats on its ECG at 3.9, 42.6, 124.4, 166.0 and 201.1 s,
% checked 2026-11-03) and many missed ones.
%
% Usage:
%   R = compare_ectopic_rules(dataPath)        % 10 copies per type
%   R = compare_ectopic_rules(dataPath, nRep)
%   R = compare_ectopic_rules(files, nRep)     % other recordings
%
% Inputs:
%   dataPath - folder with sub-033, 034, 035, 036, 038, 040, 062 and
%              070_rest.set of ds003838 made by PREP_DS003838 (missing ones
%              are skipped); the two datasets of sample_data are added.
%              Or a cell array of .set files with an 'ECG' channel.
%   nRep     - copies per recording and artifact type (default 10)
%
% Outputs:
%   R - table, one row per copy and rule: file, type, rep, rule, detected,
%       injected, fp, fp_new, rmssd_err, sdnn_err (ms, cleaned - true)
%   U - table, one row per recording and rule, the series as recorded:
%       file, rule, removed, inserted, removed_t, rmssd_err, sdnn_err
%
% Results (2026-11-05; 10 copies x 6 artifacts per recording and type).
% Three sets, made with PREP_DS003838: design (the default: the two sample
% datasets, sub-033 to 036, 040, 062, 070), held-out (sub-039, 041 to 046,
% 048 to 056) and fresh (sub-057 to 061, 063, 067 to 069, 071 to 076).
% Detected % (design / held-out / fresh), new FP per copy, |RMSSD error| ms:
%                  ratio                  adaptive_0             adaptive
%   premature .70  86.5/88.8/83.4         99.1/97.8/98.6         99.1/97.8/98.6
%     new FP       0/.06/0                .12/.16/.12            .06/.12/.03
%     RMSSD        5.7/6.5/8.9            0.6/2.9/6.6            0.6/2.9/6.5
%   premature .75  63.9/59.5/65.1         96.9/95.6/97.1         96.9/95.6/97.1
%     new FP       0/.05/.03              .08/.16/.13            .06/.11/.05
%     RMSSD        12.9/15.0/12.4         0.9/3.3/7.1            0.9/3.2/6.9
%   premature .80  12.6/14.7/13.0         93.1/92.1/96.1         93.1/92.1/96.4
%     new FP       0/.03/.03              .03/.22/.19            .01/.16/.05
%     RMSSD        22.7/19.7/19.6         1.4/3.4/6.9            1.3/3.4/6.8
%   missed         100/99.7/99.7          100/99.7/99.8          100/99.7/99.8
%     new FP       .09/.10/.07            .20/.25/.33            .20/.21/.19
%   extra          79.1/77.5/83.0         86.7/83.9/90.8         86.7/83.9/90.8
%     new FP       0/.03/.02              .10/.11/.19            .10/.11/.14
%     RMSSD        9.3/12.2/9.1           6.3/9.9/8.0            6.2/9.9/7.9
% Uncleaned, the RMSSD error was 48/47/48 ms (premature .70), 26/25/26
% (.80), 119/116/122 (missed), 56/55/57 (extra); 'none' (Steps 3-5 only)
% detected 1-26% of the premature beats. Errors on the held-out and fresh
% sets include real artifacts of the recordings. Per recording, 'ratio'
% detected 2.3-6.0 of 6 premature beats at .70 and 0.1-1.6 at .80,
% 'adaptive' 4.7-6.0 and 2.8-6.0 (fewest where RMSSD is high: sub-055,
% 134 ms, 2.8 at .80; sub-054, 105 ms, 3.9).
% As recorded (12,800 beats; ectopic step only): 'ratio' replaced 1 beat
% (a false detection, sub-061 252.75 s, 632/312 ms; 'adaptive' too);
% 'adaptive' 5 more: 2 false detections (sub-068 93.15 s, an electrode pop;
% sub-071 238.10 s, a spike), 1 sinus beat (sub-069 112.87 s, RMSSD 147
% ms), 2 unclear on the ECG (sub-033 149.48 s, noisy; sub-069 251.94 s, in
% an artifact); 'adaptive_0' 8 more, all sinus beats. sub-038 (5 abnormal
% beats on its ECG): 'ratio' replaced 4 of them, 'adaptive' 3 (not 166.0 s,
% next to a beat Step 3 removes) and 3 other beats of this noisy recording.
% 'adaptive' is the default: it detects far more premature beats, with
% more false positives ('ratio' replaced no sinus beat).
% How the rule was set: the full classification (long/short, extra and
% missed beats replaced too, as NeuroKit2 does) changed 39 beats in 6 of
% the 9 design recordings as recorded, sinus beats on the ECG, and was
% dropped. The 10% minimum shift was chosen after the 4 sinus beats the
% published rule replaced in the design and held-out sets; the fresh set
% was then run blind. The choice of the beat to move (smaller successive
% differences, instead of the larger |mRR| neighbor, which moved the beat
% after the pause when the interval the premature beat replaced was long)
% came after that run: detection -0.1 to +1.8 points; as recorded, the
% sinus beat of sub-069 was replaced too (4 -> 5 beats on the fresh set).
%
% Copyright (C) - Cedric Cannard, 2026

if nargin < 2 || isempty(nRep), nRep = 10; end
bb = fileparts(fileparts(mfilename('fullpath')));
if ~exist('pop_loadset','file'), eeglab; close; end
p = strsplit(path, pathsep);
old = p(contains(p,'BrainBeats','IgnoreCase',true) & ~startsWith(p,bb));
if ~isempty(old), rmpath(old{:}); end
addpath(fullfile(bb,'functions'));

if iscell(dataPath)
    files = dataPath;
else
    files = {fullfile(bb,'sample_data','dataset.set'), fullfile(bb,'sample_data','sub-047_rest.set')};
    for s = {'033' '034' '035' '036' '038' '040' '062' '070'}
        f = fullfile(dataPath, sprintf('sub-%s_rest.set', s{1}));
        if exist(f,'file'), files{end+1} = f; else, warning('%s not found: skipped', f); end %#ok<AGROW>
    end
end
rules = {'none' 'ratio' 'adaptive_0' 'adaptive'};
ruleArgs = {{'ectopic_rule','none'}, {'ectopic_rule','ratio'}, ...
    {'ectopic_rule','adaptive','ectopic_min_shift',0}, {'ectopic_rule','adaptive'}};
types = {'premature .70' 'premature .75' 'premature .80' 'missed' 'extra'};
prem  = [0.70 0.75 0.80 NaN NaN];
nInj  = 6;
rmssd = @(x) sqrt(mean(diff(1000*x).^2));
sdnn  = @(x) 1000*std(x);

R = table();
U = table();
for k = 1:numel(files)
    [fp, fn, ext] = fileparts(files{k});
    EEG = pop_loadset('filename',[fn ext],'filepath',fp);
    iE = find(strcmpi({EEG.chanlocs.labels},'ECG'), 1);
    [RR, RR_t] = get_RR(double(EEG.data(iE,:)), EEG.times, struct('heart_signal','ecg','fs',EEG.srate));
    fs = EEG.srate;
    bt = [RR_t(1) - RR(1); RR_t(:)];           % all beats, the first one included
    nb = numel(bt);
    fprintf('\n%s: %d beats, RMSSD %.1f ms, SDNN %.1f ms\n', fn, nb, rmssd(diff(bt)), sdnn(diff(bt)));

    % the series as recorded
    goneR = cell(1, numel(rules));  synR = goneR;
    for r = 1:numel(rules)
        [nn, nt, ii] = run_clean(bt, fs, ruleArgs{r});
        goneR{r} = bt(2:end);  goneR{r} = goneR{r}(~ismembertol(goneR{r}, nt(~ii), 1e-6, 'DataScale', 1));
        synR{r} = nt(ii);
        gone = setdiff(round(bt(2:end)*1e4), round(nt(~ii)*1e4)) / 1e4;
        fprintf('  as recorded, %-12s %d removed, %d inserted; removed at %s s\n', rules{r}, numel(gone), sum(ii), mat2str(gone', 5));
        U = [U; table({fn}, rules(r), numel(gone), sum(ii), {gone}, rmssd(nn) - rmssd(diff(bt)), sdnn(nn) - sdnn(diff(bt)), ...
            'VariableNames', {'file','rule','removed','inserted','removed_t','rmssd_err','sdnn_err'})]; %#ok<AGROW>
    end
    if contains(fn, '038'), continue; end

    % injected copies
    for ty = 1:numel(types)
        for rep = 1:nRep
            rng(100000*k + 1000*ty + rep);
            edges = round(linspace(25, nb-25, nInj+1));
            pos = arrayfun(@(e) randi([edges(e)+7, edges(e+1)-8]), 1:nInj)';   % beat before each artifact
            % site of each artifact: between the true beats around it
            site = [bt(pos) bt(pos+2)];
            in_site = @(t) any(t(:) > site(:,1)' & t(:) < site(:,2)', 2);
            bt2 = bt; te = NaN(nInj,1);
            switch types{ty}(1:5)
                case 'prema'
                    te = bt(pos) + prem(ty) * (bt(pos+1) - bt(pos));
                    bt2(pos+1) = te;
                case 'misse'
                    bt2(pos+1) = [];
                case 'extra'
                    te = bt(pos) + (0.2 + 0.6*rand(nInj,1)) .* (bt(pos+1) - bt(pos));
                    bt2 = sort([bt2; te]);
            end
            truth = bt(2:end);
            truth = truth(~in_site(truth));                 % true beats away from the artifacts
            for r = 1:numel(rules)
                [nn, nt, ii] = run_clean(bt2, fs, ruleArgs{r});
                real = nt(~ii); syn = nt(ii);
                if strcmp(types{ty}, 'missed')
                    det = sum(arrayfun(@(j) any(syn > site(j,1) & syn < site(j,2)), 1:nInj));
                else
                    det = sum(~ismembertol(te, real, 1e-6, 'DataScale', 1));
                end
                % false positives: changes away from the artifacts; new ones:
                % not made by the same rule on the series as recorded
                lost = truth(~ismembertol(truth, real, 1e-6, 'DataScale', 1));
                synOut = syn(~in_site(syn));
                fpos = numel(lost) + numel(synOut);
                fnew = sum(~ismembertol(lost, goneR{r}, 1e-6, 'DataScale', 1)) + ...
                    sum(~ismembertol(synOut, synR{r}, 0.05, 'DataScale', 1));
                R = [R; table({fn}, types(ty), rep, rules(r), det, nInj, fpos, fnew, ...
                    rmssd(nn) - rmssd(diff(bt)), sdnn(nn) - sdnn(diff(bt)), ...
                    'VariableNames', {'file','type','rep','rule','detected','injected','fp','fp_new','rmssd_err','sdnn_err'})]; %#ok<AGROW>
            end
            if rep == 1   % error of the injected series without any cleaning, for reference
                R = [R; table({fn}, types(ty), rep, {'raw'}, NaN, nInj, NaN, NaN, ...
                    rmssd(diff(bt2)) - rmssd(diff(bt)), sdnn(diff(bt2)) - sdnn(diff(bt)), ...
                    'VariableNames', R.Properties.VariableNames)]; %#ok<AGROW>
            end
        end
    end
end

%% Summary
inj = unique(R.file, 'stable');
fprintf('\n%d recordings injected (%d copies per type, %d artifacts per copy)\n', numel(inj), nRep, nInj);
fprintf('%-14s %-13s %9s %7s %9s %9s %12s %11s\n', 'type', 'rule', 'detected', '%', 'FP/copy', 'new FP', '|dRMSSD| ms', '|dSDNN| ms');
for ty = 1:numel(types)
    for r = [{'raw'} rules]
        s = strcmp(R.type, types{ty}) & strcmp(R.rule, r{1});
        if strcmp(r{1}, 'raw')
            fprintf('%-14s %-13s %9s %7s %9s %9s %12.1f %11.1f\n', types{ty}, 'uncleaned', '', '', '', '', mean(abs(R.rmssd_err(s))), mean(abs(R.sdnn_err(s))));
            continue
        end
        fprintf('%-14s %-13s %9s %7.1f %9.2f %9.2f %12.1f %11.1f\n', types{ty}, r{1}, sprintf('%d/%d', sum(R.detected(s)), sum(R.injected(s))), ...
            100*sum(R.detected(s))/sum(R.injected(s)), mean(R.fp(s)), mean(R.fp_new(s)), mean(abs(R.rmssd_err(s))), mean(abs(R.sdnn_err(s))));
    end
end
fprintf('\nPremature beats detected per recording (of %d per copy, mean over copies):\n', nInj);
for ty = 1:3
    for r = 2:numel(rules)
        m = cellfun(@(f) mean(R.detected(strcmp(R.file,f) & strcmp(R.type,types{ty}) & strcmp(R.rule,rules{r}))), inj);
        fprintf('  %-14s %-13s %s\n', types{ty}, rules{r}, sprintf('%.1f ', m));
    end
end
fprintf('  (recordings: %s)\n', strjoin(inj', ', '));
fprintf(['\nAs recorded (sub-038 excluded): beats removed + inserted; beats removed that ''none'' keeps ' ...
    '(the ectopic step: on clean recordings, false positives); |dRMSSD|, |dSDNN|\n']);
s0 = strcmp(U.rule, 'none') & ~contains(U.file, '038');
for r = 1:numel(rules)
    s = strcmp(U.rule, rules{r}) & ~contains(U.file, '038');
    extra = cellfun(@(a, b) numel(setdiff(a, b)), U.removed_t(s), U.removed_t(s0));
    fprintf('  %-13s %3d changes; ectopic step %d beats in %d of %d recordings; |dRMSSD| %.2f ms, |dSDNN| %.2f ms\n', rules{r}, ...
        sum(U.removed(s) + U.inserted(s)), sum(extra), sum(extra > 0), sum(s), ...
        mean(abs(U.rmssd_err(s))), mean(abs(U.sdnn_err(s))));
end
fprintf('compare_ectopic_rules done\n');
end

function [nn, nt, ii] = run_clean(bt, fs, args)
RR = diff(bt); RR_t = bt(2:end); pk = round(RR_t*fs) + 1;
[~, nn, nt, ~, ~, ii] = evalc('clean_rr(RR_t, RR, [], pk, args{:})');
end
