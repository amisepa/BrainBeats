function R = check_rr_alternation(dataPath)
% CHECK_RR_ALTERNATION - RR_ALTERNATION on resting recordings of healthy
% adults, after GET_RR and CLEAN_RR as BRAINBEATS_PROCESS runs them: how
% far below the 10% warning they stay (the arrhythmic side is in
% VALIDATE_ECTOPIC_MITBIH).
%
% Usage:
%   R = check_rr_alternation(dataPath)
%
% Input:
%   dataPath - folder of ds003838 sub-*_rest.set files made by
%              PREP_DS003838; the two datasets of sample_data are added.
%              Or a cell array of .set files with an 'ECG' channel.
%
% Output:
%   R - table: file, beats, alternation (%), rmssd (ms)
%
% Results (2026-11-06): 66 recordings (the two sample datasets and 64 of
% ds003838): median 0%, 95th percentile 1.0%, maximum 6.9% (sub-055, RMSSD
% 134 ms), then 3.8% (sub-038, noisy, with real ectopic beats). No warning.
%
% Copyright (C) - Cedric Cannard, 2026

bb = fileparts(fileparts(mfilename('fullpath')));
if ~exist('pop_loadset','file'), eeglab; close; end
p = strsplit(path, pathsep);
old = p(contains(p,'BrainBeats','IgnoreCase',true) & ~startsWith(p,bb));
if ~isempty(old), rmpath(old{:}); end
addpath(fullfile(bb,'functions'));

if iscell(dataPath)
    files = dataPath;
else
    d = dir(fullfile(dataPath, 'sub-*_rest.set'));
    files = [{fullfile(bb,'sample_data','dataset.set')}, fullfile({d.folder}, {d.name})];
    if ~any(contains({d.name}, 'sub-047')), files{end+1} = fullfile(bb,'sample_data','sub-047_rest.set'); end
end
R = table();
for k = 1:numel(files)
    [fp, fn, ext] = fileparts(files{k});
    EEG = pop_loadset('filename',[fn ext],'filepath',fp,'verbose','off');
    iE = find(strcmpi({EEG.chanlocs.labels},'ECG'), 1);
    [~, RR, RR_t, pk, sigF] = evalc('get_RR(double(EEG.data(iE,:)), EEG.times, struct(''heart_signal'',''ecg'',''fs'',EEG.srate))');
    [~, nn] = evalc('clean_rr(RR_t, RR, sigF(pk), pk, ''ecg_signal'', sigF, ''fs'', EEG.srate)');
    a = 100 * rr_alternation(nn);
    R = [R; table({fn}, numel(RR), a, 1000*sqrt(mean(diff(nn).^2)), ...
        'VariableNames', {'file','beats','alternation','rmssd'})]; %#ok<AGROW>
    fprintf('%-16s %4d beats  alternation %5.1f%%  RMSSD %6.1f ms\n', fn, numel(RR), a, R.rmssd(end));
end
fprintf('\n%d recordings: alternation median %.1f%%, 95th percentile %.1f%%, maximum %.1f%%; %d above 10%%\n', ...
    height(R), median(R.alternation), prctile(R.alternation, 95), max(R.alternation), sum(R.alternation > 10));
end
