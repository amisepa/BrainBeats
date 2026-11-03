function T = compare_ppg_transit(files, outfile)
% COMPARE_PPG_TRANSIT - How well the PPG pulse arrival time (PAT) correction
% recovers the ECG-locked HEP, for several datasets with EEG, ECG and PPG.
%
% For each dataset, the HEP (no EEG cleaning, all channels, -300 to 600 ms)
% is computed from the ECG R-peaks (the reference), and from the PPG beats
% shifted back by: nothing ('off'), the PAT estimated from the ECG ('auto'
% with the ECG in the file), the PAT estimated from the EEG cardiac field
% artifact ('eeg', ECG removed), a fixed 250 ms (a common literature value
% for the pulse onset at the finger), and, with 'ppg_detect_mode','peaks'
% (systolic peaks), the PAT from the ECG. A run that errors gets NaN.
% Each PPG-locked HEP is compared with the ECG-locked one.
%
% Usage:
%   T = compare_ppg_transit({'sub-01.set' 'sub-02.set'})
%   T = compare_ppg_transit(files, 'ppg_transit.csv')
%
% Output:
%   T - table, one row per dataset and correction: file, correction, pat
%       (ms), pat_iqr (ms, beat-to-beat spread of the ECG-based PAT), ratio
%       (EEG estimate: GFP peak/median), r (mean over channels of the
%       correlation with the ECG-locked HEP, 0 to 600 ms), offset (ms,
%       latency of the QRS field GFP peak minus that of the ECG-locked HEP,
%       -100 to 150 ms), peak (QRS field GFP peak relative to the ECG-locked
%       one: < 1 when beat-to-beat PAT variation smears it), epochs
%
% Copyright (C) - Cedric Cannard, 2026

if nargin < 2, outfile = ''; end
if ~exist('pop_loadset','file'), eeglab; close; end
bb = fileparts(fileparts(mfilename('fullpath')));
p = strsplit(path, pathsep);
old = p(contains(p,'BrainBeats','IgnoreCase',true) & ~startsWith(p,bb));
if ~isempty(old), rmpath(old{:}); end
addpath(bb, fullfile(bb,'functions'));

common = {'analysis','hep','clean_eeg',false,'hep_tf_freqs',[25 26],'hep_surrogates',0, ...
    'vis_cleaning',false,'vis_outputs',false,'save',false,'gong',false};
ppg = {'heart_signal','PPG','heart_channels',{'PPG'}};
files = cellstr(files);
rows = {};
for k = 1:numel(files)
    [fp, fn, ext] = fileparts(files{k});
    EEG = pop_loadset('filename',[fn ext],'filepath',fp);
    E = brainbeats_process(pop_select(EEG,'nochannel',{'PPG'}), common{:}, 'heart_signal','ECG','heart_channels',{'ECG'});
    runs = { 'off',              EEG,                                  {'ppg_transit','off'}
             'ECG (auto)',       EEG,                                  {}
             'EEG (no ECG)',     pop_select(EEG,'nochannel',{'ECG'}), {'ppg_transit','eeg'}
             'literature 250 ms', pop_select(EEG,'nochannel',{'ECG'}), {'ppg_transit',250}
             'ECG, PPG peaks',   EEG,                                  {'ppg_detect_mode','peaks'} };
    for r = 1:size(runs,1)
        try
            P = brainbeats_process(runs{r,2}, common{:}, ppg{:}, runs{r,3}{:});
        catch ME   % e.g. no clear cardiac field for the EEG estimate
            fprintf('%s %-20s failed: %s\n', fn, runs{r,1}, ME.message);
            rows(end+1,:) = {[fn ext], runs{r,1}, NaN, NaN, NaN, NaN, NaN, NaN, 0}; %#ok<AGROW>
            continue
        end
        pat = NaN; iqrP = NaN; ratio = NaN;
        if isfield(P.brainbeats.preprocessings,'ppg_transit')
            pt = P.brainbeats.preprocessings.ppg_transit;
            if isfield(pt,'median'), pat = pt.median; elseif isfield(pt,'pat'), pat = pt.pat; end
            if isfield(pt,'iqr'), iqrP = pt.iqr; end
            if isfield(pt,'ratio'), ratio = pt.ratio; end
        end
        [rr, off, pk] = similarity(E, P);
        rows(end+1,:) = {[fn ext], runs{r,1}, pat, iqrP, ratio, rr, off, pk, P.trials}; %#ok<AGROW>
        fprintf('%s %-20s PAT %5.0f ms  r = %.2f  offset %4.0f ms  peak %.2f\n', fn, runs{r,1}, pat, rr, off, pk);
    end
end
T = cell2table(rows, 'VariableNames', {'file' 'correction' 'pat' 'pat_iqr' 'ratio' 'r' 'offset' 'peak' 'epochs'});
disp(T)
if ~isempty(outfile), writetable(T, outfile); end
end

function [r, offset, peak] = similarity(E, P)
% PPG-locked HEP against the ECG-locked HEP (common channels)
[lab, iE, iP] = intersect({E.chanlocs.labels}, {P.chanlocs.labels}, 'stable');
a = mean(E.data(iE,:,:),3); b = mean(P.data(iP,:,:),3);
a = a - mean(a,2); b = b - mean(b,2);     % no cleaning: remove each channel's offset
t = E.times;
w = t >= 0 & t < 600;
r = mean(arrayfun(@(c) corr(a(c,w)', b(c,w)'), 1:numel(lab)));
g = @(x) std(x - mean(x,1), 0, 1);      % global field power (average reference)
q = t >= -100 & t <= 150; tq = t(q);
ga = g(a); gb = g(b);
[ma, ia] = max(ga(q)); [mb, ib] = max(gb(q));
offset = tq(ib) - tq(ia);
peak = mb / ma;
end
