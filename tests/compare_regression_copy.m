function compare_regression_copy(files)
% COMPARE_REGRESSION_COPY - With 'heart_removal','ecg_regression', the ECG is
% regressed out of the EEG but not out of the copy the HEP ICA is fitted on,
% so the ICA (and ICLabel) see a cardiac field the data no longer have. Does
% regressing the copy too change the result? Runs each recording twice
% (clean_eeg, Picard): as BrainBeats does, and with a variant of
% BRAINBEATS_PROCESS (written to tempdir) that also regresses the copy.
% Prints, for each: the ICs removed, the heart-locked amplitude left after
% the cleaning (-50 to 100 ms and the late window), and the difference of
% the ROI HEPs (RMS over 0-600 ms, relative to the HEP's RMS).
%
% Usage:
%   compare_regression_copy({'sample_data/sub-047_rest.set'})
%
% Results (2026-11-05; sub-047, dataset.set, ds003838 sub-062 and sub-033):
% the same number of ICs removed (1-3), the same ones on 1 of 4 recordings
% (the decompositions differ); heart-locked amplitude after cleaning within
% 0.01 uV; ROI HEP 0-600 ms differs by 0.2-9.3% RMS, at most 0.06 uV.
% Regressing the copy too is not needed.
%
% Copyright (C) - Cedric Cannard, 2026

if ~exist('pop_loadset','file'), eeglab; close; end
bb = fileparts(fileparts(mfilename('fullpath')));
p = strsplit(path, pathsep);
old = p(contains(p,'BrainBeats','IgnoreCase',true) & ~startsWith(p,bb));
if ~isempty(old), rmpath(old{:}); end
addpath(bb, fullfile(bb,'functions'));
% variant: the copy regressed too (same function, same ECG)
src = fileread(fullfile(bb,'brainbeats_process.m'));
key = '[EEG, EEG.brainbeats.preprocessings.heart_regression] = remove_heart_regression(EEG, CARDIO, Rpeaks, params);';
assert(numel(strfind(src, key)) == 1, 'regression call not found')
src = strrep(src, key, [key newline 'if isfield(params,''ica_cont''), params.ica_cont = remove_heart_regression(params.ica_cont, CARDIO, Rpeaks, params); end']);
src = strrep(src, 'function [EEG, com] = brainbeats_process(', 'function [EEG, com] = brainbeats_process_regcopy(');
d = fullfile(tempdir, 'brainbeats_regcopy'); if ~exist(d,'dir'), mkdir(d); end
fid = fopen(fullfile(d,'brainbeats_process_regcopy.m'), 'w'); fwrite(fid, src); fclose(fid);
addpath(d);
args = {'analysis','hep','heart_signal','ECG','heart_channels',{'ECG'},'clean_eeg',true,'icamethod',1, ...
    'heart_removal','ecg_regression','hep_surrogates',0,'vis_cleaning',0,'vis_outputs',0,'save',0,'gong',0};
files = cellstr(files);
for k = 1:numel(files)
    [fp, fn, ext] = fileparts(files{k});
    EEG = pop_loadset('filename',[fn ext],'filepath',fp);
    EEG = pop_select(EEG,'nochannel',intersect({EEG.chanlocs.labels}, {'PPG'}));
    A = brainbeats_process(EEG, args{:});
    B = brainbeats_process_regcopy(EEG, args{:});
    t = A.times >= 0 & A.times < 600;
    ra = A.brainbeats.roi.hep; rb = B.brainbeats.roi.hep;
    rel = sqrt(mean((ra(t) - rb(t)).^2)) / sqrt(mean(ra(t).^2));
    pa = A.brainbeats.preprocessings; pb = B.brainbeats.preprocessings;
    fprintf(['REGCOPY %s: ICs removed %d (copy not regressed) vs %d (regressed), %d in common; heart-locked ' ...
        'amplitude after cleaning %.2f vs %.2f uV (late %.2f vs %.2f); ROI HEP 0-600 ms differs by %.1f%% RMS ' ...
        '(max %.2f uV)\n'], fn, numel(pa.removed_eeg_components), numel(pb.removed_eeg_components), ...
        numel(intersect(pa.removed_eeg_components, pb.removed_eeg_components)), pa.cfa_after, pb.cfa_after, ...
        pa.cfa_late_after, pb.cfa_late_after, 100*rel, max(abs(ra(t) - rb(t))));
end
end
