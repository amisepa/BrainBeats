function T = compare_heart_removal(files, outfile)
% COMPARE_HEART_REMOVAL - Heart artifact left in the HEP after each removal
% method, for several datasets.
%
% For each dataset (EEG + an 'ECG' channel, continuous), runs the HEP
% analysis with the EEG cleaning (Picard ICA) and 'heart_removal' = 'none',
% 'ica' and 'ecg_regression', and the 'rm_heart' analysis (ICA with the ECG
% channel, ICLabel heart threshold .65). Reports the heart-locked EEG
% amplitude (CFA_AMPLITUDE: RMS over channels of the heartbeat-locked
% average, -50 to 100 ms and 200 to 500 ms) before any removal and left
% after each method, and the number of HEP epochs.
%
% Usage:
%   T = compare_heart_removal({'sub-01.set' 'sub-02.set'})
%   T = compare_heart_removal(files, 'heart_removal.csv')   % also saves T
%
% Inputs:
%   files   - cell array of .set files (full paths)
%   outfile - optional .csv file for the table
%
% Output:
%   T - table, one row per dataset and method: file, method, qrs_before,
%       qrs_after, late_before, late_after (uV), qrs_reduction (%), epochs
%
% Notes:
% - 'before' is the cleaned EEG (filter, reference, bad channels) before
%   heart removal; 'after' is the EEG after all cleaning. The 200-500 ms
%   amplitude also holds any neural HEP, so it bounds the residual artifact
%   from above.
% - The methods can also remove part of a heartbeat-locked brain response:
%   check that separately (e.g. a synthetic response added after each
%   R-peak, as for REMOVE_HEART_REGRESSION).
%
% Copyright (C) - Cedric Cannard, 2026

if ~exist('pop_loadset','file'), eeglab; close; end
bb = fileparts(fileparts(mfilename('fullpath')));
p = strsplit(path, pathsep);
old = p(contains(p,'BrainBeats','IgnoreCase',true) & ~startsWith(p,bb));
if ~isempty(old), rmpath(old{:}); end
addpath(bb, fullfile(bb,'functions'));

files = cellstr(files);
common = {'heart_signal','ECG','heart_channels',{'ECG'},'clean_eeg',true,'icamethod',1, ...
    'vis_cleaning',false,'vis_outputs',false,'save',false,'gong',false};
rows = {};
for k = 1:numel(files)
    [fp, fn, ext] = fileparts(files{k});
    EEG = pop_loadset('filename',[fn ext],'filepath',fp);
    labels = {EEG.chanlocs.labels};
    other = labels(contains(lower(labels), {'ppg' 'pleth' 'eog' 'emg' 'gsr' 'eda'}));
    if ~isempty(other), EEG = pop_select(EEG,'nochannel',other); end

    for m = {'none' 'ica' 'ecg_regression'}
        H = brainbeats_process(EEG, 'analysis','hep', common{:}, 'heart_removal',m{1});
        pre = H.brainbeats.preprocessings;
        qb = pre.cfa_before; lb = pre.cfa_late_before;
        if strcmp(m{1},'ecg_regression')   % amplitudes before the regression
            qb = pre.heart_regression.cfa_before; lb = pre.heart_regression.cfa_late_before;
        end
        rows(end+1,:) = {[fn ext], m{1}, qb, pre.cfa_after, lb, pre.cfa_late_after, H.trials}; %#ok<AGROW>
    end
    H = brainbeats_process(EEG, 'analysis','rm_heart', common{:}, 'conf_thresh',.65);
    pre = H.brainbeats.preprocessings;
    if isfield(pre,'cfa_before')
        rows(end+1,:) = {[fn ext], 'rm_heart', pre.cfa_before, pre.cfa_after, pre.cfa_late_before, pre.cfa_late_after, NaN}; %#ok<AGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'file' 'method' 'qrs_before' 'qrs_after' 'late_before' 'late_after' 'epochs'});
T.qrs_reduction = 100 * (1 - T.qrs_after ./ T.qrs_before);
disp(T)
if nargin > 1 && ~isempty(outfile)
    writetable(T, outfile);
end
