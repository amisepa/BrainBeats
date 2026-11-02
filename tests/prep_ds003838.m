function files = prep_ds003838(datapath, subjects, outpath)
% PREP_DS003838 - Resting-state EEG + ECG datasets from OpenNeuro ds003838,
% made like the sample dataset (see PREP_SAMPLEDATA) but without the
% injected artifacts: EEG and heart channels merged, 10-20 channel
% locations, downsampled to 250 Hz.
%
% Usage:
%   files = prep_ds003838(datapath, {'sub-033' 'sub-034'}, outpath)
%
% Inputs:
%   datapath - folder of ds003838 (sub-XXX/eeg/sub-XXX_task-rest_eeg.set and
%              sub-XXX/ecg/sub-XXX_task-rest_ecg.set)
%   subjects - subject folders, e.g. {'sub-033'}
%   outpath  - folder for the merged datasets (sub-XXX_rest.set)
%
% Output:
%   files    - full paths of the datasets written
%
% Copyright (C) - Cedric Cannard, 2026

if ~exist('pop_loadset','file'), eeglab; close; end
locFile = fullfile(fileparts(which('dipfitdefs.m')),'standard_BEM','elec','standard_1005.elc');
if ~exist(outpath,'dir'), mkdir(outpath); end
subjects = cellstr(subjects);
files = cell(size(subjects));
for k = 1:numel(subjects)
    s = subjects{k};
    EEG = pop_loadset('filename',sprintf('%s_task-rest_eeg.set',s),'filepath',fullfile(datapath,s,'eeg'));
    EEG = pop_chanedit(EEG,'lookup',locFile);
    CARDIO = pop_loadset('filename',sprintf('%s_task-rest_ecg.set',s),'filepath',fullfile(datapath,s,'ecg'));
    n = min(EEG.pnts, CARDIO.pnts);
    if EEG.pnts ~= CARDIO.pnts
        warning('%s: EEG and heart signals differ by %g samples: both cut to %g.', s, EEG.pnts - CARDIO.pnts, n)
    end
    EEG = pop_select(EEG,'point',[1 n]);
    CARDIO = pop_select(CARDIO,'point',[1 n]);
    assert(EEG.srate == CARDIO.srate, '%s: different sampling rates', s)
    EEG.data(end+1:end+CARDIO.nbchan,:) = CARDIO.data;
    EEG.nbchan = EEG.nbchan + CARDIO.nbchan;
    for iChan = 1:CARDIO.nbchan
        EEG.chanlocs(end+1).labels = CARDIO.chanlocs(iChan).labels;
    end
    EEG = eeg_checkset(EEG);
    EEG = pop_resample(EEG, 250);
    EEG.setname = sprintf('%s rest', s);
    files{k} = fullfile(outpath, sprintf('%s_rest.set', s));
    pop_saveset(EEG,'filename',sprintf('%s_rest.set', s),'filepath',outpath);
end
