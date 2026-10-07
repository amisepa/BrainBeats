function regen_267_reference()
% REG_267_REFERENCE - regenerate the full PPG HEP reference with save=1,
% exactly the hep_ppg test params of run_tutorial_headless.
bbroot = 'C:\Users\ccann\Documents\MATLAB\BrainBeats';
if ~exist('pop_loadset','file'), eeglab; close; end
p = strsplit(path, pathsep);
old = p(contains(p,'BrainBeats','IgnoreCase',true) & ~startsWith(p,bbroot));
if ~isempty(old), rmpath(old{:}); end
addpath(bbroot, fullfile(bbroot,'functions'));
fprintf('brainbeats_process resolved to: %s\n', which('brainbeats_process'));
EEG = pop_loadset('filename','dataset.set','filepath',fullfile(bbroot,'sample_data'));
args = {'analysis','hep','heart_signal','PPG','heart_channels',{'PPG'}, ...
    'clean_eeg',true,'linenoise',50,'ref','infinity','highpass',.5, ...
    'lowpass',20,'filttype','causal','detectMethod','median','icamethod',1, ...
    'vis_cleaning',0,'vis_outputs',0,'gong',0, ...
    'save',true,'out_path',fullfile(bbroot,'sample_data')};
HEP = brainbeats_process(EEG, args{:});
fprintf('output: %dx%dx%d, brainbeats keys: %s\n', size(HEP.data,1), ...
    size(HEP.data,2), size(HEP.data,3), strjoin(fieldnames(HEP.brainbeats),','));
if isfield(HEP.brainbeats,'hrsp')
    fprintf('hrsp: %dx%dx%d nBeats %d\n', size(HEP.brainbeats.hrsp.hrsp,1), ...
        size(HEP.brainbeats.hrsp.hrsp,2), size(HEP.brainbeats.hrsp.hrsp,3), ...
        HEP.brainbeats.hrsp.nBeats);
end
save(fullfile(bbroot,'sample_data','regen_reference_result.mat'),'HEP');
fprintf('REGEN DONE\n');
end