function R = ica_compare(files, outdir, conds)
% ICA_COMPARE - Compare ICA settings of the HEP cleaning on real recordings:
% number of dipolar ICs (residual variance < 15%, standard BEM), ICLabel
% classes, ICA time, eye activity left, and the ROI HEP.
%
% Each recording x setting runs brainbeats_process (ECG HEP, 'clean_eeg',
% 'linenoise' 50, 100 surrogates, time-frequency limited to 8-10 Hz to save
% time) and saves one result file in outdir (done ones are skipped, so the
% run can be resumed).
%
% Usage:
%   R = ica_compare(files, outdir)            % all settings below
%   R = ica_compare(files, outdir, {'default' 'lp30_epochs'})
%
% Settings (ICA method: Picard 'P_' or Infomax 'I_'):
%   default       continuous 1-Hz copy, no low-pass, ASR 80 (the HEP default)
%   lp30          the same low-passed at 30 Hz
%   lp30_noasr    30 Hz, no ASR
%   lp30_100Hz    30 Hz, resampled to 100 Hz
%   epochs        broadband, epoched like the HEP data (bad epochs removed)
%   lp30_epochs_noasr  30 Hz, epoched, no ASR (close to the v1.6 default)
%
% The README numbers (2026-11-03; 6 recordings: the tutorial dataset and
% sub-033, 035, 045, 047, 053 of ds003838 'derived') come from an earlier
% version of this script, whose settings were the same but named after the
% defaults of that time.
%
% Runtime: 3-15 min per recording and setting (ICA, dipole fits).
%
% Copyright (C) - Cedric Cannard, 2026

bbroot = fileparts(fileparts(mfilename('fullpath')));
if ~exist('pop_loadset','file'), eeglab; close; end
p = strsplit(path, pathsep);
old = p(contains(p,'BrainBeats','IgnoreCase',true) & ~startsWith(p,bbroot));
if ~isempty(old), rmpath(old{:}); end
addpath(bbroot, fullfile(bbroot,'functions'));
if ~exist(outdir,'dir'), mkdir(outdir); end

S = struct('default', {{}}, 'lp30', {{'ica_lowpass',30}}, 'lp30_noasr', {{'ica_lowpass',30,'ica_asr',0}}, ...
    'lp30_100Hz', {{'ica_lowpass',30,'ica_srate',100}}, 'epochs', {{'ica_epochs',true}}, ...
    'lp30_epochs_noasr', {{'ica_lowpass',30,'ica_asr',0,'ica_epochs',true}});
if nargin < 3 || isempty(conds), conds = fieldnames(S)'; end
methods = {'P_' 1; 'I_' 2};
dip = fullfile(fileparts(which('dipfitdefs')),'standard_BEM');
R = [];
for k = 1:numel(files)
    [fp, fn, ext] = fileparts(files{k});
    for c = 1:numel(conds)
        for m = 1:2
            cname = [methods{m,1} conds{c}];
            out = fullfile(outdir, sprintf('icacmp_%s_%s.mat', fn, cname));
            if exist(out,'file'), r = load(out); R = [R r.r]; continue, end %#ok<AGROW>
            rng(1)
            EEG = pop_loadset('filename',[fn ext],'filepath',fp);
            EEG = pop_select(EEG,'nochannel',intersect({'PPG'},{EEG.chanlocs.labels}));
            t0 = tic;
            H = brainbeats_process(EEG,'analysis','hep','heart_signal','ECG','heart_channels',{'ECG'}, ...
                'clean_eeg',true,'linenoise',50,'icamethod',methods{m,2},'hep_surrogates',100,'hep_tf_freqs',[8 10], ...
                S.(conds{c}){:},'vis_cleaning',false,'vis_outputs',false,'save',false,'gong',false);
            r = struct('file',fn,'setting',cname,'total_s',toc(t0),'trials',H.trials);
            pre = H.brainbeats.preprocessings;
            r.ica = pre.ica; r.removed = pre.removed_eeg_components;
            L = pre.ica.iclabel; [~, cls] = max(L, [], 2);
            r.n_ic = size(L,1); r.n_brain = sum(cls == 1); r.n_eye = sum(cls == 3); r.n_muscle = sum(cls == 2);
            % dipolarity: standard BEM, electrodes of standard_1005
            try
                D = eeg_emptyset; D.chanlocs = H.chanlocs; D.nbchan = numel(H.chanlocs);
                D.data = zeros(D.nbchan, 10); D.pnts = 10; D.srate = H.srate; D.trials = 1;
                D.icawinv = double(pre.ica.winv); D.icaweights = pinv(D.icawinv); D.icasphere = eye(D.nbchan);
                D.icachansind = 1:D.nbchan;
                D = eeg_checkset(D);
                D = pop_dipfit_settings(D,'hdmfile',fullfile(dip,'standard_vol.mat'),'mrifile',fullfile(dip,'standard_mri.mat'), ...
                    'chanfile',fullfile(dip,'elec','standard_1005.elc'),'coordformat','MNI','coord_transform',[0 0 0 0 0 0 1 1 1], ...
                    'chansel',1:D.nbchan);
                D = pop_multifit(D, 1:size(D.icawinv,2), 'threshold', 100, 'dipplot','off');
                r.rv = [D.dipfit.model.rv]';
            catch ME
                warning('dipfit: %s', ME.message); r.rv = nan(r.n_ic,1);
            end
            r.n_dip15 = sum(r.rv < .15); r.n_brain_dip15 = sum(r.rv < .15 & cls == 1);
            % eye activity left: frontal-pole peak-to-peak of the kept epochs
            fpc = ismember({H.chanlocs.labels}, {'Fp1' 'Fp2' 'AF7' 'AF8' 'AF3' 'AF4' 'Fpz'});
            x = squeeze(mean(H.data(fpc,:,:),1));
            r.fp_p2p_median = median(max(x) - min(x));
            t = H.brainbeats.roi.hep_times;
            r.roi_hep_390_430 = mean(H.brainbeats.roi.hep(1, t >= 390 & t <= 430));
            save(out, 'r')
            R = [R r]; %#ok<AGROW>
            fprintf('%s %-22s ICA %4.0f s, rank %d | dipolar %2d (brain %2d) | brain %2d eye %d muscle %2d | removed %d | Fp p2p %.1f uV\n', ...
                fn, cname, r.ica.seconds, r.ica.rank, r.n_dip15, r.n_brain_dip15, r.n_brain, r.n_eye, r.n_muscle, ...
                numel(r.removed), r.fp_p2p_median);
        end
    end
end
