function files = cache_coupling_data(inFiles, outPath, linenoise)
% CACHE_COUPLING_DATA - Clean EEG + heart datasets once and save what the
% brain-heart coupling analysis needs, so that its options can be compared
% without cleaning again: the cleaned EEG on its original time line (NaN
% where ASR removed data), the ECG R-peaks, NN intervals, the filtered ECG,
% and the PPG pulse onsets and filtered PPG. Same steps as BRAINBEATS_PROCESS
% ('coupling', 'clean_eeg' true, Picard ICA).
%
% Usage:
%   files = cache_coupling_data(inFiles, outPath)
%   files = cache_coupling_data(inFiles, outPath, linenoise)
%
% Inputs:
%   inFiles   - cell array of .set files with EEG, 'ECG' and 'PPG' channels
%               (e.g. made by PREP_DS003838)
%   outPath   - folder for the caches (<name>_coupling_cache.mat)
%   linenoise - line frequency (default 50)
%
% Recordings already cached are skipped. Requires EEGLAB on the path.
%
% Copyright (C) - Cedric Cannard, 2026

if nargin < 3, linenoise = 50; end
if ~exist(outPath, 'dir'), mkdir(outPath); end
inFiles = cellstr(inFiles);
files = cell(size(inFiles));
for k = 1:numel(inFiles)
    [fp, nm] = fileparts(inFiles{k});
    files{k} = fullfile(outPath, [nm '_coupling_cache.mat']);
    if exist(files{k}, 'file'), fprintf('%s: cached already\n', nm); continue; end
    try
        t0 = tic;
        EEG0 = pop_loadset('filename', [nm '.set'], 'filepath', fp);
        labels = {EEG0.chanlocs.labels};
        iE = find(strcmpi(labels,'ECG'), 1); iP = find(strcmpi(labels,'PPG'), 1);
        fs = EEG0.srate;
        ecg = double(EEG0.data(iE,:)); ppg = [];
        if ~isempty(iP), ppg = double(EEG0.data(iP,:)); end

        % heartbeats (as BRAINBEATS_PROCESS)
        hp = struct('fs', fs, 'heart_signal', 'ecg');
        [RR, RR_t, Rp, sig, sig_t] = get_RR(ecg, EEG0.times, hp);
        [NN, NN_t, Npeaks, ~, idx_interp] = clean_rr(RR_t, RR, sig(Rp), Rp, ...
            'interpolate_missing', true, 'ecg_signal', sig, 'sig_t', sig_t, 'fs', fs);
        beats = Npeaks(~idx_interp);
        pp = struct('onsets', [], 'sig', []);
        if ~isempty(ppg)
            [~, ~, pOn, psig] = get_RR(ppg, EEG0.times, struct('fs', fs, 'heart_signal', 'ppg'));
            pp = struct('onsets', pOn(:)', 'sig', single(psig(:)'));
        end

        % EEG cleaning, steps 0 and 1
        EEG = pop_select(EEG0, 'nochannel', iP);       % run_checks wants the heart channel
        params = getparams_cmd('analysis','coupling','heart_signal','ecg','heart_channels',{'ECG'}, ...
            'clean_eeg',true,'icamethod',1,'linenoise',linenoise,'vis_cleaning',false, ...
            'vis_outputs',false,'save',false,'gong',false);
        [EEG, params, err] = run_checks(EEG, params, true);
        if err, error('run_checks failed'); end
        EEG = pop_select(EEG, 'nochannel', {'ECG'});
        nPts0 = EEG.pnts;
        params.clean_eeg_step = 0; params.orichanlocs = EEG.chanlocs;
        [EEG, params] = clean_eeg(EEG, params);
        [EEG, params] = clean_eeg(EEG, params);
        keep = true(1, nPts0);
        if isfield(params,'removed_eeg_segments')
            for s = 1:size(params.removed_eeg_segments,1)
                keep(params.removed_eeg_segments(s,1):params.removed_eeg_segments(s,2)) = false;
            end
        end
        assert(sum(keep) == EEG.pnts, 'cleaned EEG does not match the samples kept')
        X = nan(EEG.nbchan, nPts0, 'single');
        X(:, keep) = EEG.data;
        nIC = 0;
        if isfield(params,'removed_eeg_components'), nIC = numel(params.removed_eeg_components); end
        C = struct('name', nm, 'fs', fs, 'chanlocs', EEG.chanlocs, 'X', X, 'beats', beats(:)', ...
            'NN', NN(:)', 'NN_t', NN_t(:)', 'ecg_sig', single(sig(:)'), 'ecg', single(ecg), ...
            'ppg', single(ppg), 'ppg_onsets', pp.onsets, 'ppg_sig', pp.sig, ...
            'removed_components', nIC, 'kept_share', mean(keep));
        save(files{k}, '-struct', 'C', '-v7.3');
        fprintf('%s: cached in %.0f s (%.1f%% kept, %d ICs removed)\n', nm, toc(t0), 100*mean(keep), C.removed_components);
    catch e
        fprintf('%s: FAILED: %s\n', nm, e.message);
        files{k} = '';
    end
    close all
end
