% BRAINBEATS_PROCESS - Process an EEGLAB dataset containing EEG and
% cardiovascular (ECG or PPG) signals recorded simultaneously.
%
% Usage:
%   [EEG, com] = brainbeats_process(EEG);                 % GUI (also: EEGLAB menu 'BrainBeats')
%   [EEG, com] = brainbeats_process(EEG, 'key', val, ...) % command line
%
% Main inputs:
%   'analysis'       - 'hep'       heartbeat-evoked potentials and oscillations
%                      'features'  EEG and HRV features
%                      'coherence' brain-heart coherence (experimental, not
%                                  for interpretation: see the README)
%   'heart_signal'   - 'ecg', 'ppg', 'rr' (beat times detected elsewhere), or
%                      'off' (EEG features only)
%   'heart_channels' - cell array of heart channel labels, e.g. {'ECG'}. With
%                      several channels, the one with the fewest abnormal
%                      heartbeats (then the best SQI) is used.
%   'beat_latencies' - ('rr' only) beat times in s from the start of the data.
%                      Peak detection, SQI and RR cleaning are skipped.
%   'clean_eeg'      - preprocess the EEG (default false): FIR filter (0.5-30 Hz
%                      for 'hep', 1-40 Hz otherwise),
%                      re-reference, remove and
%                      interpolate bad channels, remove bad epochs ('hep') or
%                      segments with ASR, then remove artifactual ICA components
%                      with ICLabel (or GEDAI instead, see 'clean_method').
%                      'hep' without it: the EEG must be high-pass filtered
%                      already (the HEP is not baseline-corrected), else an
%                      error.
%   'vis_cleaning'   - plot the preprocessing steps (default true)
%   'vis_outputs'    - plot the outputs (default true)
%   'save'           - save the output dataset next to the input file, and the
%                      features in a .mat file (default true)
%   'parpool'        - parallel computing for the EEG features and GEDAI
%                      (default false)
%   'gpu'            - GPU computing for the EEG features and cleaning
%                      (default false)
%   'gong'           - play a sound when done (default true)
%   Sampling rate: at most 1024 Hz, 512 Hz recommended (higher rates only cost
%   memory and time). Above 1024 Hz the command line stops with an error
%   (EEG = pop_resample(EEG, 512) first); above 512 Hz it warns. The GUI asks.
%
% Heart options:
%   'ppg_detect_mode' - 'valleys' (default, pulse-wave onsets) or 'peaks'
%                      (systolic peaks). An inverted PPG (light intensity, as
%                      most raw device outputs) is detected and flipped first
%                      ('ppg_polarity' to force it). PPG pulses are detected
%                      with a prominence relative to the local pulse amplitude,
%                      and their quality is scored on the filtered PPG (bad =
%                      beats of unacceptable quality, Li & Clifford, 2012).
%   'ecg_peakthresh', 'ecg_refperiod', 'ecg_searchback', 'ecg_bandpass',
%   'ecg_polarity', 'ecg_adaptive_pol', 'ppg_bandpass', 'ppg_height_method',
%   'ppg_polarity'   - see GET_RR
%   'keep_heart'     - keep the heart channel in the output (default false)
%   RR artifacts are corrected by CLEAN_RR (implausible beats removed, missing
%   beats inserted in gaps, ectopic beats replaced).
%   'ectopic_rule'   - 'adaptive' (default: Lipponen & Tarvainen, 2019,
%                      thresholds scaled to the recording's beat-to-beat
%                      variability), 'ratio' (fixed 20% rule of v1.6) or 'none'.
%                      A warning flags heart rhythms the cleaning cannot fix:
%                      frequent ectopic beats or atrial fibrillation
%                      (RR_ALTERNATION above 10%).
%   For ECG, RPEAK_QA reports beats detected on the wrong wave (e.g. T-wave).
%
% HEP options ('hep'):
%   'hep_window'       - epoch window [start end] in ms (default [-300 600]), or
%                        'adaptive': the end is set from this recording's heart
%                        rate (within-subject analyses; see RUN_HEP). Heartbeats
%                        followed by the next one before the end + 50 ms are
%                        rejected.
%   'ppg_transit'      - (PPG) pulse arrival time (PAT), to shift the PPG beats
%                        back to the heartbeats. 'auto' (default): measured from
%                        an ECG channel of the file (ESTIMATE_PAT); without an
%                        ECG, an error asks for it: a delay in ms (e.g. measured
%                        with an ECG on the same device and sensor site), 'eeg'
%                        (estimated from the cardiac field artifact of the EEG,
%                        ESTIMATE_PAT_EEG: check its figure; on 6 test subjects
%                        it was right on 3, wrong on 1 and not clear enough on
%                        2), the label of an ECG channel, or 'off'.
%   'hep_baseline_win' - window in ms of each heartbeat's exported baseline
%                        voltage (default [-150 -50], which ends before the QRS
%                        onset). The HEP is not baseline-corrected (Steinfath et
%                        al., 2026); the baselines are for a regression-based
%                        correction in the group-level model (Alday, 2019).
%                        ('hep_baseline' was removed in v1.7.)
%   'hep_roi'          - channel labels of the ROI of the time-frequency
%                        measures and plots (default frontocentral: F1-F4, Fz,
%                        FC1-FC4, FCz, C1-C4, Cz; those present are used)
%   'hep_level'        - where the HEP, HRSP and HRPC (and surrogate control)
%                        are computed: 'channels' (default: all channels, in
%                        EEG.brainbeats.hrsp and .surrogate, and the ROI of the
%                        plots, .roi), 'ics' (all independent components, in
%                        EEG.brainbeats.ics with their maps and ICLabel classes;
%                        no IC is selected) or 'both'. The HEP epochs of all
%                        channels are the output dataset in every case.
%   'hep_tf_freqs'     - frequency range in Hz (default [4 30]) of the HRSP
%                        (heartbeat-related spectral perturbations) and HRPC
%                        (heartbeat-related phase consistency: pairwise phase
%                        consistency across heartbeats); see COMPUTE_HEP_TF
%   'hep_surrogates'   - number of surrogate heartbeat trains to test the
%                        heartbeat locking of the HEP, HRSP and HRPC, e.g. 100
%                        (default 0 = none; 100 in the GUI). Results in
%                        EEG.brainbeats.surrogate (HEP, all channels) and
%                        EEG.brainbeats.roi or .ics (ROI/IC HEP, HRSP, HRPC).
%                        One surrogate beat at a random position in each
%                        inter-beat interval, so at random cardiac phases
%                        ('hep_surrogate_mode','rigid' was removed in v1.7).
%   'hep_surr_stat'    - points of the HEP, HRSP and HRPC shown as different
%                        from the surrogates (p < .05): 'tmax' (default: max
%                        statistic, family-wise over all points), 'tfce'
%                        (threshold-free cluster enhancement along time,
%                        family-wise; slower, computed only when chosen) or
%                        'fdr' (per point, FDR-corrected). tmax and FDR are
%                        always stored. The HEP uses at least 1000
%                        surrogates (HRSP and HRPC: 'hep_surrogates').
%   'hep_tf_power'     - HRSP plotted: 'induced' (default: the HEP, cardiac field
%                        artifact included, subtracted from each heartbeat
%                        before the power is taken; Cohen & Donner, 2013) or
%                        'total'. Both are stored (.hrsp, .hrsp_induced).
%   'hep_cfa'          - QRS map of the cardiac field artifact and whether it
%                        changed during the recording (default true; CFA_MAP):
%                        EEG.brainbeats.cfa, with notes on how to use it; with
%                        'save', also <file>_HEP_CFA.png and _HEP_CFA.txt.
%   Each heartbeat's baseline voltage ('hep_baseline_win') and band power
%   before the QRS are stored in EEG.brainbeats.prestim (and, with 'save', the
%   ROI values in <file>_HEP_beats.csv), for group-level covariates or
%   single-trial models.
%
% EEG preprocessing options (with 'clean_eeg'):
%   'clean_method'   - 'asr_ica' (default: bad epochs ('hep') or ASR, then ICA +
%                      ICLabel) or 'gedai' (GEDAI plugin, installed if missing;
%                      Ros et al., 2025): removes artifacts from the continuous
%                      data without ICA
%   'ref'            - 'average' (default), 'infinity' (REST), 'csd' (surface
%                      Laplacian, applied after artifact removal; channels must
%                      have 10-05 labels) or 'off'
%   'highpass' (default 0.5 Hz for 'hep', 1 Hz otherwise; ICA is fitted on a
%   1-Hz high-passed copy), 'lowpass' (Hz), 'filttype' ('noncausal' default,
%   or 'causal'), 'linenoise' (50 or 60 Hz),
%   'ica_asr', 'ica_lowpass' ('hep': the copy the ICA is fitted on is
%   continuous, cleaned with ASR at 'ica_asr' SD, default 80, 0 = off, and
%   low-passed at 'ica_lowpass' Hz, default 0 = none, line noise notched;
%   ICLabel runs on it), 'ica_epochs' (true: epoch that copy like the HEP
%   data), 'ica_srate' (resample it),
%   'flatline', 'corrThresh', 'maxBad', 'asr_cutoff', 'detectMethod'
%   (bad epochs: 'grubbs' default, 'median', 'mean'), 'icamethod' or
%   'ica_method' (1 = Picard, fast; 2 = Infomax, default; 3 = Infomax with
%   lrate 1e-5, slow but replicable). See CLEAN_EEG.
%
% Features options ('features'):
%   'hrv_features'   - any of {'time' 'frequency' 'nonlinear'} (default all), or 'off'
%   'hrv_spec'       - 'LombScargle' (default, ms^2), 'LombScargle_norm', 'welch', 'fft'
%   'hrv_norm', 'hrv_overlap', 'prsa_thresh' - see GET_HRV_FEATURES
%   'eeg'            - 'off' to turn off all EEG operations
%   'eeg_features'   - any of {'time' 'frequency' 'nonlinear'} (default all), or 'off'
%   'eeg_frange', 'eeg_wintype', 'eeg_winlen', 'eeg_winoverlap', 'eeg_norm',
%   'eeg_freqbounds' ('conventional' default, or 'individualized': alpha set
%   from the individual alpha peak), 'asy_norm' - see GET_EEG_FEATURES
%   Nonlinear HRV and EEG features: fuzzy entropy and fractal dimension come
%   from the ASCENT plugin (2.1 or later, installed if missing; NaN without
%   it). Fuzzy entropy is higher than with BrainBeats 1.6 and earlier (see
%   BB_FUZZY_FRACTAL): do not pool it with older results.
%
% Heart components:
%   'heart_removal'  - ('hep' with 'clean_eeg') 'none' (default and
%                      recommended: see the README), 'ica' (heart components
%                      found by ICLabel, 'conf_thresh'), 'ica_ecg' (components
%                      whose time course correlates with the ECG, |r| >=
%                      'ecg_corr', default .6), 'ica_caracas' (cardiac
%                      components found by CARACAS, Champetier et al., in the
%                      SASICA plugin; no ECG needed), 'ecg_regression' (ECG shifted
%                      by -20 to +20 ms regressed out of each channel, see
%                      REMOVE_HEART_REGRESSION) or 'ssp' (signal-space
%                      projection of the 'ssp_n' first topographies of the
%                      heartbeat-locked average, default 2, see
%                      REMOVE_HEART_SSP). The last three need an ECG.
%   'conf_thresh'    - minimum ICLabel heart probability to remove a component
%                      (default 0.75)
%
% Outputs:
%   EEG - processed dataset. EEG.brainbeats holds the parameters used
%         (.parameters), the preprocessing outputs (.preprocessings, e.g. NN
%         intervals, SQI, removed channels, beats, components, PPG transit),
%         and the features (.features), coherence (.coherence), or for HEP the
%         ROI and IC results (.roi / .ics), HRSP/HRPC of all channels (.hrsp),
%         surrogate control (.surrogate), QRS map (.cfa) and pre-R-peak
%         values of each heartbeat (.prestim).
%   com - command line reproducing the call (EEGLAB history)
%
% Examples:
%   EEG = brainbeats_process(EEG, 'analysis','hep', 'heart_signal','ecg', ...
%       'heart_channels',{'ECG'}, 'clean_eeg',true, 'hep_surrogates',100);
%   EEG = brainbeats_process(EEG, 'analysis','hep', 'heart_signal','ppg', ...
%       'heart_channels',{'PPG'}, 'clean_eeg',true, 'hep_level','both');
%
% See brainbeats_tutorial.m for more examples.
%
% Reference:
%   Cannard, C., Wahbeh, H., & Delorme, A. (2024). BrainBeats as an
%   Open-Source EEGLAB Plugin to Jointly Analyze EEG and Cardiovascular
%   Signals. Journal of Visualized Experiments, (206).
%
% Copyright (C) - Cedric Cannard, 2023

function [EEG, com] = brainbeats_process(EEG, varargin)

Features = [];
com = '';

% Basic checks on inputs
if nargin < 1
    help brainbeats_process; return;
end

% Add path to subfunctions
mainpath = fileparts(which('eegplugin_BrainBeats.m'));
addpath(fullfile(mainpath, 'functions'));

% GUI: parameters (a dataset can be loaded from the main window)
if nargin == 1
    [params, abort, EEG] = getparams_gui(EEG);
    if abort
        disp('Aborted.'); return
    end
end

% Basic checks on EEG data (failures: an error on the command line, a dialog
% in the GUI)
isCmd = nargin > 1;
if isempty(EEG) || isempty(EEG.data)
    bb_fail(isCmd, 'Empty EEG dataset: load a dataset with EEG and ECG/PPG channels first.');
    return
end
if isempty(EEG.chanlocs(1).labels)
    bb_fail(isCmd, 'No channel labels.');
    return
end

% Command line parameters
if nargin > 1
    params = getparams_cmd(varargin{:});
end

% Routine checks, defaults, plugin installs, parallel pool
[EEG, params, err] = run_checks(EEG,params,isCmd);
if err, return; end  % errordlg does not stop execution
params.version = eegplugin_BrainBeats;

% Separate the heart channels (CARDIO) from the EEG
if strcmpi(params.heart_signal, 'rr')
    CARDIO = [];  % no heart channel when using pre-detected beats
else
    CARDIO = pop_select(EEG,'channel',params.heart_channels);
    EEG = pop_select(EEG,'nochannel',params.heart_channels);
end

% ECG channel used to estimate the PPG pulse arrival time ('ppg_transit'):
% the one named, or with 'auto' an ECG channel of the file if there is one
ECGREF = [];
if isfield(params,'ppg_transit') && (ischar(params.ppg_transit) || isstring(params.ppg_transit))
    params.ppg_transit = char(params.ppg_transit);
    labels = {EEG.chanlocs.labels};
    if strcmpi(params.ppg_transit, 'auto')
        ecgLabel = labels(contains(lower(labels), {'ecg' 'ekg'}));
    elseif any(strcmpi(params.ppg_transit, {'eeg' 'off'}))
        ecgLabel = {};
    else
        ecgLabel = labels(strcmpi(labels, params.ppg_transit));
        if isempty(ecgLabel)
            error("ECG channel '%s' ('ppg_transit') not found.", params.ppg_transit)
        end
    end
    if ~isempty(ecgLabel)
        ECGREF = pop_select(EEG,'channel',ecgLabel(1));
        EEG = pop_select(EEG,'nochannel',ecgLabel);
    end
end

% Other auxiliary channels would bias EEG preprocessing and features
% (detection based on channel labels only)
idx = contains(lower({EEG.chanlocs.labels}), {'ecg' 'ekg' 'ppg' 'aux' 'gsr' 'eda' 'eog' 'emg'});
if any(idx)
    auxChan = strjoin({EEG.chanlocs(idx).labels});
    opts.Interpreter = 'tex';
    opts.Default = 'Yes';
    quest = sprintf("The following channels may not be EEG or ECG/PPG and may cause serious inaccuracies: \n\n %s \n\n They are about to be removed, please confirm", auxChan);
    if exist('batchStartupOptionUsed') > 0 && batchStartupOptionUsed %#ok<EXIST>
        answer = 'Yes';   % no one to answer a dialog under 'matlab -batch'
        fprintf('Removing channels that may not be EEG or ECG/PPG: %s\n', auxChan);
    else
        answer = questdlg(quest,'Potentially undesirable electrodes','Yes','No','Cancel',opts);
    end
    if strcmp(answer,'Cancel')
        disp('Aborted.'); return
    elseif strcmp(answer,'Yes')
        EEG = pop_select(EEG,'nochannel',{EEG.chanlocs(idx).labels});
    end
end

%%%%% Heartbeats: RR, SQI, NN (and HRV features) %%%%%
if ~strcmpi(params.heart_signal,'off')

    if strcmp(params.analysis, 'hep') || isfield(params,'hrv_features')

      if strcmpi(params.heart_signal, 'rr')
        % Pre-detected beats: convert seconds to sample indices. Peak
        % detection, SQI and clean_rr are skipped (the beats are taken as
        % clean; clean_rr can distort intervals that were already cleaned).
        fprintf('Using pre-detected beat latencies (%g beats provided).\n', length(params.beat_latencies));
        beat_sec = params.beat_latencies(:)';
        Rpeaks = round(beat_sec * EEG.srate) + 1;   % EEGLAB sample 1 is t = 0
        Rpeaks(Rpeaks < 1) = [];
        Rpeaks(Rpeaks > EEG.pnts) = [];
        RR = diff(Rpeaks) / EEG.srate;         % sec
        fprintf('RR intervals: mean=%.0f ms, std=%.0f ms, n=%d\n', mean(RR)*1000, std(RR)*1000, length(RR));

        % For HRV features: the provided beats are taken as clean (NN = RR)
        NN   = RR;
        NN_t = (Rpeaks(2:end) - 1) / EEG.srate;   % sec, time of the beat ending each interval
        HR   = 60 ./ NN;
        EEG.brainbeats.preprocessings.NN = NN;
        EEG.brainbeats.preprocessings.NN_times = NN_t;

      else  % ECG or PPG

        % Heart and EEG data must share the sampling rate for HEP
        if  strcmp(params.analysis, 'hep') && EEG.srate~=CARDIO.srate
            fprintf('Resampling cardiovascular data to match EEG sampling rate. \n')
            CARDIO = pop_resample(CARDIO,EEG.srate);
        end

        % Flag bad portions of the cardiovascular signal before peak detection
        % ------------------------------------------------------------------
        % DETECT_ECG_ARTIFACTS marks the stretches get_RR should never see
        % (muscle bursts, cable movement, electrode pull, non-finite samples).
        % Whatever is cut from CARDIO must be cut from EEG over the same time
        % range, or the R-peak markers inserted later land on the wrong EEG
        % samples. Optional, and it shortens the recording, so it is left
        % commented: uncomment to use.
        %
        % % 1) Flag. HFthresh 15 catches gross artifacts only (5 is stricter);
        % %    AmpThresh 6 adds slow movement, Inf leaves that mask off.
        % [~, bad_seg] = detect_ecg_artifacts(CARDIO, 'HFthresh',15, ...
        %     'AmpThresh',6, 'Plot',true);
        %
        % if ~isempty(bad_seg)
        %     % 2) The same time ranges in EEG samples. Identical to bad_seg in
        %     %    HEP-mode, where CARDIO was resampled to EEG.srate above.
        %     bad_sec = [bad_seg(:,1)-1, bad_seg(:,2)] / CARDIO.srate;
        %     seg_eeg = [round(bad_sec(:,1)*EEG.srate)+1, round(bad_sec(:,2)*EEG.srate)];
        %     seg_eeg(:,2) = min(seg_eeg(:,2), EEG.pnts);
        %
        %     % 3) Record which samples survive, so VIS_ARTIFACTS can put the
        %     %    cleaned data back on the original time axis to compare.
        %     EEG_before = EEG;
        %     EEG.etc.clean_sample_mask = true(1,EEG.pnts);
        %     for iSeg = 1:size(seg_eeg,1)
        %         EEG.etc.clean_sample_mask(seg_eeg(iSeg,1):seg_eeg(iSeg,2)) = false;
        %     end
        %
        %     % 4) Cut both. pop_select shifts event latencies and drops the
        %     %    events that fall inside the removed ranges.
        %     EEG    = pop_select(EEG, 'nopoint', seg_eeg);
        %     CARDIO = pop_select(CARDIO, 'nopoint', bad_seg);
        %
        %     % 5) Inspect: the removed stretches show as gaps against the
        %     %    original. Press 'd' in the figure for the difference view.
        %     vis_artifacts(EEG, EEG_before);
        % end

        % RR and NN intervals for each heart channel. Outputs are stored in
        % structures (one field per channel) because beat counts differ
        % across channels; the best channel is kept afterwards.
        signal = CARDIO.data;
        sqi = [];
        nElec = size(signal,1);
        for iElec = 1:nElec
            elec = sprintf('elec%g',iElec);
            fprintf('Detecting R peaks from cardiovascular time series %g (%s)... \n', iElec, CARDIO.chanlocs(iElec).labels)
            [RR.(elec), RR_t.(elec), Rpeaks.(elec), sig(iElec,:), sig_t(iElec,:)] = get_RR(signal(iElec,:), CARDIO.times, params);

            % Sample indices are shared with the EEG below
            factor = EEG.srate / CARDIO.srate;
            if factor ~= 1
                error('BrainBeats:input', 'Your EEG and cardiovascular data must have the same sampling rate.')
            end

            % Signal quality index (SQI; Vest et al., 2018)
            % ECG: agreement of two QRS detectors per window; bad = SQI < .9.
            % PPG: beat-template correlation (Li & Clifford, 2012) on the
            % filtered PPG used for detection (the raw signal's drift lowered
            % it: 33% of the clean sample PPG was < .9, 6% once filtered); bad
            % = beats Li's rule calls unacceptable ('Q'), since the direct
            % template comparison also penalizes normal long beats (sinus
            % arrhythmia), so a .9 cut flags clean recordings.
            SQIthresh = .9; % minimum SQI recommended
            maxThresh = 20; % maximum % of bad windows or abnormal beats recommended
            if strcmpi(params.heart_signal, 'ecg')
                sqi.(elec) = get_sqi_ecg(Rpeaks.(elec), signal(iElec,:), params.fs);
                badRatio = sum(sqi.(elec) < SQIthresh) / length(sqi.(elec)) * 100;
                badWhat = 'of the ECG windows have an SQI below .9';
            elseif strcmpi(params.heart_signal, 'ppg')
                [sqi.(elec),~,annot] = get_sqi_ppg(Rpeaks.(elec), sig(iElec,:)', params.fs);
                badRatio = sum(strcmp(annot,'Q')) / numel(annot) * 100;
                badWhat = 'of the PPG beats are of unacceptable quality (Li & Clifford, 2012)';
            else
                error("Heart channel should be either 'ECG' or 'PPG' ")
            end
            % per-electrode structs: beat counts differ across electrodes, so
            % SQI values cannot share one numeric array
            SQI_mu.(elec) = round(mean(sqi.(elec), 'omitnan'),2);
            SQI_badRatio.(elec) = round(badRatio,1);
            if SQI_mu.(elec) < .9
                warning("Mean signal quality index (SQI): %g. Minimum recommended SQI = .9. See Vest et al. (2018) for more detail. \n",  SQI_mu.(elec))
            else
                fprintf("Mean SQI is within recommendations (>0.9): %g\n",  SQI_mu.(elec));
            end
            msg = sprintf("%g%% %s, before correction of RR artifacts. Maximum recommended: %g%%.", SQI_badRatio.(elec), badWhat, maxThresh);
            if SQI_badRatio.(elec) > maxThresh
                warning('%s', msg)
                bb_warndlg(msg,'Signal quality warning')
            else
                fprintf("%s \n", msg)
            end

            % Correct RR artifacts (e.g. ectopic beats, missed or false
            % detections) to obtain the NN series
            disp("Correcting abnormal RR intervals...")
            rrOpts = {};
            if isfield(params,'ectopic_rule'), rrOpts = {'ectopic_rule', params.ectopic_rule}; end
            [nn, nn_t, nPeaks, idx_bad.(elec), idx_interp.(elec)] = clean_rr(RR_t.(elec), RR.(elec),  sig(iElec, Rpeaks.(elec)), Rpeaks.(elec), ...
                'interpolate_missing', true, 'ecg_signal', sig(iElec,:), 'sig_t', sig_t(iElec,:), 'fs', params.fs, rrOpts{:});
            % Abnormal beats = removed + synthetic (inserted) + gaps left unfilled
            nRemoved = numel(RR.(elec)) - (numel(nn) - sum(idx_interp.(elec)));
            nFlagged.(elec) = nRemoved + sum(idx_interp.(elec)) + sum(idx_bad.(elec));
            flaggedRatio.(elec) = nFlagged.(elec) / numel(RR.(elec)) * 100;
            Npeaks.(elec) = nPeaks;
            NN.(elec) = nn;  NN_t.(elec) = nn_t;
        end

        % Keep the heart channel with the fewest abnormal heartbeats (ties,
        % e.g. none flagged on any channel, go to the best mean SQI). Note
        % that get_RR already ignores the first heartbeat.
        fr = cell2mat(struct2cell(flaggedRatio));
        sq = cell2mat(struct2cell(SQI_mu));
        cand = find(fr == min(fr));
        [~,j] = max(sq(cand));
        best_elec = cand(j);
        elec = sprintf('elec%g',best_elec);
        if nElec > 1
            fprintf('Using heart channel %s (fewest abnormal heartbeats, then best SQI). \n', CARDIO.chanlocs(best_elec).labels);
        end
        EEG.brainbeats.preprocessings.heart_channel_used = CARDIO.chanlocs(best_elec).labels;
        flaggedRatio = flaggedRatio.(elec);
        nFlagged = nFlagged.(elec);
        idx_bad = idx_bad.(elec);
        idx_interp = idx_interp.(elec);
        sig_t = sig_t(best_elec,:);
        sig = sig(best_elec,:);
        RR = RR.(elec);
        RR_t = RR_t.(elec);
        Rpeaks = Rpeaks.(elec);
        Npeaks = Npeaks.(elec);
        NN_t = NN_t.(elec);
        NN = NN.(elec);
        HR = 60 ./ NN;      % heart rate from the clean NN series (not the raw RR)
        SQI_mu = SQI_mu.(elec);
        SQI_badRatio = SQI_badRatio.(elec);
        if flaggedRatio > 0
            fprintf('Portion of abnormal heartbeats corrected: %g/%g (%.2f%%). \n', nFlagged, length(RR), flaggedRatio);
        end
        if flaggedRatio > maxThresh
            warning("%g%% of the RR series on your best electrode was flagged as artifact and corrected. Maximum recommendation is 20%%. You may want to check for abnormal sections (e.g. electrode disconnections for long periods of time) in your cardiovascular signal and try BrainBeats again. ", round(flaggedRatio,2));
        end
        % Frequent ectopic beats or atrial fibrillation: RR cleaning corrects
        % only occasional ectopic beats and leaves a short-long alternation
        % (RR_ALTERNATION; tests/validate_ectopic_mitbih.m)
        altRatio = 100 * rr_alternation(NN);
        EEG.brainbeats.preprocessings.rr_alternation = altRatio;
        if altRatio > 10
            if strcmp(params.analysis,'hep')
                what = 'HEP epochs include ectopic beats';
            else
                what = 'HRV measures are unreliable';
            end
            msg = sprintf(['%.1f%% of the heartbeats alternate between short and long intervals after RR cleaning: ' ...
                'frequent ectopic beats, or an arrhythmia such as atrial fibrillation. RR cleaning corrects only ' ...
                'occasional ectopic beats, so %s. Check the RR series (vis_cleaning).'], altRatio, what);
            warning('BrainBeats:irregularRR', '%s', msg)
            bb_warndlg(msg, 'Irregular heart rhythm')
        end

        % Print average SQI
        fprintf("Overall signal quality index (SQI): %g. Note: SQI > 0.9 is considered good.\n", SQI_mu)

        % R-peak misdetection check (beats locked on the T or S wave), which
        % SQI and RR cleaning cannot see. Report only.
        if strcmpi(params.heart_signal, 'ecg')
            qa = rpeak_qa(sig, Rpeaks, params.fs);
            EEG.brainbeats.preprocessings.rpeak_qa = qa;
            fprintf('R-peak QA: %s (median beat-template r = %.2f; %.1f%% low-correlation beats; %.1f%% low-amplitude beats; RR lag-1 = %.2f). \n', ...
                qa.flag, qa.MedCorr, 100*qa.FracLowCorr, 100*qa.FracLowAmp, qa.RRlag1);
            if ~any(strcmp(qa.flag, {'ok' 'too-few-beats'}))
                warning(['R-peak QA flagged "%s": some beats may be detected on the wrong wave (e.g. T-wave). ' ...
                    'Check the R-peaks in the ECG plot (vis_cleaning) before interpreting HEPs.'], qa.flag)
            end
        end

        % Plot the filtered heart signal with its beats, and the RR and NN
        % series of the channel kept
        if params.vis_cleaning
            plot_NN(sig_t, sig, RR_t, RR, Rpeaks, NN_t, NN, Npeaks, params.heart_signal)
            pause(1)  % draw it now rather than after the EEG preprocessing
        end

        % Preprocessing outputs
        EEG.brainbeats.preprocessings.bad_heartbeats = idx_bad; % gaps clean_rr could not fill
        EEG.brainbeats.preprocessings.interpolated_heartbeats = idx_interp; % synthetic beats inserted in gaps
        EEG.brainbeats.preprocessings.abnormal_heartbeats_percent = flaggedRatio; % removed + inserted + unfilled
        if exist('sqi','var')
            EEG.brainbeats.preprocessings.heart_SQI_mean = SQI_mu;
            EEG.brainbeats.preprocessings.heart_SQI_badportion = SQI_badRatio;
        end
        EEG.brainbeats.preprocessings.NN = NN;
        EEG.brainbeats.preprocessings.NN_times = NN_t;

        % HEP is time-locked to the beats clean_rr kept (not the ones it
        % removed, nor the synthetic ones it inserted in gaps); all detected
        % beats are kept for RUN_HEP's next-QRS rule (a removed premature
        % beat still has its QRS)
        params.detected_beats = Rpeaks;
        Rpeaks = Npeaks(~idx_interp);

        % PPG: shift the pulse beats back to the heartbeats by the pulse
        % arrival time (PAT): from an ECG channel, a delay given by the user,
        % or, if asked ('eeg'), the cardiac field artifact of the EEG. No
        % silent estimate without an ECG: on 6 subjects of ds003838 (pulse
        % onsets, PAT 92-112 ms from the ECG), the EEG estimate was right on 3,
        % wrong on 1 (404 ms, peak/median GFP 2.05, just above the 2 needed)
        % and not clear enough on 2, and a fixed 250 ms was wrong on all 6.
        if strcmpi(params.heart_signal,'ppg') && strcmp(params.analysis,'hep')
            if ~isfield(params,'ppg_transit') || isempty(params.ppg_transit)
                params.ppg_transit = 'auto';
            end
            pat = NaN; patInfo = struct();
            if isnumeric(params.ppg_transit)
                pat = params.ppg_transit;
                patInfo = struct('median',pat,'method','user');
            elseif strcmpi(params.ppg_transit,'off')
                warning(['HEP time-locked to the PPG pulses, which follow the heartbeats by ~100 ms (onsets) to ~400 ms (systolic peaks) ' ...
                    '(pulse arrival time): ''ppg_transit'' is off.'])
            elseif ~isempty(ECGREF)
                if ECGREF.srate ~= EEG.srate
                    ECGREF = pop_resample(ECGREF, EEG.srate);
                end
                fprintf('Estimating the PPG pulse arrival time from ECG channel %s... \n', ECGREF.chanlocs(1).labels)
                ecgParams = params;
                ecgParams.heart_signal = 'ecg';
                [~, ~, rEcg] = get_RR(ECGREF.data(1,:), ECGREF.times, ecgParams);
                [pat, ~, patInfo] = estimate_pat(rEcg, Rpeaks, EEG.srate);
                patInfo.method = sprintf('estimated from ECG channel %s', ECGREF.chanlocs(1).labels);
                fprintf('Pulse arrival time: median %.0f ms (IQR %.0f ms, %g/%g beats paired). \n', ...
                    pat, patInfo.iqr, patInfo.nPaired, patInfo.nPPG)
            else
                % No ECG: the QRS field in the EEG precedes each pulse by the PAT
                [patEEG, patInfo] = estimate_pat_eeg(EEG.data, EEG.srate, Rpeaks);
                patInfo.method = 'estimated from the EEG cardiac field artifact';
                if params.vis_outputs && ~isempty(patInfo.gfp)
                    figure('Color','w','Name','PPG pulse arrival time from the EEG','NumberTitle','Off');
                    plot(patInfo.times, patInfo.gfp, 'k', 'LineWidth', 1.5); hold on
                    if ~isnan(patInfo.pat), xline(-patInfo.pat, 'r--', sprintf('%.0f ms', patInfo.pat)); end
                    xline(0, 'b:', 'PPG beat'); box on
                    xlabel('Latency from the PPG beats (ms)'); ylabel('GFP of the median EEG (\muV)')
                    title('The QRS field should be a single clear peak before the PPG beats: check it')
                    finish_figure(gcf)
                end
                if ~strcmpi(params.ppg_transit,'eeg')
                    sugg = 'no clear cardiac field peak';
                    if ~isnan(patEEG), sugg = sprintf('%.0f ms, peak/median GFP %.1f', patEEG, patInfo.ratio); end
                    bb_fail(isCmd, sprintf(['No ECG channel to measure the PPG pulse arrival time, which is needed to ' ...
                        'time-lock the HEP to the heartbeats. Set ''ppg_transit'' to the delay in ms (e.g. measured with ' ...
                        'an ECG on the same device and sensor site), to ''eeg'' to use the estimate from the cardiac field ' ...
                        'artifact of the EEG (here: %s; check its figure: it was wrong on 1 of 6 test subjects), or to ''off''.'], sugg));
                    return
                end
                if isnan(patEEG)
                    bb_fail(isCmd, sprintf(['The cardiac field artifact of the EEG is not clear enough to estimate the ' ...
                        'PPG pulse arrival time (peak/median GFP = %.1f, 2 needed): set ''ppg_transit'' to the delay in ms.'], patInfo.ratio));
                    return
                end
                pat = patEEG;
                warning(['Pulse arrival time from the EEG cardiac field artifact: %.0f ms (peak/median GFP %.1f, %g beats). ' ...
                    'Check the figure: this estimate was wrong on 1 of 6 test subjects.'], pat, patInfo.ratio, patInfo.nBeats)
            end
            if ~isnan(pat)
                Rpeaks = Rpeaks - round(pat/1000*EEG.srate);
                Rpeaks(Rpeaks < 1) = [];
                params.detected_beats = params.detected_beats - round(pat/1000*EEG.srate);
                params.detected_beats(params.detected_beats < 1) = [];
                EEG.brainbeats.preprocessings.ppg_transit = patInfo;
                fprintf('PPG beats shifted by -%.0f ms to estimate the heartbeat times for HEP (%s). \n', pat, patInfo.method)
            end
        end

      end  % ECG/PPG


      %%%%% HRV features %%%%%
      if ~isempty(params.hrv_features) && params.hrv_features ~= 0

          % The whole series is used, to allow ULF and VLF when long enough
          file_length = floor(EEG.xmax)-1;
          if file_length < 300
              warning('File length is less than 5 minutes! The minimum recommended is 300 s for estimating reliable HRV metrics.')
              bb_warndlg('File length is less than 5 minutes! The minimum recommended is 300 s for estimating reliable HRV metrics.')
          end
          params.file_length = file_length;

          % Extract HRV measures
          [features_hrv, params] = get_hrv_features(NN, NN_t, params);

          Features.HRV = features_hrv;
          Features.HRV.time.heart_rate = round(mean(HR,'omitnan'),1);

          if strcmp(params.analysis,'features') && ~params.eeg_features
              EEG.brainbeats.features = Features;
              disp('Done processing cardiovascular signals')
          end
      end
    end
end

% EEG preprocessing, step 0: filter, re-reference, bad channels
if params.clean_eeg
    params.clean_eeg_step = 0;
    params.orichanlocs = EEG.chanlocs;
    if strcmp(params.analysis,'hep') && isfield(params,'clean_method') && strcmpi(params.clean_method,'gedai')
        params.cfa_beats = Rpeaks;   % heart-locked amplitude before and after GEDAI
    end
    [EEG, params] = clean_eeg(EEG, params);

    % Preprocessing outputs
    EEG.brainbeats.preprocessings.removed_eeg_channels = params.removed_eeg_channels;
    if isfield(params,'cfa_gedai')
        EEG.brainbeats.preprocessings.cfa_gedai = params.cfa_gedai;
        params = rmfield(params, {'cfa_gedai' 'cfa_beats'});
    end

    % HEP: cardiac field artifact removed by regressing the ECG out of the EEG
    if strcmp(params.analysis,'hep') && isfield(params,'heart_removal') && strcmp(params.heart_removal,'ecg_regression')
        if strcmpi(params.heart_signal,'ecg') && ~isempty(CARDIO)
            [EEG, EEG.brainbeats.preprocessings.heart_regression] = remove_heart_regression(EEG, CARDIO, Rpeaks, params);
        else
            warning('ECG regression needs an ECG channel: the heart artifact is not removed.')
        end
    end
    % HEP: cardiac field artifact projected out (SSP of the heartbeat-locked average)
    if strcmp(params.analysis,'hep') && isfield(params,'heart_removal') && strcmp(params.heart_removal,'ssp')
        [EEG, EEG.brainbeats.preprocessings.heart_ssp] = remove_heart_ssp(EEG, Rpeaks, params);
        % the same projection on the ICA copy: an ICA fitted with the projected
        % topographies would add them back when its artifact ICs are removed
        if isfield(params,'ica_cont')
            U = EEG.brainbeats.preprocessings.heart_ssp.patterns;
            if size(U,1) == params.ica_cont.nbchan
                C = double(params.ica_cont.data);
                params.ica_cont.data = cast(C - U * (U' * C), 'like', params.ica_cont.data);
                clear C
            else
                warning('SSP: the ICA copy has other channels than the data; its projection is skipped.')
            end
        end
    end
end

%%%%% MODE 1: Heartbeat-evoked potentials (HEP) %%%%%
if strcmp(params.analysis,'hep')
    EEG = run_HEP(EEG, CARDIO, params, Rpeaks); % CARDIO: to add the heart channel back ('keep_heart')
    if isfield(params,'ica_cont'), params = rmfield(params,'ica_cont'); end   % ICA copy (data), used by RUN_HEP

%%%%% MODE 2: EEG features %%%%%
elseif isfield(params,'eeg_features') && ~isempty(params.eeg_features) && params.eeg_features~=0

    % EEG preprocessing, step 1: bad segments (ASR) and bad components (ICLabel)
    if params.clean_eeg
        [EEG, params] = clean_eeg(EEG, params);
        if isfield(params,'removed_eeg_segments')
            EEG.brainbeats.preprocessings.removed_eeg_segments = params.removed_eeg_segments;
        end
        if isfield(params,'removed_eeg_components')
            EEG.brainbeats.preprocessings.removed_eeg_components = params.removed_eeg_components;
        end
        if isfield(params,'ica_report')
            EEG.brainbeats.preprocessings.ica = params.ica_report;
            params = rmfield(params,'ica_report');
        end
        if strcmp(params.ref,'csd')
            EEG = apply_csd(EEG);   % surface Laplacian, after artifact removal
        end
    end

    params.chanlocs = EEG.chanlocs;
    [features_eeg, params] = get_eeg_features(EEG.data, params);
    Features.EEG = features_eeg;

%%%%% MODE 3: Brain-heart coherence %%%%%
elseif strcmpi(params.analysis,'coherence')
    disp("Computing brain-heart coupling...")
    warning(['Brain-heart coherence is experimental: it relates the EEG at 0-40 Hz to the NN series, whose ' ...
        'variations (HRV) are below 0.5 Hz and are high-passed at 0.5 Hz here. Do not interpret these values.'])

    if length(params.heart_channels)>1
        error("Sorry, this method only supports one hear channel at the moment. Please select only one of your heart channels and try again.")
    end

    % Coherence is computed with the NN series (not the raw ECG/PPG), which
    % is unevenly sampled: keep the EEG samples within the time range of the
    % NN series and interpolate it (cubic spline) at those samples, so both
    % share one time grid (beat times fall between samples since 2026-11-09)
    idx = EEG.times/1000 >= NN_t(1) & EEG.times/1000 <= NN_t(end);
    t_resamp = EEG.times(idx);   % ms
    CARDIO = eeg_emptyset();
    CARDIO.chanlocs.labels = 'HRV';
    CARDIO.data = interp1(NN_t(:), NN(:), t_resamp(:)/1000, 'spline')';
    CARDIO.times = t_resamp;
    CARDIO.srate = EEG.srate;
    CARDIO.pnts = size(CARDIO.data,2);
    CARDIO.xmax = CARDIO.times(end)/1000;
    CARDIO = eeg_checkset(CARDIO);
    CARDIO.data = bsxfun(@minus, CARDIO.data, trimmean(CARDIO.data,20,2)); % demean to remove offset
    if params.vis_cleaning
        pop_eegplot(CARDIO,1,1,1);
        finish_figure(gcf)
    end

    % Keep the time range covered by the NN series
    EEG.times = t_resamp;
    EEG.data = EEG.data(:,idx);
    EEG.pnts = size(EEG.data,2);
    EEG.xmax = EEG.times(end)/1000;
    EEG = eeg_checkset(EEG);

    % EEG preprocessing, step 1: bad segments (ASR) and bad components (ICLabel)
    if params.clean_eeg
        [EEG, params] = clean_eeg(EEG, params);

        % Remove the same segments from the NN series
        if isfield(params,'removed_eeg_segments')
            CARDIO = pop_select(CARDIO,'nopoint', params.removed_eeg_segments);
        end

        if isfield(params,'removed_eeg_segments')
            EEG.brainbeats.preprocessings.removed_eeg_segments = params.removed_eeg_segments;
        end
        if isfield(params,'removed_eeg_components')
            EEG.brainbeats.preprocessings.removed_eeg_components = params.removed_eeg_components;
        end
        if isfield(params,'ica_report')
            EEG.brainbeats.preprocessings.ica = params.ica_report;
            params = rmfield(params,'ica_report');
        end
        if strcmp(params.ref,'csd')
            EEG = apply_csd(EEG);   % surface Laplacian, after artifact removal
        end
    end

    params.chanlocs = EEG.chanlocs;   % after clean_eeg, which interpolates removed channels

    % Filter the NN series and rescale it to the EEG amplitude range
    CARDIO = pop_eegfiltnew(CARDIO,'locutoff',0.5,'hicutoff',40);
    disp("Rescaling heart signal to match EEG scale...")
    for iChan = 1:CARDIO.nbchan
        CARDIO.data(iChan,:) = rescale( CARDIO.data(iChan,:), quantile(mean(EEG.data,1),.01)*2, quantile(mean(EEG.data,1),.99)*2 );
    end

    % Append the NN series to the EEG as an extra channel
    if EEG.pnts == CARDIO.pnts
        EEG.data(end+1:end+CARDIO.nbchan,:) = CARDIO.data;
        EEG.nbchan = EEG.nbchan + CARDIO.nbchan;
        for iChan = 1:CARDIO.nbchan
            EEG.chanlocs(end+1).labels = params.heart_channels{iChan};
        end
        EEG = eeg_checkset(EEG);
    else
        error('Trying to merge EEG and Cardiovascular signals back together, but they have different lengths')
    end

    EEG.brainbeats.coherence = compute_brainheart_coherence(EEG,params);
    disp("Brain-heart coherence outputs can be found in: EEG.brainbeats.coherence")

end

% Store parameters in EEG structure for reporting (without the beats passed
% to RUN_HEP)
if isfield(params,'detected_beats'), params = rmfield(params,'detected_beats'); end
EEG.brainbeats.parameters = params;

% Store, save and plot features
if strcmp(params.analysis,'features')
    EEG.brainbeats.features = Features;
    disp("Features are stored in the EEG.brainbeats.features field")

    % Save next to the input file
    if params.save
        [name, folder] = output_file(EEG, '_features.mat');
        outputPath = fullfile(folder, name);
        fprintf("Exporting all features in a .mat file at this location: %s \n", outputPath);
        save(outputPath,'Features','params');
    end

    % Plot features
    if params.vis_outputs
        plot_features(Features,params)
    end
end

% Command line reproducing this call (EEGLAB history, 'eegh')
if strcmpi(params.heart_signal,'rr')
    heartArgs = sprintf('''heart_signal'',''rr'',''beat_latencies'',%s', mat2str(params.beat_latencies(:)',8));
elseif strcmpi(params.heart_signal,'off')
    heartArgs = '''heart_signal'',''off''';
else
    heartArgs = sprintf('''heart_signal'',''%s'',''heart_channels'',{%s}', params.heart_signal, ...
        strjoin(strcat('''', params.heart_channels, ''''), ' '));
end
common = sprintf(',''vis_cleaning'',%g,''vis_outputs'',%g,''save'',%g', params.vis_cleaning, params.vis_outputs, params.save);
heartOpts = argstr(params, {'ecg_peakthresh' 'ecg_refperiod' 'ecg_searchback' 'ecg_bandpass' 'ecg_polarity' ...
    'ecg_adaptive_pol' 'ppg_bandpass' 'ppg_height_method' 'ppg_polarity' 'ppg_detect_mode' 'ectopic_rule'});
cleanOpts = sprintf(',''clean_eeg'',%g', params.clean_eeg);
if params.clean_eeg
    cleanOpts = [cleanOpts argstr(params, {'clean_method' 'ref' 'highpass' 'lowpass' 'linenoise' ...
        'filttype' 'flatline' 'corrThresh' 'maxBad' 'detectMethod' 'asr_cutoff' 'icamethod' 'heart_removal'})];
    if strcmp(params.analysis,'hep')
        cleanOpts = [cleanOpts argstr(params, {'conf_thresh' 'ecg_corr' 'ssp_n' 'ica_asr' 'ica_lowpass' 'ica_epochs' 'ica_srate'})];
    end
elseif strcmp(params.analysis,'hep') && isfield(params,'hep_level') && any(strcmp(params.hep_level, {'ics' 'both'}))
    cleanOpts = [cleanOpts argstr(params, {'highpass' 'filttype'})];   % the ICA of the IC measures uses them
end
switch params.analysis
    case 'hep'
        hepArgs = argstr(params, {'hep_level' 'hep_roi' 'hep_window' ...
            'hep_baseline_win' 'hep_tf_freqs' 'hep_surrogates' 'hep_surr_stat' 'hep_tf_power' ...
            'hep_cfa' 'ppg_transit' 'keep_heart'});
        com = sprintf('EEG = brainbeats_process(EEG,''analysis'',''hep'',%s%s%s%s%s);', ...
            heartArgs, heartOpts, cleanOpts, hepArgs, common);
    case 'features'
        opt = {'time' 'frequency' 'nonlinear'};
        if params.hrv_features
            idx = logical([params.hrv_time params.hrv_frequency params.hrv_nonlinear]);
            hrv_features = ['{' sprintf(' ''%s'' ', opt{idx}) '}'];
        else
            hrv_features = '''off''';
        end
        if params.eeg_features
            idx = logical([params.eeg_time params.eeg_frequency params.eeg_nonlinear]);
            eeg_features = ['{' sprintf(' ''%s'' ', opt{idx}) '}'];
        else
            eeg_features = '''off''';
        end
        featOpts = argstr(params, {'hrv_spec' 'hrv_norm' 'hrv_overlap' 'prsa_thresh' 'eeg_frange' 'eeg_wintype' ...
            'eeg_winlen' 'eeg_winoverlap' 'eeg_freqbounds' 'eeg_norm' 'asy_norm'});
        com = sprintf('EEG = brainbeats_process(EEG,''analysis'',''features'',%s%s%s,''hrv_features'',%s,''eeg_features'',%s%s,''parpool'',%g,''gpu'',%g%s);', ...
            heartArgs, heartOpts, cleanOpts, hrv_features, eeg_features, featOpts, params.parpool, params.gpu, common);
    otherwise   % 'coherence'
        com = sprintf('EEG = brainbeats_process(EEG,''analysis'',''%s'',%s%s%s%s%s);', ...
            params.analysis, heartArgs, heartOpts, cleanOpts, argstr(params, {'conf_thresh' 'keep_heart'}), common);
end

% Downsampled in the GUI (RUN_CHECKS): part of the history
if isfield(params,'resampled')
    com = sprintf('EEG = pop_resample(EEG, %g); %s', params.resampled(2), com);
end

% Final message with ref to cite
fprintf('\n')
fprintf("Done! Thank you for using the BrainBeats toolbox! Please cite: \n");
fprintf("Cannard, C., Wahbeh, H., Delorme, A. BrainBeats as an Open-Source EEGLAB Plugin to Jointly Analyze EEG and Cardiovascular Signals. J. Vis. Exp. (2024). \n")
fprintf("https://www.jove.com/t/65829/brainbeats-as-an-open-source-eeglab-plugin-to-jointly-analyze-eeg \n")

if params.gong
    gong
end


function out = argstr(params, names)
% 'key',value pairs of the parameters present in params, as text for the
% command line (EEGLAB history)
out = '';
for i = 1:numel(names)
    if ~isfield(params, names{i}) || isempty(params.(names{i})), continue; end
    v = params.(names{i});
    if ischar(v) || isstring(v)
        txt = sprintf('''%s''', char(v));
    elseif iscellstr(v)
        txt = ['{' strjoin(strcat('''', v(:)', ''''), ' ') '}'];
    elseif (isnumeric(v) || islogical(v)) && isscalar(v)
        txt = sprintf('%g', v);
    elseif isnumeric(v)
        txt = mat2str(v);
    else
        continue
    end
    out = [out sprintf(',''%s'',%s', names{i}, txt)]; %#ok<AGROW>
end
