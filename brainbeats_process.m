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
%                      'rm_heart'  remove heart components from the EEG (ICA + ICLabel)
%                      'coherence' brain-heart coherence (beta)
%   'heart_signal'   - 'ecg', 'ppg', 'rr' (beat times detected elsewhere), or
%                      'off' (EEG features only)
%   'heart_channels' - cell array of heart channel labels, e.g. {'ECG'}. With
%                      several channels, the one with the fewest abnormal
%                      heartbeats (then the best SQI) is used.
%   'beat_latencies' - ('rr' only) beat times in s from the start of the data.
%                      Peak detection, SQI and RR cleaning are skipped.
%   'clean_eeg'      - preprocess the EEG (default false): FIR filter (0.5-30 Hz
%                      for 'hep', 1-30 Hz otherwise), re-reference, remove and
%                      interpolate bad channels, remove bad epochs ('hep') or
%                      segments with ASR, then remove artifactual ICA components
%                      with ICLabel (or GEDAI instead, see 'clean_method').
%   'vis_cleaning'   - plot the preprocessing steps (default true)
%   'vis_outputs'    - plot the outputs (default true)
%   'save'           - save the output dataset next to the input file, and the
%                      features in a .mat file (default true)
%
% Heart options:
%   'ppg_detect_mode' - 'valleys' (default, pulse-wave onsets) or 'peaks'. PPG
%                      pulses are detected with a prominence relative to the
%                      local pulse amplitude, and their quality is scored on the
%                      filtered PPG (bad = beats of unacceptable quality, Li &
%                      Clifford, 2012).
%   'ecg_peakthresh', 'ecg_refperiod', 'ecg_bandpass', 'ecg_polarity',
%   'ecg_adaptive_pol', 'ppg_bandpass', 'ppg_height_method' - see GET_RR
%   'keep_heart'     - keep the heart channel in the output (default false)
%   RR artifacts are corrected by CLEAN_RR (implausible beats removed, missing
%   beats inserted in gaps). For ECG, RPEAK_QA reports beats detected on the
%   wrong wave (e.g. T-wave).
%
% HEP options ('hep'):
%   'hep_window'       - epoch window [start end] in ms (default [-300 600]), or
%                        'adaptive': the end is set from this recording's heart
%                        rate (within-subject analyses; see RUN_HEP). Heartbeats
%                        followed by the next one before the end + 50 ms are
%                        rejected.
%   'ppg_transit'      - (PPG) pulse arrival time (PAT), to shift the PPG beats
%                        back to the heartbeats. 'auto' (default): measured from
%                        an ECG channel of the file if there is one (ESTIMATE_PAT),
%                        else from the cardiac field artifact of the EEG
%                        (ESTIMATE_PAT_EEG), else the literature value (250 ms to
%                        the pulse onset, 350 ms to the peak, finger PPG at rest).
%                        Or: the label of an ECG channel, 'eeg', a delay in ms,
%                        or 'off'.
%   'hep_baseline'     - 'none' (default; recommended by Steinfath et al., 2026:
%                        report the HEP without baseline correction and test the
%                        baseline window) or 'regression': regression-based
%                        baseline correction (Alday, 2019; see BASELINE_REGRESSION).
%                        The corrected epochs are stored in the output EEG.data.
%   'hep_baseline_win' - baseline window in ms (default [-150 -50], which ends
%                        before the QRS onset)
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
%   'hep_surrogate_mode' - 'shuffle' (default: shuffled inter-beat intervals, so
%                        surrogate beats fall at random cardiac phases) or 'rigid'
%                        (whole train shifted by up to +/-500 ms, for contrasts
%                        between conditions); see COMPUTE_HEP_TF
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
%   'flatline', 'corrThresh', 'maxBad', 'asr_cutoff', 'asr_mem', 'detectMethod'
%   (bad epochs: 'grubbs' default, 'median', 'mean'), 'icamethod' (1 = Picard,
%   fast; 2 = Infomax, default; 3 = Infomax with lrate 1e-5, slow but replicable).
%   See CLEAN_EEG.
%
% Features options ('features'):
%   'hrv_features'   - any of {'time' 'frequency' 'nonlinear'} (default all), or 'off'
%   'hrv_spec'       - 'LombScargle_norm' (default), 'LombScargle', 'welch', 'fft'
%   'hrv_norm', 'hrv_overlap' - see GET_HRV_FEATURES
%   'eeg'            - 'off' to turn off all EEG operations
%   'eeg_features'   - any of {'time' 'frequency' 'nonlinear'} (default all), or 'off'
%   'eeg_frange', 'eeg_wintype', 'eeg_winlen', 'eeg_winoverlap', 'eeg_norm',
%   'asy_norm'       - see GET_EEG_FEATURES
%   'parpool'        - use parallel computing (default false)
%   'gpu'            - use GPU computing (default false)
%
% Heart components:
%   'heart_removal'  - ('hep' with 'clean_eeg') 'ica' (default): the heart
%                      components found by ICLabel are removed with the other
%                      artifact components; 'ecg_regression': the ECG, shifted
%                      by -20 to +20 ms, is regressed out of each EEG channel
%                      (ECG only; see REMOVE_HEART_REGRESSION); or 'none'
%   'conf_thresh'    - minimum ICLabel heart probability to remove a component
%                      (default 0.75 for 'hep', 0.9 for 'rm_heart')
% rm_heart ('analysis','rm_heart', command line only; in the GUI, heart
% components are removed as part of the HEP analysis):
%   The heart-locked EEG amplitude (cardiac field artifact) is reported before
%   and after removal (see REMOVE_HEARTCOMP).
%
% Outputs:
%   EEG - processed dataset. EEG.brainbeats holds the parameters used
%         (.parameters), the preprocessing outputs (.preprocessings, e.g. NN
%         intervals, SQI, removed channels, beats, components, PPG transit),
%         and the features (.features), coherence (.coherence), or for HEP the
%         ROI and IC results (.roi / .ics), HRSP/HRPC of all channels (.hrsp)
%         and surrogate control (.surrogate).
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

% Basic checks on EEG data
if isempty(EEG) || isempty(EEG.data)
    errordlg('Empty EEG dataset: load a dataset with EEG and ECG/PPG channels first.');
    return
end
if isempty(EEG.chanlocs(1).labels)
    errordlg('No channel labels.');
    return
end

% Command line parameters
if nargin > 1
    params = getparams_cmd(varargin{:});
end

% Routine checks, defaults, plugin installs, parallel pool
[EEG, params, err] = run_checks(EEG,params);
if err, return; end  % errordlg does not stop execution

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
    quest = sprintf("The following channels may not be EEG or ECG/PPG and may cause serious innacuracies: \n\n %s \n\n They are about to be removed, please confirm", auxChan);
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
                errordlg("Your EEG and cardiovascular data must have the same sampling rate.")
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
                warndlg(msg,'Signal quality warning 1')
            else
                fprintf("%s \n", msg)
            end

            % Correct RR artifacts (e.g. ectopic beats, missed or false
            % detections) to obtain the NN series
            disp("Correcting abnormal RR intervals...")
            [nn, nn_t, nPeaks, idx_bad.(elec), idx_interp.(elec)] = clean_rr(RR_t.(elec), RR.(elec),  sig(iElec, Rpeaks.(elec)), Rpeaks.(elec), ...
                'interpolate_missing', true, 'ecg_signal', sig(iElec,:), 'sig_t', sig_t(iElec,:), 'fs', params.fs);
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
        % removed, nor the synthetic ones it inserted in gaps)
        Rpeaks = Npeaks(~idx_interp);

        % PPG: shift the pulse beats back to the heartbeats by the pulse
        % arrival time (PAT): from an ECG channel when there is one, else
        % from the cardiac field artifact of the EEG, else a literature value
        if strcmpi(params.heart_signal,'ppg') && strcmp(params.analysis,'hep')
            if ~isfield(params,'ppg_transit') || isempty(params.ppg_transit)
                params.ppg_transit = 'auto';
            end
            pat = NaN; patInfo = struct();
            if isnumeric(params.ppg_transit)
                pat = params.ppg_transit;
                patInfo = struct('median',pat,'method','user');
            elseif strcmpi(params.ppg_transit,'off')
                warning(['HEP time-locked to the PPG pulses, which follow the heartbeats by ~200-450 ms ' ...
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
                fprintf('Estimating the PPG pulse arrival time from the cardiac field artifact of the EEG... \n')
                [pat, patInfo] = estimate_pat_eeg(EEG.data, EEG.srate, Rpeaks);
                patInfo.method = 'estimated from the EEG cardiac field artifact';
                if isnan(pat)
                    % Finger PPG at rest, R-peak to pulse foot ~250 ms (Mukkamala
                    % et al., 2015: 180-260 ms; Kortekaas et al., 2018; Block
                    % et al., 2020), to the systolic peak ~350 ms (Charlier et
                    % al., 2026). Between-subject SD ~30 ms, more with device
                    % delays: the sample data give 424 ms (ECG).
                    if isfield(params,'ppg_detect_mode') && strcmpi(params.ppg_detect_mode,'peaks')
                        pat = 350;
                    else
                        pat = 250;
                    end
                    patInfo.pat = pat;
                    patInfo.method = 'literature default (no ECG, cardiac field artifact not clear in the EEG)';
                    warning(['No ECG channel, and the cardiac field artifact of the EEG is not clear enough ' ...
                        '(peak/median GFP = %.1f, 2 needed): using the literature default PAT of %g ms. ' ...
                        'Set ''ppg_transit'' to the delay in ms if you know it.'], patInfo.ratio, pat)
                else
                    fprintf('Pulse arrival time: %.0f ms (EEG cardiac field peak/median = %.1f, %g beats). \n', ...
                        pat, patInfo.ratio, patInfo.nBeats)
                end
            end
            if ~isnan(pat)
                Rpeaks = Rpeaks - round(pat/1000*EEG.srate);
                Rpeaks(Rpeaks < 1) = [];
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
              warndlg('File length is less than 5 minutes! The minimum recommended is 300 s for estimating reliable HRV metrics.')
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
    [EEG, params] = clean_eeg(EEG, params);

    % Preprocessing outputs
    EEG.brainbeats.preprocessings.removed_eeg_channels = params.removed_eeg_channels;

    % HEP: cardiac field artifact removed by regressing the ECG out of the EEG
    if strcmp(params.analysis,'hep') && isfield(params,'heart_removal') && strcmp(params.heart_removal,'ecg_regression')
        if strcmpi(params.heart_signal,'ecg') && ~isempty(CARDIO)
            [EEG, EEG.brainbeats.preprocessings.heart_regression] = remove_heart_regression(EEG, CARDIO, Rpeaks, params);
        else
            warning('ECG regression needs an ECG channel: the heart artifact is not removed.')
        end
    end
end

%%%%% MODE 1: Heartbeat-evoked potentials (HEP) %%%%%
if strcmp(params.analysis,'hep')
    EEG = run_HEP(EEG, CARDIO, params, Rpeaks); % CARDIO: to add the heart channel back ('keep_heart')

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

    if length(params.heart_channels)>1
        error("Sorry, this method only supports one hear channel at the moment. Please select only one of your heart channels and try again.")
    end

    % Coherence is computed with the NN series (not the raw ECG/PPG), which
    % is unevenly sampled: interpolate it at the EEG sampling rate
    [NN_resamp, t_resamp] = resample_NN(NN_t,NN,EEG.srate,'cub');
    t_resamp = t_resamp*1000; % ms
    CARDIO = eeg_emptyset();
    CARDIO.chanlocs.labels = 'HRV';
    CARDIO.data = NN_resamp;
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
    idx = EEG.times >= t_resamp(1) & EEG.times <= t_resamp(end);
    EEG.times = EEG.times(idx);
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

% Store parameters in EEG structure for reporting
EEG.brainbeats.parameters = params;

% Store, save and plot features
if strcmp(params.analysis,'features')
    EEG.brainbeats.features = Features;
    disp("Features are stored in the EEG.brainbeats.features field")

    % Save next to the input file
    if params.save
        outputPath = fullfile(EEG.filepath, sprintf('%s_features.mat', EEG.filename(1:end-4)));
        fprintf("Exporting all features in a .mat file at this location: %s \n", outputPath);
        save(outputPath,'Features');
    end

    % Plot features
    if params.vis_outputs
        plot_features(Features,params)
    end
end

%%%%% MODE 4: Remove heart (cardiac field artifact) components from the EEG %%%%%
if strcmp(params.analysis,'rm_heart')

    % Preprocess EEG
    if params.clean_eeg

        % Step 1: bad segments (ASR) and non-heart artifact components
        % (step 0 was done above)
        params.clean_eeg_step = 1;
        [EEG, params] = clean_eeg(EEG, params);

        % Filter the heart signal (default 1-20 Hz)
        if ~isfield(params,'highpass_ecg')
            params.highpass_ecg = 1;
        end
        if ~isfield(params,'lowpass_ecg')
            params.lowpass_ecg = 20;
        end
        CARDIO = pop_eegfiltnew(CARDIO,'locutoff',params.highpass_ecg,'hicutoff',params.lowpass_ecg);

        % Remove the same segments from the heart signal
        CARDIO = pop_select(CARDIO,'nopoint', params.removed_eeg_segments);

        EEG.brainbeats.preprocessings.removed_eeg_segments = params.removed_eeg_segments;
        EEG.brainbeats.preprocessings.removed_eeg_components = params.removed_eeg_components;
    end

    % Add the heart channel back: ICA uses it to isolate the heart components
    EEG.data(end+1:end+CARDIO.nbchan,:) = CARDIO.data;
    EEG.nbchan = EEG.nbchan + CARDIO.nbchan;
    for iChan = 1:CARDIO.nbchan
        EEG.chanlocs(end+1).labels = params.heart_channels{iChan};
    end
    EEG = eeg_checkset(EEG);

    EEG = remove_heartcomp(EEG, params);

    % Surface Laplacian after the heart components are removed (the ICA
    % weights kept in EEG describe the average-referenced data)
    if params.clean_eeg && strcmp(params.ref,'csd')
        EEG = apply_csd(EEG, params.heart_channels);
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
heartOpts = argstr(params, {'ecg_peakthresh' 'ecg_refperiod' 'ecg_searchback' 'ppg_detect_mode'});
cleanOpts = sprintf(',''clean_eeg'',%g', params.clean_eeg);
if params.clean_eeg
    cleanOpts = [cleanOpts argstr(params, {'clean_method' 'ref' 'highpass' 'lowpass' 'linenoise' ...
        'filttype' 'corrThresh' 'detectMethod' 'asr_cutoff' 'icamethod' 'heart_removal'})];
    if strcmp(params.analysis,'hep'), cleanOpts = [cleanOpts argstr(params, {'conf_thresh'})]; end
end
switch params.analysis
    case 'hep'
        hepArgs = argstr(params, {'hep_level' 'hep_roi' 'hep_window' 'hep_baseline' ...
            'hep_baseline_win' 'hep_tf_freqs' 'hep_surrogates' 'hep_surrogate_mode' 'ppg_transit' 'keep_heart'});
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
        featOpts = argstr(params, {'hrv_spec' 'hrv_norm' 'hrv_overlap' 'eeg_frange' 'eeg_wintype' ...
            'eeg_winlen' 'eeg_winoverlap' 'eeg_freqbounds' 'eeg_norm' 'asy_norm'});
        com = sprintf('EEG = brainbeats_process(EEG,''analysis'',''features'',%s%s%s,''hrv_features'',%s,''eeg_features'',%s%s,''parpool'',%g,''gpu'',%g%s);', ...
            heartArgs, heartOpts, cleanOpts, hrv_features, eeg_features, featOpts, params.parpool, params.gpu, common);
    otherwise   % 'rm_heart', 'coherence'
        com = sprintf('EEG = brainbeats_process(EEG,''analysis'',''%s'',%s%s%s%s%s);', ...
            params.analysis, heartArgs, heartOpts, cleanOpts, argstr(params, {'conf_thresh' 'keep_heart'}), common);
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
