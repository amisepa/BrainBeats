% RUN_HEP - Epoch EEG data around heartbeats for heartbeat-evoked potential
% (HEP), heartbeat-related spectral perturbation (HRSP) and phase coupling
% (HRPC) analyses.
%
% Steps:
%   1. Epoch window: params.hep_window in ms (default [-300 600]), or
%      'adaptive' to end the epoch 50 ms before the next QRS of the 5%
%      shortest inter-beat intervals of this recording (400-1000 ms; for
%      within-subject analyses, since the window then differs across subjects).
%   2. Reject heartbeats followed by the next one before the epoch end + 50 ms
%      (the next QRS and its cardiac field artifact would fall in the epoch),
%      then outlier inter-beat intervals (Grubbs test).
%   3. Insert 'R-peak' events and epoch around them with 650 ms of padding
%      (boundary events are kept, so epochs spanning a data gap are dropped).
%   4. If params.clean_eeg, remove bad epochs (amplitude and high-frequency
%      noise over the analysis window, the latter on the broadband copy the
%      ICA is fitted on) and artifactual ICA components (see CLEAN_EEG; ICA
%      is fitted on a continuous 1-Hz high-passed copy), then crop the
%      epochs to the window.
%   5. The continuous data get the same cleaning as the epochs (linear map
%      estimated from the epochs before and after it), for the time-frequency
%      measures and surrogates, which use the continuous data at the
%      heartbeats of the final epochs.
%   6. params.hep_level: 'channels' (default), 'ics' or 'both'. On the
%      independent components ('ics', 'both'): HEP, HRSP and HRPC (and the
%      surrogate control) of all ICs: those of the cleaning ICA (after the
%      artifact components are removed), else of the ICA in the dataset,
%      else of an ICA run here (Picard). No component is selected: whether
%      a heartbeat-locked response can be separated from the cardiac field
%      artifact is not known. Results in HEP.brainbeats.ics, with each IC's
%      scalp map and ICLabel classification; .ic holds the IC numbers of the
%      decomposition, .removed_ics those removed by the cleaning.
%   7. With 'ref','csd', the surface Laplacian is applied to the epochs and
%      the continuous data (after the IC step, which uses the average-
%      referenced data).
%   8. QRS map of the cardiac field artifact (params.hep_cfa, default true;
%      CFA_MAP): its latency, scalp map and whether it stayed the same over
%      the recording (a warning if it changed: the head or body position,
%      or the electrodes, may have moved). HEP.brainbeats.cfa, with notes to
%      interpret it; with params.save, also <filename>_HEP_CFA.png and .txt.
%   9. Each heartbeat's voltage in the baseline window (params.hep_baseline_win,
%      default [-150 -50] ms) and its induced band power (theta, alpha, beta,
%      within params.hep_tf_freqs) before the QRS (epoch start to -50 ms),
%      for group-level or single-trial models (e.g. the baseline as a
%      regressor, Alday, 2019, or ongoing power as a covariate): stored in
%      HEP.brainbeats.prestim, not used to correct the HEP. With params.save,
%      the ROI values of each heartbeat are also written to
%      <filename>_HEP_beats.csv.
%   10. No baseline correction (as recommended by Steinfath et al., 2026),
%      since the pre-R-peak window contains activity from the previous
%      cardiac cycle. A regression-based correction within one recording
%      (Alday, 2019) would not change the average HEP; for data that pool
%      conditions, see BASELINE_REGRESSION.
%   11. On the channels ('channels', 'both'): HRSP and HRPC (COMPUTE_HEP_TF)
%      of all channels (HEP.brainbeats.hrsp) and of the ROI used for the
%      plots (params.hep_roi, default frontocentral: F1-F4, Fz, FC1-FC4,
%      FCz, C1-C4, Cz; their average signal; HEP.brainbeats.roi). With
%      params.hep_surrogates > 0, surrogate heartbeat control (surrogate
%      beats at random cardiac phases, their windows away from the epochs
%      rejected by the cleaning) of the HEP, HRSP and HRPC of all channels (HEP.brainbeats.surrogate) and of the ROI.
%      The HRSP plotted is that of the induced power, with the HEP (and the
%      cardiac field artifact) subtracted from each heartbeat (Cohen & Donner,
%      2013; params.hep_tf_power 'induced', default), or of all the power
%      ('total'); both are stored.
%   12. With params.keep_heart, the heart channel(s) are added to the epochs
%      (in their own units, after the EEG cleaning).
%   13. Plot (params.vis_outputs: IBI distribution, HEP of all channels with
%      scalp maps at 0 ms and at the largest ROI deflection at 150-350 ms and
%      after 350 ms, each labeled with its spatial correlation with the
%      cardiac field artifact at the QRS (r_CFA), ROI single-trial image, ROI HEP against the surrogates
%      with HRSP and HRPC; on the ICs, the map and HEP against the
%      surrogates of each IC, PLOT_HEP_ICS; QRS map and its stability) and
%      save (params.save, as <filename>_HEP.set).
%
% Usage:
%   HEP = run_HEP(EEG, CARDIO, params, Rpeaks)
%
% Inputs:
%   EEG    - continuous EEGLAB dataset (EEG channels only)
%   CARDIO - EEGLAB dataset with the heart channel(s), added back to the
%            output if params.keep_heart (can be empty)
%   params - BrainBeats parameters (see BRAINBEATS_PROCESS)
%   Rpeaks - heartbeat sample indices (in EEG samples)
%
% Output:
%   HEP    - epoched EEGLAB dataset time-locked to the heartbeats
%
% References:
%   Candia-Rivera et al. (2021). The role of EEG reference in the assessment
%   of functional brain-heart interplay: From methodology to user guidelines.
%   Journal of Neuroscience Methods.
%   Park & Blanke (2019). Heartbeat-evoked cortical responses: Underlying
%   mechanisms, functional roles, and methodological considerations.
%   NeuroImage.
%   Alday (2019). How much baseline correction do we need in ERP research?
%   Psychophysiology.
%   Steinfath et al. (2026). Heartbeat-evoked responses in M/EEG: a
%   systematic review of methods with suggestions for analysis and
%   reporting. Psychophysiology.
%   Lee et al. (2024). Heartbeat-related spectral perturbation of
%   electroencephalogram reflects dynamic interoceptive attention states in
%   the trial-by-trial classification analysis. NeuroImage (HRSP, 5-20 Hz).
%   Cohen & Donner (2013). Midfrontal conflict-related theta-band power
%   reflects neural oscillations that predict behavior. Journal of
%   Neurophysiology.
%   Dirlich et al. (1997). Cardiac field effects on the EEG.
%   Electroencephalography and Clinical Neurophysiology.
%
% Copyright (C) - Cedric Cannard, 2023

function HEP = run_HEP(EEG, CARDIO, params, Rpeaks)

% Epoch window (ms). The pre-R-peak part is fixed; the end is fixed
% (default) or set from this recording's heart rate ('adaptive').
Rpeaks = Rpeaks(:)';
IBI = [diff(Rpeaks) NaN] / EEG.srate * 1000;    % ms, from each beat to the next (NaN for the last)
if ~isfield(params,'hep_window') || isempty(params.hep_window)
    params.hep_window = [-300 600];
end
qrsMargin = 50;   % ms: the QRS starts ~50 ms before the R-peak
if ischar(params.hep_window) || isstring(params.hep_window)
    % 'adaptive': end the epoch 50 ms before the next QRS of the 5% shortest
    % inter-beat intervals, so that few beats are lost at any heart rate
    winEnd = floor((prctile(IBI,5) - qrsMargin)/10)*10;
    winEnd = min(max(winEnd, 400), 1000);
    epochWin = [-300 winEnd];
    fprintf('Adaptive HEP window: %g to %g ms (5th percentile of the inter-beat intervals: %.0f ms). \n', epochWin, prctile(IBI,5))
    warning(['The adaptive HEP window depends on the heart rate of this recording. It is meant for ' ...
        'within-subject analyses: for group analyses, use the same fixed window for all subjects ' ...
        '(e.g. ''hep_window'',[-300 600]) or crop all subjects to the shortest one.'])
else
    epochWin = params.hep_window;
end
params.hep_window = epochWin;   % exported (ms)

% Reject beats followed by the next one before the end of the epoch (+ the
% QRS margin), so that no epoch contains the next heartbeat's QRS and
% cardiac field artifact (see Candia-Rivera et al. 2021; Park & Blanke 2019)
minIBI = epochWin(2) + qrsMargin;
keep = ~(IBI < minIBI);
if any(~keep)
    fprintf('Removing %g/%g heartbeats followed by the next one within %g ms (epoch end + %g ms). \n', ...
        sum(~keep), length(IBI), minIBI, qrsMargin)
end

% Remove remaining outlier inter-beat intervals (Grubbs test on those kept)
kept = find(keep & ~isnan(IBI));
outliers = kept(isoutlier(IBI(kept),'grubbs'));
if ~isempty(outliers)
    fprintf('Removing %g outlier trials with the following interbeat intervals (IBI): \n', length(outliers))
    fprintf('   - IBI = %g (ms) \n', IBI(outliers))
    keep(outliers) = false;
end
fprintf('%g/%g heartbeats kept for HEP epochs. \n', sum(keep), length(Rpeaks))
nShort = sum(IBI < minIBI);
if nShort > .25*length(Rpeaks)
    warning(['%.0f%% of the heartbeats are followed by the next one within %g ms: the HEP comes from the slower ' ...
        'heartbeats of this recording, a subset that depends on the heart rate.'], 100*nShort/length(Rpeaks), minIBI)
end
beatCounts = struct('total',length(Rpeaks), 'short_ibi',nShort, 'outlier_ibi',length(outliers), ...
    'median_ibi_all',median(IBI,'omitnan'), 'median_ibi_kept',median(IBI(keep),'omitnan'));
allBeats = Rpeaks;   % whole heartbeat train, for the surrogates
Rpeaks = Rpeaks(keep);
IBI = IBI(keep);

% Add the heartbeat events after the existing ones
nEv = length(EEG.event);
urevents = num2cell(nEv+1:nEv+length(Rpeaks));
evt = num2cell(Rpeaks);
types = repmat({'R-peak'},1,length(evt));
[EEG.event(1,nEv+1:nEv+length(Rpeaks)).latency] = evt{:};
[EEG.event(1,nEv+1:nEv+length(Rpeaks)).type] = types{:};
[EEG.event(1,nEv+1:nEv+length(Rpeaks)).urevent] = urevents{:};
beatNum = num2cell(1:length(Rpeaks));
[EEG.event(1,nEv+1:nEv+length(Rpeaks)).beat] = beatNum{:};   % to find each epoch's heartbeat again
EEG = eeg_checkset(EEG, 'eventconsistency');   % sort events by latency

% Inter-beat interval distribution against the rejection threshold
if params.vis_outputs
    figure('color','w'); histfit(IBI(~isnan(IBI))); hold on
    plot([1 1]*minIBI,ylim,'--r','linewidth',2)
    title('Interbeat intervals (IBI) of the heartbeats kept');
    xlabel('Time (ms)'); ylabel('Number of IBIs')
    legend('', '', sprintf('Minimum IBI (epoch end + %g ms)', qrsMargin))
    try icadefs; set(gcf, 'color', BACKCOLOR); catch; end     % EEGLAB background color
    set(gcf,'Name','Inter-beat intervals (IBI) distribution','NumberTitle','Off','Toolbar','none','Menu','none')
    set(findall(gcf,'type','axes'),'fontSize',11,'fontweight','bold');
    finish_figure(gcf)
end

% Analysis settings
heartLabels = {};
if isfield(params,'heart_channels') && ~isempty(params.heart_channels)
    heartLabels = cellstr(params.heart_channels);
end
level = 'channels';
if isfield(params,'hep_level') && ~isempty(params.hep_level), level = params.hep_level; end
doChans = any(strcmp(level, {'channels' 'both'}));
doICs = any(strcmp(level, {'ics' 'both'}));
useCSD = params.clean_eeg && isfield(params,'ref') && strcmp(params.ref,'csd');
useGEDAI = isfield(params,'clean_method') && strcmpi(params.clean_method,'gedai');
nSurr = 0;
if isfield(params,'hep_surrogates') && ~isempty(params.hep_surrogates), nSurr = params.hep_surrogates; end
surrStat = 'tmax';   % statistic shown in the plots (tmax and FDR always stored; TFCE when chosen)
if isfield(params,'hep_surr_stat') && ~isempty(params.hep_surr_stat), surrStat = lower(params.hep_surr_stat); end
tfFreqs = 4:30;
if isfield(params,'hep_tf_freqs') && ~isempty(params.hep_tf_freqs)
    tfFreqs = params.hep_tf_freqs(1):params.hep_tf_freqs(end);
end
tfPower = 'induced';   % HRSP plotted (both are stored)
if isfield(params,'hep_tf_power') && ~isempty(params.hep_tf_power), tfPower = params.hep_tf_power; end
doCFA = ~isfield(params,'hep_cfa') || isempty(params.hep_cfa) || params.hep_cfa;
blWin = [-150 -50];    % ms: ends before the QRS onset
if isfield(params,'hep_baseline_win') && ~isempty(params.hep_baseline_win), blWin = params.hep_baseline_win; end

% Epoch around the R-peaks only (not around other events in the file), with
% 650 ms of padding on each side so that the same heartbeats can be used for
% the time-frequency measures (3 SD of a 5-cycle wavelet at 4 Hz = 600 ms)
pad = 650;   % ms
HEPwide = pop_epoch(EEG,{'R-peak'},(epochWin + [-pad pad])/1000,'epochinfo','yes');

% ICA is fitted on a 1-Hz high-passed copy when the data were high-passed
% below 1 Hz (slow drifts degrade ICA and ICLabel)
hasICA = ~isempty(HEPwide.icaweights);
needICA = (params.clean_eeg && ~useGEDAI) || (doICs && (useGEDAI || ~params.clean_eeg) && ~hasICA);
if needICA && isfield(params,'ica_cont')
    % The continuous copy made by CLEAN_EEG, for the ICA fit only: a lax ASR
    % ('ica_asr', cutoff in SD, default 80; 0 = off) removes large bursts
    % but keeps the eye activity; optional epoching and resampling
    C = params.ica_cont;
    if params.clean_eeg
        % its epochs at the same heartbeats (before ASR), for the high-
        % frequency metric of the bad-epoch detection (FIND_BADTRIALS): the
        % data are low-passed at 30 Hz, the copy is not
        Ce = C; Ce.event = EEG.event; Ce = eeg_checkset(Ce, 'eventconsistency');
        Ce = pop_epoch(Ce,{'R-peak'},(epochWin + [-pad pad])/1000,'epochinfo','yes');
        if Ce.trials == HEPwide.trials && Ce.pnts == HEPwide.pnts && Ce.nbchan == HEPwide.nbchan
            params.hf_source = Ce.data;
        end
        clear Ce
    end
    if isfield(params,'ica_asr') && params.ica_asr > 0
        C0 = C;
        try
            C = clean_asr(C, params.ica_asr);
            fprintf('ASR (cutoff %g SD) on the ICA copy: %.1f%% of the samples changed. \n', params.ica_asr, ...
                100*mean(any(abs(C.data - C0.data) > 1e-6*max(abs(C0.data(:))), 1)))
        catch ME
            warning('ASR on the ICA copy failed (%s): ICA is fitted on the copy without ASR.', ME.message)
            C = C0;
        end
        clear C0
    end
    if isfield(params,'ica_epochs') && params.ica_epochs
        % epoched like the HEP data, so that the bad epochs are removed from it too
        C.event = EEG.event; C = eeg_checkset(C, 'eventconsistency');
        C = pop_epoch(C,{'R-peak'},(epochWin + [-pad pad])/1000,'epochinfo','yes');
    end
    if isfield(params,'ica_srate') && params.ica_srate < C.srate
        C = pop_resample(C, params.ica_srate);
    elseif isfield(params,'ica_srate')
        fprintf('''ica_srate'' (%g Hz) is not below the sampling rate (%g Hz): the ICA copy is not resampled. \n', params.ica_srate, C.srate);
    end
    params.ica_source = C;
    clear C
elseif needICA && isfield(params,'highpass') && ~isempty(params.highpass) && params.highpass < 1
    causal = isfield(params,'filttype') && strcmpi(params.filttype,'causal');
    EEG1 = pop_eegfiltnew(EEG,'locutoff',1,'minphase',causal);
    params.ica_source = pop_epoch(EEG1,{'R-peak'},(epochWin + [-pad pad])/1000,'epochinfo','yes');
    clear EEG1
end

% Remove bad epochs, run ICA, and remove bad components
rawWide = HEPwide;   % to apply the same cleaning to the continuous data below
if params.clean_eeg && isfield(params,'heart_removal') && strcmpi(params.heart_removal,'ica_ecg')
    params.ecg_epochs = heart_epochs(CARDIO, HEPwide, Rpeaks, params);
end
if params.clean_eeg && isfield(params,'heart_removal') && strcmpi(params.heart_removal,'ica_caracas')
    params.caracas_source = EEG;   % CARACAS detects R waves on continuous IC time courses
end
if params.clean_eeg
    [HEPwide, params] = clean_eeg(HEPwide,params);
    HEPwide.brainbeats.preprocessings.removed_eeg_trials = params.removed_eeg_trials;
    HEPwide.brainbeats.preprocessings.removed_eeg_components = params.removed_eeg_components;
    if isfield(params,'ica_report')
        HEPwide.brainbeats.preprocessings.ica = params.ica_report;
        params = rmfield(params,'ica_report');
    end
elseif doICs && ~hasICA
    % No cleaning: ICA (Picard) and ICLabel only for the IC-level analysis
    src = HEPwide;
    if isfield(params,'ica_source'), src = params.ica_source; end
    if ~exist('picard','file'), plugin_askinstall('picard','picard',1); end
    dataRank = data_rank(src.data);
    src = pop_runica(src,'icatype','picard','maxiter',400,'mode','standard','pca',dataRank);
    if ~all(isfinite([src.icaweights(:); src.icasphere(:)]))
        warning('Picard did not converge (non-finite weights): ICA rerun with extended Infomax.')
        src = pop_runica(src,'icatype','runica','extended',1,'pca',dataRank);
    end
    HEPwide.icaweights = src.icaweights; HEPwide.icasphere = src.icasphere;
    HEPwide.icawinv = src.icawinv; HEPwide.icachansind = src.icachansind;
    HEPwide = eeg_checkset(HEPwide);
    clear src
end
if ~isfield(params,'removed_eeg_trials'), params.removed_eeg_trials = []; end
for f = {'ica_source' 'ica_cont' 'ecg_epochs' 'caracas_source' 'hf_source'}
    if isfield(params,f{1}), params = rmfield(params,f{1}); end
end

% HEP epochs: the padded epochs cropped to the analysis window
iWin = find(HEPwide.times >= epochWin(1) & HEPwide.times < epochWin(2));
HEP = pop_select(HEPwide, 'point', [iWin(1) iWin(end)]);

% Only keep the time-locking R-peak (latency 0) when an epoch holds several
% events (e.g. another event of the file)
count = 0;
for iEv = 1:length(HEP.epoch)
    lat = HEP.epoch(iEv).eventlatency;
    if iscell(lat) && length(lat) > 1
        [~, i0] = min(abs(cell2mat(lat)));
        HEP.epoch(iEv).event = HEP.epoch(iEv).event(i0);
        HEP.epoch(iEv).eventlatency = HEP.epoch(iEv).eventlatency(i0);
        HEP.epoch(iEv).eventduration = HEP.epoch(iEv).eventduration(i0);
        HEP.epoch(iEv).eventtype = HEP.epoch(iEv).eventtype(i0);
        count = count+1;
    end
end
if count > 0
    fprintf('%g epochs contained more than 1 event: only their time-locking R-peak was kept in HEP.epoch. \n',count)
end
if isfield(HEP.epoch,'eventurevent')
    HEP.epoch = rmfield(HEP.epoch,'eventurevent');
end
HEP = eeg_checkset(HEP);
HEP.brainbeats.preprocessings.hep_window = epochWin;   % ms
beatCounts.bad_epochs = numel(params.removed_eeg_trials);
beatCounts.kept = HEP.trials;
HEP.brainbeats.preprocessings.hep_beats = beatCounts;

% Continuous sample of each epoch's heartbeat, and data discontinuities
beatIdx = arrayfun(@(e) HEP.event(e.event(1)).beat, HEP.epoch);
hepBeats = Rpeaks(beatIdx);
bnd = [];
if ~isempty(EEG.event) && any(strcmp({EEG.event.type}, 'boundary'))
    bnd = [EEG.event(strcmp({EEG.event.type}, 'boundary')).latency];
end
% Windows of the epochs rejected by the cleaning (bad data), and the data
% before the first and after the last epoch (never checked by the cleaning):
% no surrogate window falls on them. (The tutorial dataset has 5,000-uV
% artifacts in its first 2 s, before the first epoch: surrogates that
% sometimes included them widened the null.)
wSmp = round(epochWin/1000*EEG.srate);
excl = zeros(0,2);
if ~isempty(params.removed_eeg_trials)
    rb = arrayfun(@(e) epoch_beat(rawWide, e), params.removed_eeg_trials(:));
    excl = Rpeaks(rb(:))' + wSmp;                               % [start end] samples
end
ab = arrayfun(@(e) epoch_beat(rawWide, e), 1:rawWide.trials);
excl = [excl; 1 Rpeaks(min(ab))+wSmp(1)-1; Rpeaks(max(ab))+wSmp(2)+1 EEG.pnts];

% Continuous data with the same (linear) cleaning as the epochs: channel
% interpolation and ICA component removal, estimated from the epochs
% before and after cleaning
needCont = doChans || doICs;
if needCont
    Xc = double(EEG.data);
    if params.clean_eeg
        keptTrials = setdiff(1:rawWide.trials, params.removed_eeg_trials);
        Xr = reshape(double(rawWide.data(:,:,keptTrials)), rawWide.nbchan, []);
        Yc = reshape(double(HEPwide.data), HEPwide.nbchan, []);
        T = (Yc*Xr') * pinv(Xr*Xr');
        err = norm(Yc - T*Xr, 'fro') / norm(Yc, 'fro');
        if err > 1e-3
            warning('The EEG cleaning of the epochs is not reproduced exactly on the continuous data (relative error %.1e).', err)
        end
        before = cfa_amplitude(Xc, hepBeats, EEG.srate);
        Xc = T * Xc;
        clear Xr Yc
        if ~useGEDAI    % GEDAI cleaned the continuous data: reported by CLEAN_EEG
            S = cfa_report(before, cfa_amplitude(Xc, hepBeats, EEG.srate), 'EEG cleaning');
            for fn = fieldnames(S)', HEP.brainbeats.preprocessings.(fn{1}) = S.(fn{1}); end
        end
    end
end

% HEP, HRSP, HRPC (and surrogate control) of all independent components
if doICs
    removedICs = [];
    if isfield(params,'removed_eeg_components'), removedICs = params.removed_eeg_components; end
    HEP.brainbeats.ics = hep_all_ics(HEPwide, Xc, hepBeats, epochWin, EEG.srate, struct('freqs',tfFreqs, ...
        'nSurr',nSurr, 'allBeats',allBeats, 'boundaries',bnd, 'exclude',excl, 'hep_times',HEP.times, ...
        'stat',surrStat), removedICs);
end

% Surface Laplacian (after the ICA-based steps, which need average-
% referenced data), on the epochs and the continuous data alike
if useCSD
    [HEP, C] = apply_csd(HEP, heartLabels);
    if ~isempty(C) && needCont
        Xc = C * Xc;
    end
end

% Units of the EEG after the reference (surface Laplacian: uV/cm^2)
unitTxt = 'uV';
if ischar(HEP.ref) && strcmpi(HEP.ref, 'csd-transform'), unitTxt = 'uV/cm^2'; end
beatTime = (hepBeats(:)' - 1) / EEG.srate;      % s; sample 1 is t = 0

% QRS map of the cardiac field artifact, and whether it changed during the
% recording (in the output's reference, after the cleaning). Not with PPG
% beats that were not shifted by the pulse arrival time: the QRS is then
% outside the epochs' first 100 ms.
pre = struct();
if isfield(EEG,'brainbeats') && isfield(EEG.brainbeats,'preprocessings'), pre = EEG.brainbeats.preprocessings; end
if doCFA && strcmpi(params.heart_signal,'ppg') && ~isfield(pre,'ppg_transit')
    fprintf('QRS map not computed: the epochs are time-locked to the PPG pulses without the pulse arrival time. \n');
    doCFA = false;
end
if doCFA
    % reference and heart-artifact removal as applied
    hr = 'none';
    if params.clean_eeg && isfield(params,'heart_removal') && ~isempty(params.heart_removal)
        hr = params.heart_removal;
        if (strcmp(hr,'ecg_regression') && ~isfield(pre,'heart_regression')) || (strcmp(hr,'ssp') && ~isfield(pre,'heart_ssp'))
            hr = sprintf('none (%s was not applied)', hr);
        end
    end
    if params.clean_eeg
        info = {sprintf('Reference: %s. Heart-artifact removal: %s%s.', char(string(HEP.ref)), hr, ...
            repmat(' (so this is the field left after the removal)', 1, ~startsWith(hr,'none')))};
    else
        info = {sprintf('No EEG cleaning (reference of the input data: %s).', char(string(HEP.ref)))};
    end
    try
        cfa = cfa_map(HEP.data, HEP.times, beatTime, struct('info',{info},'unit',unitTxt));
        cfa.channels = {HEP.chanlocs.labels};
        HEP.brainbeats.cfa = cfa;
        fprintf('QRS map (cardiac field artifact): %g ms, GFP %.2f %s. %s \n', cfa.latency, cfa.gfp, unitTxt, cfa.summary);
        if cfa.changed
            warning('BrainBeats:cfa_changed', 'QRS map (cardiac field artifact) %s See HEP.brainbeats.cfa.notes.', cfa.summary)
        end
    catch ME
        warning('QRS map not computed: %s', ME.message)
        doCFA = false;
    end
end

% Each heartbeat's voltage in the baseline window (Alday's regressor; band
% power before the R-peak is added with the time-frequency measures below)
iBl = HEP.times >= blWin(1) & HEP.times <= blWin(2);
if ~any(iBl)
    warning('The baseline window [%g %g] ms is outside the epochs: no baseline values exported.', blWin)
end
prestim = struct('baseline_window', blWin, ...
    'baseline', reshape(mean(double(HEP.data(:,iBl,:)), 2), HEP.nbchan, []), ...
    'channels', {{HEP.chanlocs.labels}}, 'beat_time', beatTime, 'ibi', IBI(beatIdx), 'unit', unitTxt);

% Channel ROI (default frontocentral, where HEP effects are usually reported)
labels = {HEP.chanlocs.labels};
roiLabels = {'F1' 'F2' 'F3' 'F4' 'Fz' 'FC1' 'FC2' 'FC3' 'FC4' 'FCz' 'C1' 'C2' 'C3' 'C4' 'Cz'};
if isfield(params,'hep_roi') && ~isempty(params.hep_roi), roiLabels = cellstr(params.hep_roi); end
roi = find(ismember(lower(labels), lower(roiLabels)));
if isempty(roi)
    roi = find(ismember(lower(labels), {'fz' 'cz'}), 1);
    if isempty(roi), roi = 1; end
    warning('None of the ROI channels (%s) is in the data: using %s.', strjoin(roiLabels, ' '), labels{roi})
end
HEP.brainbeats.preprocessings.hep_roi = labels(roi);

% Channels: time-frequency measures (HRSP, HRPC) and surrogate control, of
% the ROI average (plots) and of all channels
tfOpts = struct('freqs',tfFreqs, 'nSurr',nSurr, 'allBeats',allBeats, ...
    'boundaries',bnd, 'exclude',excl, 'hep_times',HEP.times, 'stat',surrStat);
tfMain = []; surrMain = [];
if doChans
    sig = mean(Xc(roi,:), 1);
    if nSurr > 0, fprintf('Surrogate heartbeat control (%g surrogates)... \n', nSurr); end
    [tfMain, surrMain] = compute_hep_tf(sig, EEG.srate, hepBeats, epochWin, tfOpts);
    HEP.brainbeats.roi = struct('channels',{labels(roi)}, 'hep',tfMain.hep, 'hep_times',tfMain.hep_times, ...
        'tf',rmfield(tfMain,'hep'), 'surrogate',surrMain);
end
if doChans
    fprintf('Computing HRSP and HRPC on %g channels... \n', size(Xc,1));
    % with each heartbeat's induced band power before the QRS (epoch start
    % to -50 ms; theta is smeared over the cardiac cycle by its long wavelets)
    tfOpts.prestim = struct('win',[epochWin(1) -50], 'bands',[4 8; 8 13; 13 31], 'names',{{'theta' 'alpha' 'beta'}});   % [low high) Hz
    [tf, surr] = compute_hep_tf(Xc, EEG.srate, hepBeats, epochWin, tfOpts);
    tf.channels = labels;
    [ok, loc] = ismember(hepBeats, tf.beats);
    P = nan(size(tf.prestim.power,1), size(tf.prestim.power,2), numel(hepBeats));
    P(:,:,ok) = tf.prestim.power(:,:,loc(ok));
    prestim.power_window = tf.prestim.win;
    prestim.bands = tf.prestim.names;
    prestim.band_limits = tf.prestim.bands;
    prestim.power = P;            % channels x bands x heartbeats (uV^2)
    HEP.brainbeats.hrsp = rmfield(tf, {'hep' 'prestim'});
    if nSurr > 0
        surr.channels = tf.channels;
        surr.times = tf.times; surr.hep_times = tf.hep_times; surr.freqs = tf.freqs;
        surr.hep.real = tf.hep;   % HEP of the same heartbeats (from the continuous data)
        HEP.brainbeats.surrogate = surr;
        nT = NaN; if isfield(surr.hep,'p_tfce'), nT = sum(surr.hep.p_tfce(:) < .05); end
        fprintf(['Surrogate control (%d surrogate HEPs): %g/%g HEP points (channels x latencies) differ from the surrogates ' ...
            '(p < .05; max statistic: %g, FDR: %g, TFCE: %g). \n'], surr.nSurrHep, sum(surr.hep.(['p_' surrStat])(:) < .05), ...
            numel(surr.hep.p_fdr), sum(surr.hep.p_tmax(:) < .05), sum(surr.hep.p_fdr(:) < .05), nT)
    end
end
clear Xc
prestim.notes = {['Values of each heartbeat before its QRS, for group-level or single-trial models; not used to ' ...
    'correct the HEP (within one recording, a regressor centered on its mean cannot change the average HEP).'] ...
    sprintf('baseline: mean voltage at %g to %g ms (%s, channels x heartbeats): the regressor of a regression-based baseline correction (Alday, 2019).', blWin, unitTxt) ...
    ['power: induced power (the HEP subtracted from each heartbeat) averaged over the latencies of power_window and ' ...
    'the frequencies of each band ((' unitTxt ')^2, channels x bands x heartbeats; Morlet wavelets of 5 cycles: theta is ' ...
    'smeared over the cardiac cycle). NaN: band outside hep_tf_freqs, or heartbeat too close to a data edge.'] ...
    'beat_time: time of each heartbeat in the recording (s); ibi: interval to the next heartbeat (ms).'};
HEP.brainbeats.prestim = prestim;

% Heart channel(s), in their own units, epoched at the same heartbeats
nEEG = HEP.nbchan;
if isfield(params,'keep_heart') && params.keep_heart && ~isempty(CARDIO)
    hepIdx = round(HEP.times/1000*EEG.srate);
    Hd = zeros(CARDIO.nbchan, HEP.pnts, HEP.trials);
    for e = 1:HEP.trials
        Hd(:,:,e) = CARDIO.data(:, hepBeats(e) + hepIdx);
    end
    HEP.data(end+1:end+CARDIO.nbchan,:,:) = Hd;
    HEP.nbchan = HEP.nbchan + CARDIO.nbchan;
    for iChan = 1:CARDIO.nbchan
        HEP.chanlocs(end+1).labels = params.heart_channels{iChan};
    end
    HEP = eeg_checkset(HEP);
end


% Save next to the input file (before the plots: a plotting error must not
% lose the results)
if params.save
    [name, folder] = output_file(HEP, '_HEP.set');
    saved = pop_saveset(HEP,'filename',name,'filepath',folder);
    base = fullfile(folder, name(1:end-4));
    % ROI values of each heartbeat (pre-R baseline and band power)
    try
        P = HEP.brainbeats.prestim;
        sfx = '_uV'; if ~strcmp(P.unit,'uV'), sfx = '_uVcm2'; end
        T = table((1:HEP.trials)', P.beat_time(:), P.ibi(:), mean(P.baseline(roi,:),1)', ...
            'VariableNames', {'heartbeat' 'time_s' 'ibi_ms' ['baseline_roi' sfx]});
        if isfield(P,'power')
            for b = 1:numel(P.bands)
                T.([P.bands{b} '_roi' sfx '_sq']) = reshape(mean(P.power(roi,b,:),1), [], 1);
            end
        end
        writetable(T, [base '_beats.csv']);
    catch ME
        warning('Heartbeat table not saved: %s', ME.message)
    end
    if doCFA
        write_lines([base '_CFA.txt'], HEP.brainbeats.cfa.notes);
    end
end
% QRS map figure (saved with the outputs)
if doCFA && (params.vis_outputs || params.save)
    try
        fig = plot_cfa_map(HEP.brainbeats.cfa, HEP.chanlocs(1:nEEG), params.vis_outputs, strrep(unitTxt, 'uV', '\muV'));
        if params.save
            exportgraphics(fig, [base '_CFA.png'], 'Resolution', 150, 'BackgroundColor', 'white');
        end
        if ~params.vis_outputs, close(fig); end
    catch ME
        warning('QRS map figure not made: %s', ME.message)
    end
end

%% Plots. HEP effects are usually reported 200-600 ms after the R-peak, over
% frontocentral electrodes.

if params.vis_outputs

    EEGonly = pop_select(HEP, 'channel', 1:nEEG);
    unit = strrep(unitTxt, 'uV', '\muV');
    if numel(roi) > 1
        sigName = sprintf('ROI (%d channels)', numel(roi));
    else
        sigName = labels{roi};
    end

    if params.vis_cleaning
        pop_eegplot(HEP,1,1,1);
        set(gcf,'Name','Final output','NumberTitle','Off','Toolbar','none','Menu','none');
        finish_figure(gcf)
    end

    % Trimmed-mean HEP at each electrode (click on one to enlarge it)
    options = { 'frames' EEGonly.pnts 'limits' [EEGonly.xmin EEGonly.xmax 0 0]*1000 ...
        'title' 'Heartbeat-evoked potentials (HEP, 20% trimmed mean)' 'chans' 1:EEGonly.nbchan ...
        'chanlocs' EEGonly.chanlocs 'ydir' 1 'legend' {unitTxt 'Time (ms)'}};
    figure
    plottopo( trimmean(EEGonly.data,20,3), options{:} );
    try icadefs; set(gcf, 'color', BACKCOLOR); catch; end
    set(findall(gcf,'type','axes'),'fontSize',11,'fontweight','bold');
    set(gcf,'Toolbar','none','Menu','none');
    finish_figure(gcf)

    % HEP of all electrodes with scalp maps (+ the heart signal, scaled to
    % the EEG range, if kept), and single heartbeats of the ROI/IC signal
    figure
    try icadefs; set(gcf, 'color', BACKCOLOR); catch; end
    subplot(2,1,1)
    if HEP.nbchan > nEEG
        % Butterfly plot with the heart signal (thick black, scaled to the
        % EEG range) to check that the heartbeats are aligned
        hold on
        erp = trimmean(EEGonly.data,20,3);
        plot(HEP.times, erp);
        yl = max(abs(erp(:)));
        for iChan = nEEG+1:HEP.nbchan
            h = trimmean(HEP.data(iChan,:,:),20,3);
            h = h - median(h);
            plot(HEP.times, h / max(abs(h)) * yl * 0.9, 'k', 'LineWidth', 2.5);
        end
        xline(0,'k--'); xlim(HEP.times([1 end])); box on
        xlabel('Latency (ms)'); ylabel(sprintf('Potential (%s)', unit));
        title(sprintf('Heartbeat-evoked potentials (HEP, 20%% trimmed mean) - all electrodes, and %s (black, scaled)', ...
            strjoin({HEP.chanlocs(nEEG+1:end).labels}, ', ')));
    else
        % Scalp maps at the R-peak (cardiac field artifact) and at the largest
        % deflection of the ROI HEP in an early (150-350 ms) and a late (350 ms
        % to the epoch end) window
        roiHEP = trimmean(mean(EEGonly.data(roi,:,:),1),20,3);
        mapTimes = 0;
        for w = [150 350; 350 EEGonly.times(end)]'
            iW = find(EEGonly.times >= w(1) & EEGonly.times <= w(2));
            if isempty(iW), continue, end
            [~, iMax] = max(abs(roiHEP(iW)));
            mapTimes(end+1) = EEGonly.times(iW(iMax)); %#ok<AGROW>
        end
        pop_timtopo(EEGonly, [EEGonly.times(1) EEGonly.times(end)], mapTimes, ...
            'Heartbeat-evoked potentials (HEP) - all electrodes','verbose','off');
        % Label each map with its spatial correlation with the map at the
        % largest global field power from -30 to 60 ms (the cardiac field
        % artifact at the QRS): how CFA-like the scalp distribution is. Same
        % average as the maps drawn (mean of the epochs).
        erpAll = mean(double(EEGonly.data), 3);
        iQ = find(EEGonly.times >= -30 & EEGonly.times <= 60);
        if numel(iQ) > 1 && size(erpAll,1) > 2
            [~, kQ] = max(std(erpAll(:,iQ), 0, 1));
            cfaMap = erpAll(:, iQ(kQ));
            labs = findall(gcf, 'Type', 'text');
            for iM = 1:numel(mapTimes)
                rCFA = corr(erpAll(:, EEGonly.times == mapTimes(iM)), cfaMap);
                for h = labs(:)'
                    p = get(h, 'Position');
                    if ischar(get(h,'String')) && strcmp(strtrim(get(h,'String')), num2str(mapTimes(iM),'%4.0f')) && p(1) == 0 && p(2) == 0.8
                        set(h, 'String', {sprintf('%g ms', mapTimes(iM)), sprintf('r_{CFA} = %.2f', rCFA)}, ...
                            'VerticalAlignment', 'bottom', 'Position', [0 0.62 p(3:end)]);
                    end
                end
            end
        end
    end
    subplot(2,1,2)
    roiTrials = reshape(mean(EEGonly.data(roi,:,:), 1), EEGonly.pnts, []);   % frames x heartbeats
    % color limits: 3 robust SDs (from the median absolute deviation) of the
    % 10-heartbeat moving average that erpimage draws. Its default, the
    % maximum, let a few artifact heartbeats flatten the image (and a 99th
    % percentile, on uncleaned data with artifacts: +/-32 uV for a 1-uV HEP);
    % artifact heartbeats now saturate the colors instead.
    sm = movmean(double(roiTrials), 10, 2);
    lim = 3 * 1.4826 * median(abs(sm(:) - median(sm(:))));
    cax = {};
    if lim > 0, cax = {'caxis', [-lim lim]}; end
    erpimage(roiTrials, [], HEP.times, ...
        sprintf('Single heartbeats: %s', sigName), 10, 1, 'erp', 'on', 'cbar', 'on', 'yerplabel', unit, cax{:});
    colormap("parula")
    set(findall(gcf,'type','axes'),'fontSize',10,'fontweight','bold');
    set(gcf,'Name','HEP','NumberTitle','Off')
    finish_figure(gcf)

    % ROI: HEP against the surrogates, HRSP and HRPC, with the shortest
    % intervals to the previous and next heartbeats of the beats used
    ibis = diff(allBeats) / EEG.srate * 1000;
    [~, loc] = ismember(hepBeats, allBeats);
    near = [min(ibis(loc(loc > 1) - 1)) min(ibis(loc(loc <= numel(ibis))))];
    plot_hep_tf(tfMain, surrMain, sigName, unit, near, tfPower);

    % ICs: map and HEP against the surrogates of each one
    if doICs && isfield(HEP.brainbeats,'ics')
        plot_hep_ics(HEP.brainbeats.ics);
    end
end

% The output points to its own file, not to the input file
if params.save
    HEP = saved;
end


%% Subfunctions

function b = epoch_beat(E, e)
% Heartbeat number (field 'beat') of the R-peak that time-locks epoch e
ev = E.epoch(e).event; lat = E.epoch(e).eventlatency; ty = E.epoch(e).eventtype;
if ~iscell(lat), lat = {lat}; ty = {ty}; end
lat = cell2mat(lat); lat(~strcmp(ty, 'R-peak')) = Inf;
[~, i0] = min(abs(lat));
b = E.event(ev(i0)).beat;


function H = heart_epochs(CARDIO, HEPwide, Rpeaks, params)
% Heart channel(s) filtered like the EEG and epoched like HEPwide (same
% heartbeats and latencies), for the ICs correlated with the ECG
if isempty(CARDIO)
    error('BrainBeats:input', "'heart_removal','ica_ecg' needs an ECG channel.")
end
causal = isfield(params,'filttype') && strcmpi(params.filttype,'causal');
C = pop_eegfiltnew(CARDIO,'locutoff',params.highpass,'hicutoff',params.lowpass,'minphase',causal);
lags = round(HEPwide.times/1000*HEPwide.srate);
H = zeros(C.nbchan, HEPwide.pnts, HEPwide.trials);
for e = 1:HEPwide.trials
    ev = HEPwide.epoch(e).event; lat = HEPwide.epoch(e).eventlatency; ty = HEPwide.epoch(e).eventtype;
    if ~iscell(lat), lat = {lat}; ty = {ty}; end
    isR = strcmp(ty, 'R-peak'); lat = cell2mat(lat);
    lat(~isR) = Inf; [~, i0] = min(abs(lat));
    H(:,:,e) = C.data(:, Rpeaks(HEPwide.event(ev(i0)).beat) + lags);
end


function out = hep_all_ics(HEPwide, Xc, hepBeats, epochWin, fs, opts, removed)
% HEP, HRSP, HRPC and surrogate control of the ICs left after the cleaning,
% from their activations on the continuous (cleaned) data; polarity: largest
% map value positive. removed: IC numbers removed by the cleaning.
out = [];
if isempty(HEPwide.icaweights)
    warning('No ICA decomposition: the IC measures need one (''clean_eeg'' or an ICA in the dataset).')
    return
end
chIdx = HEPwide.icachansind;
W = HEPwide.icaweights * HEPwide.icasphere;
maps = HEPwide.icawinv;
probs = [];
try probs = HEPwide.etc.ic_classification.ICLabel.classifications; catch; end
if isempty(probs) || size(probs,1) ~= size(W,1)
    try
        tmp = pop_iclabel(HEPwide, 'default');
        probs = tmp.etc.ic_classification.ICLabel.classifications;
    catch
        probs = [];
    end
end
[~, iMax] = max(abs(maps), [], 1);
sgn = sign(maps(sub2ind(size(maps), iMax, 1:size(maps,2))));
A = (W * Xc(chIdx,:)) .* sgn(:);
fprintf('Computing HEP, HRSP and HRPC of %d independent components... \n', size(A,1))
[tf, surr] = compute_hep_tf(A, fs, hepBeats, epochWin, opts);
out = struct('ic', setdiff(1:size(A,1)+numel(removed), removed)', 'removed_ics', removed(:), ...
    'maps', maps .* sgn, 'polarity', sgn(:), ...
    'chanlocs', HEPwide.chanlocs(chIdx), 'iclabel', probs, ...
    'iclabel_classes', {{'Brain' 'Muscle' 'Eye' 'Heart' 'Line Noise' 'Channel Noise' 'Other'}}, ...
    'hep', tf.hep, 'hep_times', tf.hep_times, 'tf', rmfield(tf,'hep'), 'surrogate', surr);


function [sig, name] = surr_sig(surr)
% HEP points that differ from the surrogates (p < .05) with the statistic
% chosen for the plots (surr.stat), and its name for the figure title
stat = 'fdr';
if isfield(surr,'stat') && isfield(surr.hep, ['p_' surr.stat]), stat = surr.stat; end
sig = surr.hep.(['p_' stat]) < .05;
names = struct('tfce','TFCE p < .05 (family-wise)', 'tmax','max-statistic p < .05 (family-wise)', 'fdr','FDR p < .05');
name = names.(stat);


function plot_hep_tf(tf, surr, sigName, unit, near, tfPower)
% HEP of the ROI signal against the surrogates (gray: 95% of the
% surrogates, red: p < .05 with the statistic of surr.stat, max statistic
% by default), HRSP (of the induced power by default, tfPower) and HRPC with
% the contours of p < .05 with the same statistic. The
% dashed lines mark +/-2 SD of the wavelets around the R-peak: within them,
% the time-frequency values include the QRS and its cardiac field artifact.
% The dotted lines do the same for the nearest previous and next heartbeats
% (near: shortest intervals to them, in ms).
if isempty(tf), return; end
figure('color','w','Name','HEP time-frequency','NumberTitle','Off');
try icadefs; set(gcf, 'color', BACKCOLOR); catch; end
hasS = ~isempty(surr);
t = tf.hep_times;

% HEP vs surrogates
subplot(2,1,1)
hold on
if hasS
    [sig, statName] = surr_sig(surr);
    sig = sig(1,:);
    yl = [min([surr.hep.null_lo(1,:) tf.hep(1,:)]) max([surr.hep.null_hi(1,:) tf.hep(1,:)])];
    yl = yl + [-1 1]*0.08*diff(yl);
    d = diff([false sig false]);
    on = find(d == 1); off = find(d == -1) - 1;
    for k = 1:numel(on)
        patch([t(on(k)) t(off(k)) t(off(k)) t(on(k))], yl([1 1 2 2]), [1 .85 .85], 'EdgeColor','none');
    end
    fill([t fliplr(t)], [surr.hep.null_lo(1,:) fliplr(surr.hep.null_hi(1,:))], [.8 .8 .8], 'EdgeColor','none', 'FaceAlpha', .8);
    ylim(yl)
end
plot(t, tf.hep(1,:), 'k', 'LineWidth', 2);
xline(0,'k--'); yline(0,'k:');
xlim(t([1 end])); box on
xlabel('Latency (ms)'); ylabel(sprintf('HEP (%s)', unit));
if hasS
    title(sprintf('HEP, %s: gray = 95%% of %d surrogate heartbeat trains, red = %s', ...
        sigName, surr.nSurrHep, statName))
else
    title(sprintf('HEP, %s', sigName))
end

% HRSP and HRPC
guide = 2 * tf.cycles ./ (2*pi*tf.freqs) * 1000;   % ms
hName = 'hrsp'; hLabel = 'HRSP';
if strcmp(tfPower, 'induced'), hName = 'hrsp_induced'; hLabel = 'HRSP, HEP subtracted'; end
maps = {squeeze(tf.(hName)(1,:,:)), squeeze(tf.hrpc(1,:,:))};
names = {sprintf('%s, %s (dB)', hLabel, sigName), sprintf('HRPC, %s (pairwise phase consistency)', sigName)};
for k = 1:2
    subplot(2,2,2+k); hold on
    M = maps{k};
    imagesc(tf.times, tf.freqs, M); axis xy tight
    if hasS
        % contours: p < .05 with the statistic of the HEP plot (max
        % statistic by default, family-wise over frequencies and latencies)
        if k == 1, S = surr.(hName); else, S = surr.hrpc; end
        pf = 'p_fdr';
        if isfield(surr,'stat') && isfield(S, ['p_' surr.stat]), pf = ['p_' surr.stat]; end
        contour(tf.times, tf.freqs, double(squeeze(S.(pf)(1,:,:)) < .05), [.5 .5], 'k', 'LineWidth', 1.5);
    end
    plot(-guide, tf.freqs, 'w--', guide, tf.freqs, 'w--', 'LineWidth', 1);
    plot(guide - near(1), tf.freqs, 'w:', near(2) - guide, tf.freqs, 'w:', 'LineWidth', 1.5);
    xline(0,'k--');
    if k == 1
        lim = max(abs(M(:)));
        if lim > 0, set(gca,'CLim',[-lim lim]); end
    elseif max(M(:)) > min(M(:))
        set(gca,'CLim',[min(M(:)) max(M(:))]);   % after the contour, which changes it
    end
    colorbar; box on
    xlim(tf.times([1 end])); ylim(tf.freqs([1 end]))
    xlabel('Latency (ms)'); ylabel('Frequency (Hz)'); title(names{k})
end
colormap("parula")
set(findall(gcf,'type','axes'),'fontSize',10,'fontweight','bold');
finish_figure(gcf)


function fig = plot_cfa_map(cfa, chanlocs, visible, unit)
% QRS map (cardiac field artifact), and its stability over the recording:
% map of the heartbeats before vs after each time point, and of each block
% vs the other heartbeats (r), with the same for the heartbeats in random
% order (what noise alone gives)
vis = {'off' 'on'};
fig = figure('Color','w','Name','QRS map (cardiac field artifact)','NumberTitle','Off', ...
    'Visible',vis{1+logical(visible)},'Position',[100 100 1000 480]);
subplot(1,2,1)
topoplot(cfa.map, chanlocs, 'electrodes', 'on');
lim = max(abs(cfa.map));
if lim > 0, set(gca,'CLim',[-lim lim]); end
colormap(fig, 'parula')          % as the HEP figures
colorbar
title({sprintf('QRS map: %g ms after the R-peak', cfa.latency), sprintf('GFP %.2f %s', cfa.gfp, unit)})
gray = [.55 .55 .55];
P = {cfa.split, cfa.blocks};
x = {cfa.split.minutes, arrayfun(@(b) mean(b.minutes), cfa.blocks)};
ttl = {sprintf('Stability: %s', strtok(cfa.summary, ':(')), ''};
lab = {'before vs after (min)', 'block vs the other heartbeats (min)'};
for k = 1:2
    subplot(2,2,2*k); hold on
    r = [P{k}.r]; m = [P{k}.r_random]; lo = [P{k}.r_random_lo];
    plot(x{k}, m, '-', 'Color', gray, 'LineWidth', 2);
    plot(x{k}, lo, '--', 'Color', gray, 'LineWidth', 1.5);
    plot(x{k}, r, 'ko-', 'MarkerFaceColor', 'k', 'LineWidth', 1.5);
    ylim([min([lo r .5]) - .05, 1]); box on
    xlabel(sprintf('Time in the recording: %s', lab{k})); ylabel('Map r')
    if k == 1
        legend({'random order (median)' 'random order (2.5th percentile)' 'recording order'}, ...
            'Location','southwest','Box','off')
    end
    title(ttl{k})
end
set(findall(fig,'type','axes'),'fontSize',10,'fontweight','bold');
finish_figure(fig)


function write_lines(file, lines)
% Write a cell array of text lines to a text file
fid = fopen(file, 'w');
if fid < 0
    warning('Could not write %s.', file); return
end
fprintf(fid, '%s\n', lines{:});
fclose(fid);
