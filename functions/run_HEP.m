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
%   4. If params.clean_eeg, remove bad epochs and artifactual ICA components
%      (see CLEAN_EEG; with a high-pass below 1 Hz, ICA is fitted on a 1-Hz
%      high-passed copy of the epochs), then crop the epochs to the window.
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
%      scalp map and ICLabel classification.
%   7. With 'ref','csd', the surface Laplacian is applied to the epochs and
%      the continuous data (after the IC step, which uses the average-
%      referenced data).
%   8. If params.hep_baseline is 'regression', regression-based baseline
%      correction (Alday, 2019; see BASELINE_REGRESSION). Otherwise no
%      baseline correction (default, as recommended by Steinfath et al.,
%      2026), since the pre-R-peak window contains activity from the
%      previous cardiac cycle.
%   9. On the channels ('channels', 'both'): HRSP and HRPC (COMPUTE_HEP_TF)
%      of all channels (HEP.brainbeats.hrsp) and of the ROI used for the
%      plots (params.hep_roi, default frontocentral: F1-F4, Fz, FC1-FC4,
%      FCz, C1-C4, Cz; their average signal; HEP.brainbeats.roi). With
%      params.hep_surrogates > 0, surrogate heartbeat control
%      (params.hep_surrogate_mode, 'shuffle' by default) of the HEP, HRSP
%      and HRPC of all channels (HEP.brainbeats.surrogate) and of the ROI.
%   10. With params.keep_heart, the heart channel(s) are added to the epochs
%      (in their own units, after the EEG cleaning).
%   11. Plot (params.vis_outputs: IBI distribution, HEP of all channels with
%      scalp maps, ROI/IC single-trial image, ROI/IC HEP against the
%      surrogates with HRSP and HRPC) and save (params.save, as
%      <filename>_HEP.set).
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
surrMode = 'shuffle';
if isfield(params,'hep_surrogate_mode') && ~isempty(params.hep_surrogate_mode), surrMode = params.hep_surrogate_mode; end
tfFreqs = 4:30;
if isfield(params,'hep_tf_freqs') && ~isempty(params.hep_tf_freqs)
    tfFreqs = params.hep_tf_freqs(1):params.hep_tf_freqs(end);
end

% Epoch around the R-peaks only (not around other events in the file), with
% 650 ms of padding on each side so that the same heartbeats can be used for
% the time-frequency measures (3 SD of a 5-cycle wavelet at 4 Hz = 600 ms)
pad = 650;   % ms
HEPwide = pop_epoch(EEG,{'R-peak'},(epochWin + [-pad pad])/1000,'epochinfo','yes');

% ICA is fitted on a 1-Hz high-passed copy when the data were high-passed
% below 1 Hz (slow drifts degrade ICA and ICLabel)
hasICA = ~isempty(HEPwide.icaweights);
needICA = (params.clean_eeg && ~useGEDAI) || (doICs && (useGEDAI || ~params.clean_eeg) && ~hasICA);
if needICA && isfield(params,'highpass') && ~isempty(params.highpass) && params.highpass < 1
    causal = isfield(params,'filttype') && strcmpi(params.filttype,'causal');
    EEG1 = pop_eegfiltnew(EEG,'locutoff',1,'minphase',causal);
    params.ica_source = pop_epoch(EEG1,{'R-peak'},(epochWin + [-pad pad])/1000,'epochinfo','yes');
    clear EEG1
end

% Remove bad epochs, run ICA, and remove bad components
rawWide = HEPwide;   % to apply the same cleaning to the continuous data below
if params.clean_eeg
    [HEPwide, params] = clean_eeg(HEPwide,params);
    HEPwide.brainbeats.preprocessings.removed_eeg_trials = params.removed_eeg_trials;
    HEPwide.brainbeats.preprocessings.removed_eeg_components = params.removed_eeg_components;
elseif doICs && ~hasICA
    % No cleaning: ICA (Picard) and ICLabel only for the IC-level analysis
    src = HEPwide;
    if isfield(params,'ica_source'), src = params.ica_source; end
    if ~exist('picard','file'), plugin_askinstall('picard','picard',1); end
    dataRank = sum(eig(cov(double(src.data(:,:)'))) > 1E-7);
    src = pop_runica(src,'icatype','picard','maxiter',400,'mode','standard','pca',dataRank);
    HEPwide.icaweights = src.icaweights; HEPwide.icasphere = src.icasphere;
    HEPwide.icawinv = src.icawinv; HEPwide.icachansind = src.icachansind;
    HEPwide = eeg_checkset(HEPwide);
    clear src
end
if ~isfield(params,'removed_eeg_trials'), params.removed_eeg_trials = []; end
if isfield(params,'ica_source'), params = rmfield(params,'ica_source'); end

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

% Continuous sample of each epoch's heartbeat, and data discontinuities
beatIdx = arrayfun(@(e) HEP.event(e.event(1)).beat, HEP.epoch);
hepBeats = Rpeaks(beatIdx);
bnd = [];
if ~isempty(EEG.event) && any(strcmp({EEG.event.type}, 'boundary'))
    bnd = [EEG.event(strcmp({EEG.event.type}, 'boundary')).latency];
end

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
        Xc = T * Xc;
        clear Xr Yc
    end
end

% HEP, HRSP, HRPC (and surrogate control) of all independent components
if doICs
    HEP.brainbeats.ics = hep_all_ics(HEPwide, Xc, hepBeats, epochWin, EEG.srate, struct('freqs',tfFreqs, ...
        'nSurr',nSurr, 'surrMode',surrMode, 'allBeats',allBeats, 'boundaries',bnd, 'hep_times',HEP.times));
end

% Surface Laplacian (after the ICA-based steps, which need average-
% referenced data), on the epochs and the continuous data alike
if useCSD
    [HEP, C] = apply_csd(HEP, heartLabels);
    if ~isempty(C) && needCont
        Xc = C * Xc;
    end
end

% Regression-based baseline correction (Alday, 2019). The corrected epochs
% are stored in HEP.data; the slopes and trial baselines are kept so the
% correction can be undone (see BASELINE_REGRESSION).
if isfield(params,'hep_baseline') && strcmpi(params.hep_baseline,'regression')
    if ~isfield(params,'hep_baseline_win') || isempty(params.hep_baseline_win)
        params.hep_baseline_win = [-150 -50];   % ms: ends before the QRS onset
    end
    fprintf('Regression-based baseline correction (%g to %g ms)... \n', params.hep_baseline_win)
    [HEP.data, beta, bl] = baseline_regression(HEP.data, HEP.times, params.hep_baseline_win);
    HEP.brainbeats.preprocessings.baseline_regression.window = params.hep_baseline_win;
    HEP.brainbeats.preprocessings.baseline_regression.channels = {HEP.chanlocs.labels};
    HEP.brainbeats.preprocessings.baseline_regression.beta = beta;
    HEP.brainbeats.preprocessings.baseline_regression.baseline = bl;
end

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
tfOpts = struct('freqs',tfFreqs, 'nSurr',nSurr, 'surrMode',surrMode, 'allBeats',allBeats, ...
    'boundaries',bnd, 'hep_times',HEP.times);
tfMain = []; surrMain = [];
if doChans
    sig = mean(Xc(roi,:), 1);
    if nSurr > 0, fprintf('Surrogate heartbeat control (%g surrogates, %s)... \n', nSurr, surrMode); end
    [tfMain, surrMain] = compute_hep_tf(sig, EEG.srate, hepBeats, epochWin, tfOpts);
    HEP.brainbeats.roi = struct('channels',{labels(roi)}, 'hep',tfMain.hep, 'hep_times',tfMain.hep_times, ...
        'tf',rmfield(tfMain,'hep'), 'surrogate',surrMain);
end
if doChans
    fprintf('Computing HRSP and HRPC on %g channels... \n', size(Xc,1));
    [tf, surr] = compute_hep_tf(Xc, EEG.srate, hepBeats, epochWin, tfOpts);
    tf.channels = labels;
    HEP.brainbeats.hrsp = rmfield(tf, 'hep');
    if nSurr > 0
        surr.channels = tf.channels;
        surr.times = tf.times; surr.hep_times = tf.hep_times; surr.freqs = tf.freqs;
        surr.hep.real = tf.hep;   % HEP of the same heartbeats (from the continuous data)
        HEP.brainbeats.surrogate = surr;
        fprintf('Surrogate control: %g/%g HEP points (channels x latencies) differ from the surrogates (FDR-corrected p < .05). \n', ...
            sum(surr.hep.p_fdr(:) < .05), numel(surr.hep.p_fdr))
    end
end
clear Xc

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


%% Plots. HEP effects are usually reported 200-600 ms after the R-peak, over
% frontocentral electrodes.

if params.vis_outputs

    EEGonly = pop_select(HEP, 'channel', 1:nEEG);
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
        'title' 'Heartbeat-evoked potentials (HEP)' 'chans' 1:EEGonly.nbchan ...
        'chanlocs' EEGonly.chanlocs 'ydir' 1 'legend' {'uV' 'Time (ms)'}};
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
        xlabel('Latency (ms)'); ylabel('Potential (\muV)');
        title(sprintf('Heartbeat-evoked potentials (HEP) - all electrodes, and %s (black, scaled)', ...
            strjoin({HEP.chanlocs(nEEG+1:end).labels}, ', ')));
    else
        pop_timtopo(EEGonly, [EEGonly.times(1) EEGonly.times(end)], [-25 0 250 350 450], ...
            'Heartbeat-evoked potentials (HEP) - all electrodes','verbose','off');
    end
    subplot(2,1,2)
    roiTrials = reshape(mean(EEGonly.data(roi,:,:), 1), EEGonly.pnts, []);   % frames x heartbeats
    unit = '\muV';
    erpimage(roiTrials, [], HEP.times, ...
        sprintf('Single heartbeats: %s', sigName), 10, 1, 'erp', 'on', 'cbar', 'on', 'yerplabel', unit);
    colormap("parula")
    set(findall(gcf,'type','axes'),'fontSize',10,'fontweight','bold');
    set(gcf,'Name','HEP','NumberTitle','Off')
    finish_figure(gcf)

    % ROI: HEP against the surrogates, HRSP and HRPC
    plot_hep_tf(tfMain, surrMain, sigName, unit);
end

% Save next to the input file
if params.save
    newname = sprintf('%s_HEP.set', HEP.filename(1:end-4));
    HEP = pop_saveset(HEP,'filename',newname,'filepath',HEP.filepath);   % the output then points to this file, not to the input file
end


%% Subfunctions

function out = hep_all_ics(HEPwide, Xc, hepBeats, epochWin, fs, opts)
% HEP, HRSP, HRPC and surrogate control of all ICs, from their activations
% on the continuous (cleaned) data; polarity: largest map value positive
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
out = struct('ic', (1:size(A,1))', 'maps', maps .* sgn, 'polarity', sgn(:), ...
    'chanlocs', HEPwide.chanlocs(chIdx), 'iclabel', probs, ...
    'iclabel_classes', {{'Brain' 'Muscle' 'Eye' 'Heart' 'Line Noise' 'Channel Noise' 'Other'}}, ...
    'hep', tf.hep, 'hep_times', tf.hep_times, 'tf', rmfield(tf,'hep'), 'surrogate', surr);


function plot_hep_tf(tf, surr, sigName, unit)
% HEP of the ROI signal against the surrogates (gray: 95% of the
% surrogates, red: FDR p < .05), HRSP and HRPC with their FDR contours. The
% dashed lines mark +/-2 SD of the wavelets around the R-peak: within them,
% the time-frequency values include the QRS and its cardiac field artifact.
if isempty(tf), return; end
figure('color','w','Name','HEP time-frequency','NumberTitle','Off');
try icadefs; set(gcf, 'color', BACKCOLOR); catch; end
hasS = ~isempty(surr);
t = tf.hep_times;

% HEP vs surrogates
subplot(2,1,1)
hold on
if hasS
    sig = surr.hep.p_fdr(1,:) < .05;
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
    title(sprintf('HEP, %s: gray = 95%% of %d surrogate heartbeat trains (%s), red = FDR p < .05', ...
        sigName, surr.nSurr, surr.mode))
else
    title(sprintf('HEP, %s', sigName))
end

% HRSP and HRPC
guide = 2 * tf.cycles ./ (2*pi*tf.freqs) * 1000;   % ms
maps = {squeeze(tf.hrsp(1,:,:)), squeeze(tf.hrpc(1,:,:))};
names = {sprintf('HRSP, %s (dB)', sigName), sprintf('HRPC, %s (pairwise phase consistency)', sigName)};
for k = 1:2
    subplot(2,2,2+k); hold on
    M = maps{k};
    imagesc(tf.times, tf.freqs, M); axis xy tight
    if hasS
        if k == 1, p = surr.hrsp.p_fdr; else, p = surr.hrpc.p_fdr; end
        contour(tf.times, tf.freqs, double(squeeze(p(1,:,:)) < .05), [.5 .5], 'k', 'LineWidth', 1.5);
    end
    plot(-guide, tf.freqs, 'w--', guide, tf.freqs, 'w--', 'LineWidth', 1);
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
