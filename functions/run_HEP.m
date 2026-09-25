% RUN_HEP - Epoch EEG data around heartbeats for heartbeat-evoked potential
% (HEP), heartbeat-related spectral perturbation (HRSP) and phase coupling
% (HEPC) analyses.
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
%      (see CLEAN_EEG), then crop the epochs to the window.
%   5. HRSP and HEPC (params.hep_tf: all channels, in HEP.brainbeats.hrsp;
%      otherwise only the plotted channel) and the surrogate heartbeat
%      control (params.hep_surrogates: number of surrogates, results in
%      HEP.brainbeats.surrogate), see COMPUTE_HEP_TF. They are computed from
%      the continuous data at the heartbeats of the final epochs; with
%      params.clean_eeg, the continuous data get the same cleaning as the
%      epochs (linear map estimated from the epochs before and after it).
%   6. If params.hep_baseline is 'regression', regression-based baseline
%      correction (Alday, 2019; see BASELINE_REGRESSION). Otherwise no
%      baseline correction, since the pre-R-peak window contains activity
%      from the previous cardiac cycle.
%   7. Plot (params.vis_outputs: HEP, ERP image, HRSP/HEPC and, with
%      surrogates, the HEP against the surrogate range) and save
%      (params.save, as <filename>_HEP.set).
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

% Add the heart channel back, rescaled to the EEG range (mainly to check
% visually that the R-peak events align with the ECG)
if isfield(params,'keep_heart') && params.keep_heart && ~isempty(CARDIO)
    for iChan = 1:CARDIO.nbchan
        CARDIO.data(iChan,:) = rescale(CARDIO.data(iChan,:), prctile(EEG.data(:), 5)*2,  prctile(EEG.data(:), 95)*2);
    end
    if EEG.pnts ~= CARDIO.pnts
        EEG.data(end+1:end+CARDIO.nbchan,:) = CARDIO.data(:,1:end-1); % PPG can be one sample longer
    else
        EEG.data(end+1:end+CARDIO.nbchan,:) = CARDIO.data;
    end
    EEG.nbchan = EEG.nbchan + CARDIO.nbchan;
    for iChan = 1:CARDIO.nbchan
        EEG.chanlocs(end+1).labels = params.heart_channels{iChan};
    end
    EEG = eeg_checkset(EEG);
end

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

% Epoch around the R-peaks only (not around other events in the file), with
% 650 ms of padding on each side so that the same heartbeats can be used for
% the time-frequency measures (3 SD of a 5-cycle wavelet at 4 Hz = 600 ms)
pad = 650;   % ms
HEPwide = pop_epoch(EEG,{'R-peak'},(epochWin + [-pad pad])/1000,'epochinfo','yes');

% Remove bad epochs, run ICA, and remove bad components
if params.clean_eeg
    rawWide = HEPwide;   % to apply the same cleaning to the continuous data below
    [HEPwide, params] = clean_eeg(HEPwide,params);
    HEPwide.brainbeats.preprocessings.removed_eeg_trials = params.removed_eeg_trials;
    HEPwide.brainbeats.preprocessings.removed_eeg_components = params.removed_eeg_components;
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

% Regression-based baseline correction (Alday, 2019), on the EEG channels
% only. The corrected epochs are stored in HEP.data; the slopes and trial
% baselines are kept so the correction can be undone (see BASELINE_REGRESSION).
if isfield(params,'hep_baseline') && strcmpi(params.hep_baseline,'regression')
    if ~isfield(params,'hep_baseline_win') || isempty(params.hep_baseline_win)
        params.hep_baseline_win = [-300 -100];   % ms (Park & Blanke, 2019)
    end
    eegIdx = true(1,HEP.nbchan);
    if isfield(params,'heart_channels') && ~isempty(params.heart_channels)
        eegIdx = ~ismember(lower({HEP.chanlocs.labels}), lower(params.heart_channels));
    end
    fprintf('Regression-based baseline correction (%g to %g ms)... \n', params.hep_baseline_win)
    [HEP.data(eegIdx,:,:), beta, bl] = baseline_regression(HEP.data(eegIdx,:,:), HEP.times, params.hep_baseline_win);
    HEP.brainbeats.preprocessings.baseline_regression.window = params.hep_baseline_win;
    HEP.brainbeats.preprocessings.baseline_regression.channels = {HEP.chanlocs(eegIdx).labels};
    HEP.brainbeats.preprocessings.baseline_regression.beta = beta;
    HEP.brainbeats.preprocessings.baseline_regression.baseline = bl;
end
HEP.brainbeats.preprocessings.hep_window = epochWin;   % ms

% Channel shown in the plots: Fz, else Cz, else the first channel
elecName = 'Fz';
elecNum = find(strcmpi({HEP.chanlocs.labels}, elecName));
if isempty(elecNum)
    elecName = 'Cz';
    elecNum = find(strcmpi({HEP.chanlocs.labels}, elecName));
    if isempty(elecNum)
        elecNum = 1;
        elecName = HEP.chanlocs(elecNum).labels;
    end
end

% Time-frequency measures (HRSP, HEPC) and surrogate heartbeat control,
% computed from the continuous data at the same heartbeats as the HEP epochs.
% With EEG cleaning, the continuous data get the same (linear) cleaning as
% the epochs: channel interpolation and ICA component removal, estimated
% from the epochs before and after cleaning.
doTF = isfield(params,'hep_tf') && params.hep_tf;
nSurr = 0;
if isfield(params,'hep_surrogates') && ~isempty(params.hep_surrogates), nSurr = params.hep_surrogates; end
if doTF || nSurr > 0 || params.vis_outputs
    heartIdx = false(1,HEP.nbchan);
    if isfield(params,'heart_channels') && ~isempty(params.heart_channels)
        heartIdx = ismember(lower({HEP.chanlocs.labels}), lower(params.heart_channels));
    end
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
    % continuous sample of each epoch's heartbeat
    beatIdx = arrayfun(@(e) HEP.event(e.event(1)).beat, HEP.epoch);
    hepBeats = Rpeaks(beatIdx);
    bnd = [];
    if ~isempty(EEG.event) && any(strcmp({EEG.event.type}, 'boundary'))
        bnd = [EEG.event(strcmp({EEG.event.type}, 'boundary')).latency];
    end
    if doTF || nSurr > 0
        tfChans = find(~heartIdx);        % all EEG channels
    else
        tfChans = elecNum;                % plotted channel only
    end
    tfFreqs = 4:30;
    if isfield(params,'hep_tf_freqs') && ~isempty(params.hep_tf_freqs)
        tfFreqs = params.hep_tf_freqs(1):params.hep_tf_freqs(end);
    end
    if doTF, fprintf('Computing HRSP and HEPC on %g channels... \n', numel(tfChans)); end
    if nSurr > 0, fprintf('Surrogate heartbeat control (%g surrogates)... \n', nSurr); end
    [tf, surr] = compute_hep_tf(Xc(tfChans,:), EEG.srate, hepBeats, epochWin, ...
        struct('freqs',tfFreqs, 'nSurr',nSurr, 'boundaries',bnd, 'tf', doTF || params.vis_outputs, ...
        'hep_times',HEP.times));
    tf.channels = {HEP.chanlocs(tfChans).labels};
    if doTF
        HEP.brainbeats.hrsp = rmfield(tf, 'hep');
    end
    if nSurr > 0
        surr.channels = tf.channels;
        surr.times = tf.times; surr.hep_times = tf.hep_times; surr.freqs = tf.freqs;
        surr.hep.real = tf.hep;   % HEP of the same heartbeats (from the continuous data)
        HEP.brainbeats.surrogate = surr;
        fprintf('Surrogate control: %g/%g HEP points (channels x latencies) differ from the surrogates (FDR-corrected p < .05). \n', ...
            sum(surr.hep.p_fdr(:) < .05), numel(surr.hep.p_fdr))
    end
end


%% Plot heartbeat-evoked potentials (HEP) and heartbeat-related spectral
% perturbations (HRSP). HEP effects are usually reported 200-600 ms after
% the R-peak, over frontocentral electrodes.

if params.vis_outputs

    if params.vis_cleaning
        pop_eegplot(HEP,1,1,1);
        set(gcf,'Name','Final output','NumberTitle','Off','Toolbar','none','Menu','none');
        finish_figure(gcf)
    end

    % Trimmed-mean HEP at each electrode (click on one to enlarge it)
    options = { 'frames' HEP.pnts 'limits' [HEP.xmin HEP.xmax 0 0]*1000 ...
        'title' 'Heartbeat-evoked potentials (HEP)' 'chans' 1:HEP.nbchan ...
        'chanlocs' HEP.chanlocs 'ydir' 1 'legend' {'uV' 'Time (ms)'}};
    figure
    plottopo( trimmean(HEP.data,20,3), options{:} );
    try icadefs; set(gcf, 'color', BACKCOLOR); catch; end
    set(findall(gcf,'type','axes'),'fontSize',11,'fontweight','bold');
    set(gcf,'Toolbar','none','Menu','none');
    finish_figure(gcf)

    % Average HEP of all channels
    figure
    try icadefs; set(gcf, 'color', BACKCOLOR); catch; end
    subplot(2,1,1)

    % With the heart channel kept, plot it in thick black over the EEG
    if isfield(params,'keep_heart') && params.keep_heart
        hold on
        times = HEP.times;  % in ms
        heartIdx = ismember(lower({HEP.chanlocs.labels}), lower(params.heart_channels));
        for ch = find(~heartIdx)
            plot(times, trimmean(HEP.data(ch,:,:),20,3));
        end
        for ch = find(heartIdx)
            plot(times, trimmean(HEP.data(ch,:,:),20,3), 'k', 'LineWidth', 2.5);
        end
        xline(0, 'k--', 'LineWidth', 1.5);
        xlabel('Latency (ms)'); ylabel('Potential (µV)');
        title('Grand average HEP for each channel, including heart channel (thick black)');

    else
        % Butterfly plot with scalp topographies at selected latencies
        pop_timtopo(HEP, [HEP.times(1) HEP.times(end)], [-25 0 250 350 450], ...
            'Heartbeat-evoked potentials (HEP) - all electrodes','verbose','off');
    end

    % Single-trial HEPs over time (ERP image), with the surrogate 95% envelope
    subplot(2,1,2)
    pop_erpimage(HEP,1, elecNum,[],sprintf('Heartbeat-evoked potentials (HEP) over time for channel %s',elecName), ...
        10,1,{'R-peak'},[],'','yerplabel','\muV','erp','on','cbar','on' );
    colormap("parula")
    set(findall(gcf,'type','axes'),'fontSize',10,'fontweight','bold');
    set(gcf,'Name','HEP','NumberTitle','Off')
    finish_figure(gcf)

    if nSurr > 0
        iS = find(strcmp(surr.channels, elecName));
        figure('color','w'); hold on
        fill([surr.hep_times fliplr(surr.hep_times)], [surr.hep.null_lo(iS,:) fliplr(surr.hep.null_hi(iS,:))], ...
            [.8 .8 .8], 'EdgeColor','none');
        plot(surr.hep_times, tf.hep(iS,:), 'k', 'LineWidth', 2);
        sig = surr.hep.p_fdr(iS,:) < .05;
        plot(surr.hep_times(sig), tf.hep(iS,sig), 'r.', 'MarkerSize', 12);
        xline(0,'k--'); xlabel('Latency (ms)'); ylabel('Potential (µV)');
        title(sprintf('HEP at %s vs %g surrogate heartbeat trains (gray: 95%% of surrogates; red: FDR p < .05)', elecName, nSurr));
        set(gcf,'Name','HEP surrogate control','NumberTitle','Off')
        finish_figure(gcf)
    end

    % HRSP and HEPC at the same channel: 5-cycle Morlet wavelets, 4-30 Hz,
    % power in dB relative to its mean over the cardiac cycle (the pre-R-peak
    % window is not a neutral baseline: it holds the previous cycle and, once
    % smeared by the wavelets, the QRS). Contours: FDR-corrected p < .05
    % against the surrogates, if the surrogate control was computed.
    iT = find(strcmp(tf.channels, elecName));
    figure('color','w');
    subplot(1,2,1)
    imagesc(tf.times, tf.freqs, squeeze(tf.hrsp(iT,:,:))); axis xy
    lim = max(abs(tf.hrsp(iT,:)),[],'all'); set(gca,'CLim',[-lim lim]); colorbar
    xline(0,'k--'); xlabel('Latency (ms)'); ylabel('Frequency (Hz)');
    title(sprintf('HRSP at %s (dB)', elecName))
    subplot(1,2,2)
    imagesc(tf.times, tf.freqs, squeeze(tf.hepc(iT,:,:))); axis xy; colorbar
    xline(0,'k--'); xlabel('Latency (ms)'); ylabel('Frequency (Hz)');
    title(sprintf('HEPC at %s (pairwise phase consistency)', elecName))
    if nSurr > 0
        iS = find(strcmp(surr.channels, elecName));
        subplot(1,2,1); hold on
        contour(tf.times, tf.freqs, double(squeeze(surr.hrsp.p_fdr(iS,:,:)) < .05), [.5 .5], 'k', 'LineWidth', 1.5);
        subplot(1,2,2); hold on
        contour(tf.times, tf.freqs, double(squeeze(surr.hepc.p_fdr(iS,:,:)) < .05), [.5 .5], 'k', 'LineWidth', 1.5);
    end
    colormap("parula")
    set(gcf,'Name','Heartbeat-related spectral perturbations (HRSP) and phase coupling (HEPC)','NumberTitle','Off')
    finish_figure(gcf)

end

% Save next to the input file
if params.save
    newname = sprintf('%s_HEP.set', HEP.filename(1:end-4));
    pop_saveset(HEP,'filename',newname,'filepath',HEP.filepath);
end
