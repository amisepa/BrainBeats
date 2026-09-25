% RUN_HEP - Epoch EEG data around heartbeats for heartbeat-evoked potential
% (HEP) and heartbeat-related spectral perturbation (HRSP) analyses.
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
%   5. If params.hep_baseline is 'regression', regression-based baseline
%      correction (Alday, 2019; see BASELINE_REGRESSION). Otherwise no
%      baseline correction, since the pre-R-peak window contains activity
%      from the previous cardiac cycle.
%   6. Plot (params.vis_outputs; HEP, ERP image and HRSP on the padded
%      epochs) and save (params.save, as <filename>_HEP.set).
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
end

% Epoch around the R-peaks only (not around other events in the file), with
% 650 ms of padding on each side for the time-frequency decomposition
% (HRSP: 5-cycle wavelets at 5 Hz span 1 s)
pad = 650;   % ms (half a 5-cycle wavelet at 5 Hz + margin)
HEPwide = pop_epoch(EEG,{'R-peak'},(epochWin + [-pad pad])/1000,'epochinfo','yes');

% Remove bad epochs, run ICA, and remove bad components
if params.clean_eeg
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


%% Plot heartbeat-evoked potentials (HEP) and heartbeat-related spectral
% perturbations (HRSP). HEP effects are usually reported 200-600 ms after
% the R-peak, over frontocentral electrodes.

if params.vis_outputs

    if params.vis_cleaning
        pop_eegplot(HEP,1,1,1);
        set(gcf,'Name','Final output','NumberTitle','Off','Toolbar','none','Menu','none');
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

    % Single-trial HEPs over time (ERP image) at Fz, else Cz, else the first channel
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
    subplot(2,1,2)
    pop_erpimage(HEP,1, elecNum,[],sprintf('Heartbeat-evoked potentials (HEP) over time for channel %s',elecName), ...
        10,1,{'R-peak'},[],'','yerplabel','\muV','erp','on','cbar','on' );
    colormap("parula")
    set(findall(gcf,'type','axes'),'fontSize',10,'fontweight','bold');
    set(gcf,'Name','HEP','NumberTitle','Off')

    % HRSP and inter-trial coherence at the same channel (Lee et al., 2024):
    % 5-cycle Morlet wavelets, 5-20 Hz in 1-Hz steps, computed on the padded
    % epochs so the whole window is free of edge effects. Power is expressed
    % in dB relative to its mean over the whole epoch window (cardiac cycle),
    % since the pre-R-peak window is not a neutral baseline (it holds the
    % previous cycle and, once smeared by the wavelets, the QRS).
    % Bootstrap statistics, FDR-corrected.
    figure('color','w');
    tout = epochWin(1):10:epochWin(2);
    pop_newtimef(HEPwide,1,elecNum,[HEPwide.times(1) HEPwide.times(end)],5, ...
        'freqs',[5 20],'nfreqs',16,'freqscale','linear','timesout',tout, ...
        'baseline',epochWin,'plotphase','off','padratio',2, ...
        'alpha',0.05,'mcorrect','fdr','naccu',1000, ...
        'caption',sprintf('HRSP - channel %s (p = 0.05, FDR-corrected)',elecName));
    colormap("parula")
    set(gcf,'Toolbar','none','Menu','none');
    set(gcf,'Name','Heartbeat-related spectral perturbations (HRSP) and ITC','NumberTitle','Off')

end

% Save next to the input file
if params.save
    newname = sprintf('%s_HEP.set', HEP.filename(1:end-4));
    pop_saveset(HEP,'filename',newname,'filepath',HEP.filepath);
end
