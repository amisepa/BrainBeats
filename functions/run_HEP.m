% Process data for heartbeat-evoked potentials (HEP) analysis
%
% Following guidelines by:
%   Candia-Rivera et al. (2021). The role of EEG reference in the assessment
%   of functional brain–heart interplay: From methodology to user guidelines.
%   Journal of Neuroscience Methods.
%
%   Park and Blanke (2019). Heartbeat-evoked cortical responses: Underlying
%   mechanisms, functional roles, and methodological considerations.
%
% Cedric Cannard, 2023

function HEP = run_HEP(EEG, params, Rpeaks)

% Add heartbeat markers in the signals, taking into account
% existing events
nEv = length(EEG.event);
urevents = num2cell(nEv+1:nEv+length(Rpeaks));
evt = num2cell(Rpeaks);
types = repmat({'R-peak'},1,length(evt));
[EEG.event(1,nEv+1:nEv+length(Rpeaks)).latency] = evt{:};       % assign latencies
[EEG.event(1,nEv+1:nEv+length(Rpeaks)).type] = types{:};        % assign types
[EEG.event(1,nEv+1:nEv+length(Rpeaks)).urevent] = urevents{:};  % assign event index
EEG = eeg_checkset(EEG);

% Calculate inter-beat-intervals (IBI) from EEG markers and R peaks (to
% compare)
IBI = nan(length(EEG.event),1);
for iEv = 2:length(EEG.event)
    EEG.event(iEv).duration = ( [EEG.event(iEv).latency] - [EEG.event(iEv-1).latency] )/EEG.srate *1000;
    IBI(iEv,:) = ( [EEG.event(iEv).latency] - [EEG.event(iEv-1).latency] )/EEG.srate *1000;
end
IBI(1) = []; % remove first NaN value
% figure('color','w'); histogram(IBI); hold on;
% histogram(diff(Rpeaks)/EEG.srate *1000); legend('from EEG events','from R peaks');

% Remove trials with IBI<550 ms (see Candia-Rivera et al. 2021; Park & Blanke 2019)
shortTrials = find(IBI<550);
if ~isempty(shortTrials)
    warning('Removing %g trials with interbeat intervals (IBI) < 550 ms', length(shortTrials))
    IBI(shortTrials) = [];
    [EEG.event(shortTrials).type] = deal([]);
else
    fprintf('All trial have an interbeat interval (IBI) above 550 ms (minimum recommended). \n')
end

% Remove remaining outlier trials
outliers = find(isoutlier(IBI,'grubbs'));
if ~isempty(outliers)
    warning('Removing %g outlier trials with the following interbeat intervals (IBI): ', length(outliers))
    fprintf('   - IBI = %g (ms) \n', EEG.event(outliers).duration)
    IBI(outliers) = [];
    for iEv = 1:length(outliers)
        EEG.event(outliers(iEv)).type = 'outlier';
    end
else
    fprintf('No outlier trial detected. \n')
end

% Visualize IBI distribution and lower 95% percentile that will be used as
% epoch length
if params.vis
    figure('color','w'); histfit(IBI); hold on
    plot([prctile(IBI,5) prctile(IBI,5)],ylim,'--r','linewidth',2)
    % plot([prctile(IBI,97.5) prctile(IBI,97.5)],ylim,'--r','linewidth',2)
    title('Interbeat intervals (IBI) after removal of outlier trials'); xlabel('time (ms)')
    legend('','','lower 95% percentile (epoch size)')
end

% Epoch (no baseline removal as it can bias results because of
% pre R-peak activity)
HEP = pop_epoch(EEG,{},[-.3 prctile(IBI,5)/1000],'epochinfo','yes');
% warning('%g trials (i.e., 10%%) with an interbeat interval (IBI) <500 ms were removed. See Candia-Rivera et al. (2021) and Park and Blanke (2019) for more detail) . \n', length(shortTrials),quantile(IBI,.1));
% HEP = pop_rejepoch(HEP, shortTrials, 0);  % Trial numbers have changed so this is incorrect

% Remove bad epochs, run ICA, and remove bad components
if params.clean_eeg
    [HEP, params] = clean_eeg(HEP,params);
end

% Remove epochs containing more than 1 R-peak
warning('Removing remaining epochs containing more than 1 R-peak (i.e. remaining interbeat intervals that are too short)')
idx = false(length(HEP.epoch),1);
for iEv = 1:length(HEP.epoch)
    if length(HEP.epoch(iEv).eventtype)>1
        idx(iEv,1) = true;
    end
end
HEP = pop_rejepoch(HEP, find(idx), 0);
% HEP.event = [];
% HEP = eeg_checkset(HEP);


%% Plot Heartbeat-evoked potentials (HEP) and oscillations (HEO)

% Known effects are 200-600 ms after R-peak, with peak effect at 300-450 ms
% over frontocentral electrodes (Fz, Cz, Pz), most likely in the alpha
% band for ERSP. 

if params.vis

    figure
    pop_plottopo(HEP, 1:HEP.nbchan, 'Heartbeat-evoked potentials (HEP)', 0, 'ydir',1);
    set(findall(gcf,'type','axes'),'fontSize',11,'fontweight','bold');

    % HEP all elecs + scalp topographies of window of interest
    figure
    subplot(2,1,1)
    pop_timtopo(HEP, [HEP.times(1) HEP.times(end)], [300 400], ...
        'Heartbeat-evoked potentials (HEP) - all electrodes','verbose','off');

    % ERPimage of known peak electrode to see change over time
    elecName = 'Fz';
    elecNum = find(strcmpi({HEP.chanlocs.labels}, elecName));
    if isempty(elecNum)
        % if Fz is not present, try Cz
        elecName = 'Fz';
        elecNum = find(strcmpi({HEP.chanlocs.labels}, elecName));
        if isempty(elecNum)
            % if Cz is not present either, take 1st channel
            elecNum = 1;
            elecName = HEP.chanlocs(elecNum).labels;
        end
    end
    subplot(2,1,2)
    pop_erpimage(HEP,1, elecNum,[],sprintf('Heartbeat-evoked potentials (HEP) over time for channel %s',elecName), ...
        10,1,{'R-peak'},[],'','yerplabel','\muV','erp','on','cbar','on' );
    colormap("parula") % parula hot bone sky  summer winter
    set(findall(gcf,'type','axes'),'fontSize',11,'fontweight','bold');

    % % Headplot
    % pop_headplot(HEP, 1, 0, 'HEP', [1  1], 'setup', ...
    %     { fullfile(HEP.filepath,'sample_data', sprintf('%s_HEP.spl', HEP.filename(1:end-4))), ...
    %     'meshfile','mheadnew.mat','transform',[-1.136 7.7523 11.4527 -0.027117 0.015531 -1.5455 0.91234 0.93161 0.80698] });
    % colormap("parula")

    % Event-related spectral perturbations (ERSP) and inter-trial coherence (ITC) (non-corrected)
    figure('color','w');
    pop_newtimef(HEP,1,elecNum,[HEP.times(1) HEP.times(end)],[3 0.8], ...
        'topovec',1,'elocs',HEP.chanlocs,'chaninfo',HEP.chaninfo,'freqs',[5 30], ...
        'baseline',0,'plotphase','on','padratio',2,'caption', ...
        sprintf('%s (non-corrected)',elecName));
    colormap("parula") % "parula" "hot" "bone" "sky" "summer" "winter"
    
    % ERSP and ITC (bootstrap + FDR-corrected)
    figure('color','w');
    pop_newtimef(HEP,1,elecNum,[HEP.times(1) HEP.times(end)],[3 0.8],'topovec',1,'elocs', ...
        HEP.chanlocs,'chaninfo',HEP.chaninfo, 'freqs',[2 30], ...
        'baseline',0,'plotphase','on','padratio',2, ...
        'alpha',0.05,'mcorrect','fdr','naccu',1000,'caption', ...
        sprintf('%s (FDR-corrected)',elecName));
    colormap("parula") % parula hot bone sky summer winter

end

% Save
if params.save
    newname = sprintf('%s_HEP.set', HEP.filename(1:end-4));
    pop_saveset(HEP,'filename',newname,'filepath',HEP.filepath); % FIXME: add output
end


