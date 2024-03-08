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

function HEP = run_HEP(EEG, CARDIO, params, Rpeaks)

% Remove boundary events
idx = strcmpi({EEG.event.type}, 'boundary');
if any(idx)
    EEG.event(idx) = [];
end

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

% Add back heart channel (mainly for plotting to check if R-peaks events 
% align correctly with ECG signal)
if isfield(params,'keep_heart') && params.keep_heart
    % if EEG.srate ~= CARDIO.srate
    %     CARDIO = pop_resample(CARDIO,EEG.srate);
    % end
    CARDIO = pop_eegfiltnew(CARDIO,'locutoff',1,'hicutoff',20);    
    for iChan = 1:CARDIO.nbchan
        CARDIO.data(iChan,:) = rescale(CARDIO.data(iChan,:), -100, 100);
    end
    if EEG.pnts ~= CARDIO.pnts
        EEG.data(end+1:end+CARDIO.nbchan,:) = CARDIO.data(:,1:end-1); % for PPG
    else
        EEG.data(end+1:end+CARDIO.nbchan,:) = CARDIO.data;
    end
    EEG.nbchan = EEG.nbchan + CARDIO.nbchan;
    for iChan = 1:CARDIO.nbchan
        EEG.chanlocs(end+1).labels = params.heart_channels{iChan};
    end
    EEG = eeg_checkset(EEG);
end

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

% Remove epochs with IBI<550 ms (see Candia-Rivera et al. 2021; Park & Blanke 2019)
shortTrials = find(IBI<550);
if ~isempty(shortTrials)
    warning('Removing %g trials with interbeat intervals (IBI) < 550 ms', length(shortTrials))
    IBI(shortTrials) = [];
    [EEG.event(shortTrials).type] = deal([]);
else
    fprintf('All trial have an interbeat interval (IBI) above 550 ms (minimum recommended). \n')
end

% Remove remaining outlier epochs
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
if params.vis_outputs
    figure('color','w'); histfit(IBI); hold on
    plot([prctile(IBI,5) prctile(IBI,5)],ylim,'--r','linewidth',2)
    % plot([prctile(IBI,97.5) prctile(IBI,97.5)],ylim,'--r','linewidth',2)
    title('Interbeat intervals (IBI) after removal of outliers'); 
    xlabel('Time (ms)'); ylabel('Number of IBIs')
    legend('','','lower 95% percentile (epoch size)')
    try icadefs; set(gcf, 'color', BACKCOLOR); catch; end     % eeglab background color
    set(gcf,'Name','Inter-beat intervals (IBI) distribution','NumberTitle','Off','Toolbar','none','Menu','none')
    set(findall(gcf,'type','axes'),'fontSize',11,'fontweight','bold');
end

% Epoch (no baseline removal as it can bias results because of
% pre R-peak activity)
HEP = pop_epoch(EEG,{},[-.3 prctile(IBI,5)/1000],'epochinfo','yes');
% warning('%g trials (i.e., 10%%) with an interbeat interval (IBI) <500 ms were removed. See Candia-Rivera et al. (2021) and Park and Blanke (2019) for more detail) . \n', length(shortTrials),quantile(IBI,.1));
% HEP = pop_rejepoch(HEP, shortTrials, 0);  % Trial numbers have changed so this is incorrect

% Remove bad epochs, run ICA, and remove bad components
if params.clean_eeg
    [HEP, params] = clean_eeg(HEP,params);

    % Preprocessing outputs
    EEG.brainbeats.preprocessing.removed_eeg_trials = params.removed_eeg_trials;
    EEG.brainbeats.preprocessing.removed_eeg_components = params.removed_eeg_components;
end

% Only keep 1st marker when several are present in sam eepoch
count = 0;
for iEv = 1:length(HEP.epoch)
    if length(HEP.epoch(iEv).eventtype)>1
        HEP.epoch(iEv).event = HEP.epoch(iEv).event(1);
        HEP.epoch(iEv).eventlatency = HEP.epoch(iEv).eventlatency(1);
        HEP.epoch(iEv).eventduration = HEP.epoch(iEv).eventduration(1);
        HEP.epoch(iEv).eventtype = HEP.epoch(iEv).eventtype(1);
        count = count+1;
    end
end
warning('%g epochs contained more than 1 R-peak and were adjusted by preserving only the first event of the epochs',count)
HEP.epoch = rmfield(HEP.epoch,'eventurevent');
HEP = eeg_checkset(HEP);


%% Plot Heartbeat-evoked potentials (HEP) and oscillations (HEO)

% Known effects are 200-600 ms after R-peak, with peak effect at 300-450 ms
% over frontocentral electrodes (Fz, Cz, Pz), most likely in the alpha
% band for ERSP. 

if params.vis_outputs

    if params.vis_cleaning
        pop_eegplot(HEP,1,1,1);
        set(gcf,'Name','Final output','NumberTitle','Off','Toolbar','none','Menu','none');  % remove toolbobar and menu
    end

    % Show mean HEP for each electrodes and allows clicking on them to see
    % more closely
    % figure
    % pop_plottopo(HEP, 1:HEP.nbchan, 'Heartbeat-evoked potentials (HEP)', 0,...
    %     'ydir',1);
    options = { 'frames' HEP.pnts 'limits' [HEP.xmin HEP.xmax 0 0]*1000 ...
        'title' 'Heartbeat-evoked potentials (HEP)' 'chans' 1:HEP.nbchan ...
        'chanlocs' HEP.chanlocs 'ydir' 1 'legend' {'uV' 'Time (ms)'}};
    figure
    % plottopo( HEP.data, options{:} );           % single trials (long)
    % plottopo_mod( trimmean(HEP.data,20,3), options{:} );   % modified version with y label
    plottopo( trimmean(HEP.data,20,3), options{:} );   % average across trials
    try icadefs; set(gcf, 'color', BACKCOLOR); catch; end  % eeglab color

    set(findall(gcf,'type','axes'),'fontSize',11,'fontweight','bold');
    set(gcf,'Toolbar','none','Menu','none');  % remove toolbobar and menu
    % set(gcf,'Name','Mean HEP for EEG each channel','NumberTitle','Off')  % name

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
    set(gcf,'Toolbar','none','Menu','none');  % remove toolbobar and menu
    set(gcf,'Name','HEP','NumberTitle','Off')  % name

    % % Headplot
    % pop_headplot(HEP, 1, 0, 'HEP', [1  1], 'setup', ...
    %     { fullfile(HEP.filepath,'sample_data', sprintf('%s_HEP.spl', HEP.filename(1:end-4))), ...
    %     'meshfile','mheadnew.mat','transform',[-1.136 7.7523 11.4527 -0.027117 0.015531 -1.5455 0.91234 0.93161 0.80698] });
    % colormap("parula")
    % set(gcf,'Toolbar','none','Menu','none');  % remove toolbobar and menu
    % set(gcf,'Name','HEP','NumberTitle','Off')  % name

    % ERSP and ITC (uncorrected)
    % figure('color','w');
    % pop_newtimef(HEP,1,elecNum,[HEP.times(1) HEP.times(end)],[3 0.8], ...
    %     'topovec',1,'elocs',HEP.chanlocs,'chaninfo',HEP.chaninfo,'freqs',[7 25], ...
    %     'baseline',0,'plotphase','on','padratio',2,'caption', ...
    %     sprintf('%s (uncorrected)',elecName));
    % colormap("parula") % "parula" "hot" "bone" "sky" "summer" "winter"
    % set(gcf,'Toolbar','none','Menu','none');  % remove toolbobar and menu
    % set(gcf,'Name','Heartbeat-evoked oscillations (HEO) and ITC','NumberTitle','Off')  % name

    % ERSP and ITC (bootstrap + FDR-corrected)
    figure('color','w');
    pop_newtimef(HEP,1,elecNum,[HEP.times(1) HEP.times(end)],[3 0.8],...
        'topovec',1,'elocs', HEP.chanlocs,'chaninfo',HEP.chaninfo, 'freqs',[7 25], ...
        'baseline',0,'plotphase','on','padratio',2, ...
        'alpha',0.05,'mcorrect','fdr','naccu',1000,...
        'caption',sprintf('Channel %s (p=0.05; FDR-corrected)',elecName));
    colormap("parula") % parula hot bone sky summer winter
    set(gcf,'Toolbar','none','Menu','none');  % remove toolbobar and menu
    set(gcf,'Name','Heartbeat-evoked oscillations (HEO) and ITC','NumberTitle','Off')  % name

end

% Save
if params.save
    newname = sprintf('%s_HEP.set', HEP.filename(1:end-4));
    pop_saveset(HEP,'filename',newname,'filepath',HEP.filepath);
end


