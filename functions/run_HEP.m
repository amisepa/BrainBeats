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
% figure('color','w'); histogram(IBI); hold on; 
% histogram(diff(Rpeaks)/EEG.srate *1000); legend('from EEG events','from R peaks');

% Remove trials with IBI<550 ms. see Candia-Rivera et al. (2021)  and Park & Blanke (2019)
shortTrials = find(IBI<550);
if ~isempty(shortTrials)
    warning('Removing %g trials with interbeat intervals (IBI) < 550 ms', length(shortTrials))
    IBI(shortTrials) = [];
    EEG.event(shortTrials).type = [];
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

if params.vis
    figure('color','w'); histfit(IBI); hold on
    plot([prctile(IBI,5) prctile(IBI,5)],ylim,'--r','linewidth',2)
    % plot([prctile(IBI,97.5) prctile(IBI,97.5)],ylim,'--r','linewidth',2)
    title('Interbeat intervals (IBI) after removal of outlier trials'); xlabel('time (ms)')
    legend('','','lower 95% percentile (epoch size)')
end

% Epoch
HEP = pop_epoch(EEG,{},[-.3 prctile(IBI,5)/1000],'epochinfo','yes'); 
% warning('%g trials (i.e., 10%%) with an interbeat interval (IBI) <500 ms were removed. See Candia-Rivera et al. (2021) and Park and Blanke (2019) for more detail) . \n', length(shortTrials),quantile(IBI,.1));
% HEP = pop_rejepoch(HEP, shortTrials, 0);  % Trial numbers have changed so this is incorrect

% Remove bad trials and eye/muscle components (ICLabel)
if params.clean_eeg
    [HEP, params] = clean_eeg(HEP,params);
end

if params.vis

    figure; pop_plottopo(HEP, 1:HEP.nbchan, 'Heartbeat-evoked potentials (HEP)', 0, 'ydir',1);

    try
        elecName = 'Fpz';
        elecNum = find(strcmpi({EEG.chanlocs.labels}, elecName));
    catch
        elecNum = 1;
        elecName = EEG.chanlocs(elecNum).labels;
    end

    figure; 
    subplot(2,1,1)
    pop_timtopo(HEP, [HEP.times(1) HEP.times(end)], [250 350 450], 'Heartbeat-evoked potentials (HEP) - all electrodes');

    subplot(2,1,2)
    pop_erpimage(HEP,1, elecNum,[], sprintf('Heartbeat-evoked potentials (HEP) - %s',elecName),10,1,{},[],'','yerplabel','\muV','erp','on','cbar','on','topo', { 1 EEG.chanlocs EEG.chaninfo } );
    colormap("parula") % parula hot bone sky  summer winter

    % pop_headplot(HEP, 1, 0, 'HEP', [1  1], 'setup', ...
    %     { fullfile(HEP.filepath,'sample_data', sprintf('%s_HEP.spl', HEP.filename(1:end-4))), ...
    %     'meshfile','mheadnew.mat','transform',[-1.136 7.7523 11.4527 -0.027117 0.015531 -1.5455 0.91234 0.93161 0.80698] });
    % colormap("parula") 

    % ERSP and ITC (non-corrected)
    figure('color','w'); 
    pop_newtimef(HEP,1,elecNum,[],[3 0.8],'topovec',1,'elocs', ...
        HEP.chanlocs,'chaninfo',HEP.chaninfo, 'freqs',[7 30], ...
        'baseline',0,'plotphase','on','padratio',2,'caption', ...
        sprintf('%s (non-corrected)',elecName));
    colormap("parula") % parula hot bone sky summer winter

    % ERSP and ITC (bootstrap + FDR-corrected)
    figure('color','w'); 
    pop_newtimef(HEP,1,elecNum,[],[3 0.8],'topovec',1,'elocs', ...
        HEP.chanlocs,'chaninfo',HEP.chaninfo, 'freqs',[7 30], ...
        'baseline',0,'plotphase','on','padratio',2, ...
        'alpha',0.05,'mcorrect','fdr','caption', ...
        sprintf('%s (FDR-corrected)',elecName));
    colormap("parula") % parula hot bone sky summer winter
    
end

% Save
if params.save
    newname = sprintf('%s_HEP.set', HEP.filename(1:end-4));
    pop_saveset(HEP,'filename',newname,'filepath',HEP.filepath); % FIXME: add output
end

%%%%%%%%%%%% HEO (ERSP) Same but with wider epochs %%%%%%%%%%%%%%
% HEO = pop_epoch(EEG,{},[-.3 .7],'epochinfo','yes');
% warning('Removing %g trials shorter than [-300 700] ms long for heartbeat-evoked oscillations (HEO) analysis. \n', length(shortTrials));
% HEO = pop_rejepoch(HEO, shortTrials, 0);
% 
% if params.clean_eeg
%     params.clean_eeg_step = 1;
%     [HEO, params] = clean_eeg(HEO,params);
% end
% 
% if params.vis
%     if sum(strcmpi({EEG.chanlocs.labels}, 'cz')) > 0
%         elec = find(strcmpi({EEG.chanlocs.labels}, 'cz'));
%     else
%         elec = 1;
%     end
% 
%     % ERSP
    % figure; pop_newtimef(HEO, 1, elec, [], [3 0.8], 'topovec', 1, ...
    %     'elocs', HEO.chanlocs, 'chaninfo', HEO.chaninfo, ...
    %     'caption', HEO.chanlocs(elec).labels,'baseline',[-300 -200], ...
    %     'freqs', [4 13], 'plotphase','off','padratio',1);
    % colormap("parula") % parula hot bone sky  summer winter
% end
% 
% % Save
% if params.hep_save
%     newname = sprintf('%s_HEO.set', HEO.filename(1:end-4));
%     pop_saveset(HEO,'filename',newname,'filepath',HEO.filepath);% FIXME: add output
% end

