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

function EEG = run_HEP(EEG, params, Rpeaks)

% Save HEP/HEO files
if ~isfield(params,'hep_save') % not available from GUI yet
    params.hep_save = true;
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

% Remove inter-beat-intervals (IBI) shorter than 700 ms (should be at least 
% 500 ms post R peak but extend window for ERSP). 
% Candia-Rivera, Catrambone, & Valenza (2021). The role of EEG reference 
% in the assessment of functional brain–heart interplay: From 
% methodology to user guidelines. Journal of Neuroscience Methods.
IBI = diff(Rpeaks)/EEG.srate *1000;
shortTrials = find(IBI<700);
if params.vis
    figure; histogram(IBI); 
    title('Gaps between heartbeats (minimum should be 500 ms)'); xlabel('millisecond'); ylabel('Gaps');
end
if length(shortTrials)/length(IBI) < .05
    IBI(shortTrials) = []; % remove them if <5% is very short to get min trial possible below
end

% Minimum ICI = 500 ms (but more is better for HEO/ERSP)
% if quantile(gaps,.33) < 500
%     warning("At least 33% between heartbeats are %g% long. HEP with epochs < 500 ms are strongly discouraged.")
%     warning("see Candia-Rivera et al. (2021) and Park and Blanke (2019) for more detail.")
%     EEG = pop_epoch(EEG,{},[-0.05 .5],'epochinfo','yes');
% elseif quantile(gaps,.1) >= 700 && quantile(gaps,.1) < 800
%     EEG = pop_epoch(EEG,{},[-.05 quantile(gaps,.1)/1000],'epochinfo','yes');
% else
%     EEG = pop_epoch(EEG,{},[-.05 min(gaps)/1000],'epochinfo','yes');
% end

% Epoch
HEP = pop_epoch(EEG,{},[-.05 .7],'epochinfo','yes');
warning('Removing %g trials shorter than 550 ms long, for heartbeat-evoked potential (HEP) analysis. \n', length(shortTrials));
HEP = pop_rejepoch(HEP, shortTrials, 0);

if params.clean_eeg
    [HEP, params] = clean_eeg(HEP,params);
end

if params.vis
    if sum(strcmpi({EEG.chanlocs.labels}, 'cz')) > 0
        elec = find(strcmpi({EEG.chanlocs.labels}, 'cz'));
    else
        elec = 1;
    end

    figure; 
    subplot(2,1,1)
    pop_timtopo(HEP, [HEP.times(1) HEP.times(end)], [250 350 450], 'Heartbeat-evoked potentials (HEP) - all electrodes');
    subplot(2,1,2)
    pop_erpimage(HEP,1, elec,[],'Heartbeat-evoked potentials (HEP) - Fcz',10,1,{},[],'','yerplabel','\muV','erp','on','cbar','on','topo', { 1 EEG.chanlocs EEG.chaninfo } );
    colormap("parula") % parula hot bone sky  summer winter
    figure; pop_plottopo(HEP, 1:HEP.nbchan, 'Heartbeat-evoked potentials (HEP)', 0, 'ydir',1);
    % pop_headplot(EEG, 1, 0, 'HEP', [1  1], 'setup',{fullfile(dataDir,'sample_data','sample_data2_HEP.spl'),'meshfile','mheadnew.mat','transform',[-1.136 7.7523 11.4527 -0.027117 0.015531 -1.5455 0.91234 0.93161 0.80698] });
    % colormap("parula") 
end

% Save
if params.hep_save
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
%     figure; pop_newtimef(HEO, 1, elec, [], [3 0.8], 'topovec', 1, ...
%         'elocs', HEO.chanlocs, 'chaninfo', HEO.chaninfo, ...
%         'caption', HEO.chanlocs(elec).labels,'baseline',[-300 -200], ...
%         'freqs', [4 13], 'plotphase','off','padratio',1);
%     colormap("parula") % parula hot bone sky  summer winter
% end
% 
% % Save
% if params.hep_save
%     newname = sprintf('%s_HEO.set', HEO.filename(1:end-4));
%     pop_saveset(HEO,'filename',newname,'filepath',HEO.filepath);% FIXME: add output
% end

