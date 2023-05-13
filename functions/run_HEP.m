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

% Remove trials shorter than 500 ms post R peak
% Candia-Rivera, Catrambone, & Valenza (2021). The role of EEG reference 
% in the assessment of functional brain–heart interplay: From 
% methodology to user guidelines. Journal of Neuroscience Methods.
gaps = diff(Rpeaks)/EEG.srate *1000;
shortTrials = find(gaps<550);
if params.vis
    figure; histogram(gaps); 
    title('Gaps between heartbeats (minimum should be 500 ms)'); xlabel('millisecond'); ylabel('Gaps');
end
if length(shortTrials)/length(gaps) < .05
    gaps(shortTrials) = []; % remove them if <5% is very short to get min trial possible below
end

% Epoch (-50 to + 500 ms is minimum)
% Warn if a third of trials are shorter than 500 ms long
% If 10% smallest are above 700 ms, epoch at 700 for time-frequency (200 ms more than 500 ms to allow 1 cycle when 5 Hz is lowest freq)
% if quantile(gaps,.33) < 500
%     warning("At least 33% between heartbeats are %g% long. HEP with epochs < 500 ms are strongly discouraged.")
%     warning("see Candia-Rivera et al. (2021) and Park and Blanke (2019) for more detail.")
%     EEG = pop_epoch(EEG,{},[-0.05 .5],'epochinfo','yes');
% elseif quantile(gaps,.1) >= 700 && quantile(gaps,.1) < 800
%     EEG = pop_epoch(EEG,{},[-.05 quantile(gaps,.1)/1000],'epochinfo','yes');
% else
%     EEG = pop_epoch(EEG,{},[-.05 min(gaps)/1000],'epochinfo','yes');
% end

%%%%%%%%%%%% HEP (ERP) %%%%%%%%%%%%%%%%
HEP = pop_epoch(EEG,{},[-.05 .5],'epochinfo','yes'); % FIXME: 500 ok for ERP, 700 better for ERSP
warning('Removing %g trials shorter than 550 ms long, for heartbeat-evoked potential (HEP) analysis. \n', length(shortTrials));
HEP = pop_rejepoch(HEP, shortTrials, 0);

if params.clean_eeg
    [HEP, params] = clean_eeg(HEP,params);
end

if params.vis
    if sum(strcmpi({EEG.chanlocs.labels}, 'fcz')) > 0
        elec = find(strcmpi({EEG.chanlocs.labels}, 'fcz'));
    else
        elec = 1;
    end

    % ERP
    figure; pop_erpimage(HEP,1, elec,[],'HEP (Fcz)',10,1,{},[],'','yerplabel','\muV','erp','on','cbar','on','topo', { 1 EEG.chanlocs EEG.chaninfo } );
    figure; pop_timtopo(HEP, [HEP.times(1) HEP.times(end)], [250 350 450]);
    figure; pop_plottopo(HEP, 1:HEP.nbchan, 'HEP data', 0, 'ydir',1);
    colormap("parula") % parula hot bone sky  summer winter
end

% Save
newname = sprintf('%s_HEP.set', HEP.filename(1:end-4));
% pop_saveset(HEP,'filename',newname,'filepath',HEP.filepath); % FIXME: add output

%%%%%%%%%%%% FOR HEO (ERSP) Same but with wider epochs %%%%%%%%%%%%%%
HEO = pop_epoch(EEG,{},[-.3 .7],'epochinfo','yes'); % FIXME: 500 ok for ERP, 700 better for ERSP
warning('Removing %g trials shorter than [-300 700] ms long for heartbeat-evoked oscillations (HEO) analysis. \n', length(shortTrials));
HEO = pop_rejepoch(HEO, shortTrials, 0);

if params.clean_eeg
    params.clean_eeg_step = 1;
    [HEO, params] = clean_eeg(HEO,params);
end

if params.vis
    if sum(strcmpi({EEG.chanlocs.labels}, 'cz')) > 0
        elec = find(strcmpi({EEG.chanlocs.labels}, 'cz'));
    else
        elec = 1;
    end
    
    % ERSP
    figure; pop_newtimef(HEO, 1, elec, [], [3 0.8], 'topovec', 1, ...
        'elocs', HEO.chanlocs, 'chaninfo', HEO.chaninfo, ...
        'caption', HEO.chanlocs(elec).labels,'baseline',[-300 -200], ...
        'freqs', [4 13], 'plotphase','off','padratio',1);
    colormap("parula") % parula hot bone sky  summer winter
end

% Save
newname = sprintf('%s_HEO.set', HEO.filename(1:end-4));
% pop_saveset(HEO,'filename',newname,'filepath',HEO.filepath);% FIXME: add output


