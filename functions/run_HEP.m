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


%% Plot HEP
% Known effects are 200-600 ms after R-peak, with peak effect at 300-450 ms
% over frontocentral electrodes (Fz, Cz, Pz)
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
end

%% Heartbeats-evoked oscillations (HEO; i.e. ERSP)

% Apply reflection to control for edge effects by adding a backward-version
% of the signal before and after (developed by Mike X Cohen).
nFrames = size(HEP.data,2); % epoch length before applying reflection
ref_length = nFrames*3;     % epoch length after reflection (will be 3 times longer)
ori_lats = [HEP.times(1) HEP.times(end)];   % latencies of original epoch
TMPDATA = nan(HEP.nbchan,ref_length,size(HEP.data,3));
for iChan = 1:HEP.nbchan
    for iEpoch = 1:size(HEP.data,3)
        sig = squeeze(HEP.data(iChan,:,iEpoch));
        TMPDATA(iChan,:,iEpoch) = [ sig(end:-1:1) sig sig(end:-1:1) ];
    end
end
fres = HEP.times(2) - HEP.times(1);
newlat1 = HEP.times(1) - nFrames*fres;
newlat2 = HEP.times(end) + nFrames*fres;
TMPTIMES = newlat1:fres:newlat2;
% check, epochs should be 3 times longer now
if size(TMPDATA,2) ~= size(HEP.data,2)*3
    error('Reflection failed to control for edge effects in run_HEP.m')
end
pop_eegplot(HEP,1,1,1);
HEP2 = HEP;
HEP2.data = TMPDATA;
HEP2.times = TMPTIMES;
HEP2.pnts = ref_length;
for iEv = 1:length(HEP2.epoch)
    HEP2.epoch(iEv).eventlatency = 0;  % here should all be 0
    %     if iEv==1
    %         HEP2.event(iEv).latency = HEP2.epoch(iEv).eventlatency; % but here updates for the whole length of the file
    %     else
    %         HEP2.event(iEv).latency = HEP2.event(iEv-1).latency + HEP2.event(iEv-1).duration; % but here updates for the whole length of the file
    %     end
    %     HEP2.event(iEv).type = HEP2.epoch(iEv).eventtype;
    %     HEP2.event(iEv).urevent = iEv;
    %     HEP2.event(iEv).duration = HEP2.epoch(iEv).eventduration;
    %     HEP2.event(iEv).epoch = iEv;
end
HEP2 = eeg_checkset(HEP2); % update all fields
% pop_eegplot(HEP2,1,1,1);  % it looks like htere are extra events, but
% it's not the case when running analyses (visible in HEP.event or
% HEP.epoch)
HEP = HEP2; clear HEP2

% ERSP and ITC (non-corrected)
% figure('color','w');
% pop_newtimef(HEP,1,elecNum,[],[3 0.8],'topovec',1,'elocs', ...
%     HEP.chanlocs,'chaninfo',HEP.chaninfo, 'freqs',[7 30], ...
%     'baseline',0,'plotphase','on','padratio',2,'caption', ...
%     sprintf('%s (non-corrected)',elecName));
% colormap("parula") % parula hot bone sky summer winter

% ERSP and ITC (bootstrap + FDR-corrected)
% figure('color','w');
% pop_newtimef(HEP,1,elecNum,[],[3 0.8],'topovec',1,'elocs', ...
%     HEP.chanlocs,'chaninfo',HEP.chaninfo, 'freqs',[7 30], ...
%     'baseline',0,'plotphase','on','padratio',2, ...
%     'alpha',0.05,'mcorrect','fdr','caption', ...
%     sprintf('%s (FDR-corrected)',elecName));
% pop_newtimef(HEP2,1,elecNum,[] ,[3 0.8],'topovec',1,'elocs', ...
%     HEP.chanlocs,'chaninfo',HEP.chaninfo, 'freqs',[7 30], ...
%     'baseline',0,'plotphase','on','padratio',2, ...
%     'alpha',0.05,'mcorrect','fdr','caption', ...
%     sprintf('%s (FDR-corrected)',elecName));
% colormap("parula") % parula hot bone sky summer winter

% newtimef( tmpsig(:, :), length(pointrange), [tlimits(1) tlimits(2)], EEG.srate, cycles ,'topovec',1,'elocs',struct('labels',{'AF7','TP10'},'type',{'',''},'theta',{-38.6506,118.6299},'radius',{0.53821,0.63914},'X',{68.5722,-47.0353},'Y',{54.8397,-86.1618},'Z',{-10.59,-45.869},'sph_theta',{38.6506,-118.6299},'sph_phi',{-6.8772,-25.0452},'sph_radius',{88.4403,108.3519},'urchan',{2,4},'ref',{'average','average'}),'chaninfo',struct('plotrad',{[]},'shrink',{[]},'nosedir',{'+X'},'nodatchans',{struct([])},'icachansind',{[1 2] },'filename',{'C:\\\\Users\\\\IONSLAB\\\\Documents\\\\MATLAB\\\\eeglab\\\\plugins\\\\dipfit\\\\standard_BEM\\\\elec\\\\standard_1005.elc'},'removedchans',{struct('labels',{'AF8','AUX','ECG','ECG','TP9'},'type',{'','','','',''},'theta',{38.6688,[],[],[],-118.5141},'radius',{0.53819,[],[],[],0.63961},'X',{69.6568,[],[],[],-46.5147},'Y',{-55.7433,[],[],[],85.6192},'Z',{-10.755,[],[],[],-45.707},'sph_theta',{-38.6688,[],[],[],118.5141},'sph_phi',{-6.8739,[],[],[],-25.1306},'sph_radius',{89.8613,[],[],[],107.6262},'urchan',{3,5,[],[],1},'ref',{'','','','average','average'})}),'freqs',[7 30] ,'baseline',0,'plotphase','on','padratio',2,'alpha',0.05,'mcorrect','fdr','caption','AF7 (FDR-corrected)');'
% newtimef(1x46374, 262, [-300.7812  718.7500], 256, [3 0.8],...
% 'topovec',1,'elocs',HEP.chanlocs
frames = size(HEP.data,2);
winsize = max(pow2(nextpow2(frames)-3),4);
tlimits = [HEP.times(1) HEP.times(end)];    % time limits of the epochs, not 
                                            % a sub-window to extract  
flimits = [7 30];
detrend = 'off';        % linearly detrend each epoch {default: 'off'}
itctype = 'phasecoher';
cycles = [3 0.8];
verbose = 'on';
% ntimes = [];
padratio = 2;           % FFT-length/winframes (2^k) {default: 2}
                        % Multiplies the number of output frequencies by dividing
                        % their spacing (standard FFT padding). When cycles~=0,
freqscale = 'linear';   % ['log'|'linear'] frequency scale. {default: 'linear'}
                        % Note that for obtaining 'log' spaced freqs using FFT,
                        % closest correspondent frequencies in the 'linear' space
                        % are returned.
% srate = HEP.srate;
scale = 'log';          % visualize power in log (default; dB) or absolute ('abs')
                        % scale. 
for iChan = 1:HEP.nbchan
    % Compute time-frequency
    [tf, freqs, timesout, itcvals] = timefreq(squeeze(HEP.data(iChan,:,:)), ...
        HEP.srate,'winsize', winsize,'tlimits', tlimits, 'detrend', detrend, ...
        'itctype', itctype, 'cycles', cycles, 'verbose', verbose, ...
        'padratio', padratio, 'freqs', flimits, 'freqscale', freqscale);
    
    % multiply by 2 account for negative frequencies and counteract the 
    % reduction by a factor 0.375 that occurs as a result of cosine (Hann) 
    % tapering.
    if cycles(1) == 0
        tf = 2/0.375*tf/winsize; % normalization, divide by winsize
        P  = tf.*conj(tf); % power
    else
        P  = tf.*conj(tf); % power for wavelets
    end

    % Chop off reflection periods
    idx = timesout>=ori_lats(1) & timesout<=ori_lats(2);
    timesout = timesout(idx);
    tf = tf(:,idx,:);
    itcvals = itcvals(:,idx);

    % Remove baseline
    bsln_lat = [-300 -1];             % baseline end-time (in ms). NaN --> no baseline is used. 
                            % A [min max] range may also be entered
                            % You may also enter one row per region for baseline
                            % e.g. [0 100; 300 400] considers the window 0 to 100 ms and
                            % 300 to 400 ms. This parameter validly defines all baseline types 
                            % below. {default: 0 -> all negative time values}. 
    powbase = NaN;          % baseline spectrum to log-subtract {default|NaN -> calculate from data}
    basenorm = 'on';       % 'on' normalize baseline in the power spectral
                            % average; 'off' (default) divide by the average power across 
                            % trials at each frequency (gain model)
    trialbase = 'off';      % perform baseline normalization or division above 
                            % in single trial instead of the trial average
                            % (default: off).
    if strcmpi(scale, 'log') && ~any(isnan(powbase))
        powbase = 10.^(powbase/10);
    end
    % P = newtimeftrialbaseln(P,timesout,'baseline',bsln_lat,'basenorm',basenorm,...
    %     'trialbase',trialbase);
    % [P, baseln, mbase] = newtimefbaseln(P,timesout,'baseline',bsln_lat, ...
    %     'basenorm',basenorm,'verbose',verbose,'powbase',powbase, ...
    %     'trialbase',trialbase,'singletrials','on');

    if ~any(isnan(bsln_lat)) && ~any(isnan(bsln_pwr))
        P_mu = squeeze(trimmean(P,20,3));       % trimmed mean power across trials
        idx = timesout>=bsln_lat(1) & timesout<=bsln_lat(2);
        bsln_pwr = trimmean(P_mu(:,idx),20,2);  % trimmed mean across frames
        bsln_sd = std(P_mu(:,idx),[],2);               % SD across frames

         % baseline removal set to off, divide by mean power
        if strcmpi(basenorm, 'off')
            P = bsxfun(@rdivide, P, bsln_pwr);
        
        % baseline normalization
        else 
            P = bsxfun(@rdivide, bsxfun(@minus,P, bsln_pwr), bsln_sd);
        end
    % elseif


    end

    % convert to log if necessary
    if params.norm
        P = 10 * log10(P);
        % if ~isempty(Pboot)
        %     Pboot = 10 * log10(Pboot);
        % end
    end

    % auto scalling
    % erspmax = [max(max(abs(P)))]/2;
    if strcmpi(params.scale,'abs') && strcmpi(params.baseline,'off')
        erspmax = [max(max(abs(P)))];
        if erspmax > 1
            erspmax = [1-(erspmax-1) erspmax];
        else 
            erspmax = [erspmax 1+(1-erspmax)];
        end
    end

    % plot
    if params.vis
        if ndims(P) == 3
            P = squeeze(P(2,:,:,:));
            itcvals = squeeze(itcvals(2,:,:,:));
            mbase = squeeze(mbase(2,:));
            ERP = mean(squeeze(data(1,:,:)),2);
        else
            ERP = mean(data,2);
        end
        if strcmpi(g.plottype, 'image')
            plottimef(P, itcvals, Pboot, Rboot, ERP, freqs, timesout, mbase, maskersp, maskitc, g);
        else
            plotallcurves(P, itcvals, Pboot, Rboot, ERP, freqs, timesout, mbase, g);
        end
    end

    % output
    if strcmpi(g.outputformat, 'old')
        itcvals = abs(itcvals); % convert coherence vector to magnitude
        if strcmpi(g.scale, 'log'), mbase = 10.^(mbase/10); end
    end
    if strcmpi(g.verbose, 'on')
        disp('Note: Add output variables to command line call in history to');
        disp('      retrieve results and use the tftopo function to replot them');
    end
    mbase = mbase';


end

% Save
if params.save
    newname = sprintf('%s_HEP.set', HEP.filename(1:end-4));
    pop_saveset(HEP,'filename',newname,'filepath',HEP.filepath); % FIXME: add output
end

%%%%%%%%%%%% HEO (ERSP) Same but with wider epochs %%%%%%%%%%%%%%
% HEO = pop_epoch(HEP,{},[-.3 .7],'epochinfo','yes');
% warning('Removing %g trials shorter than [-300 700] ms long for heartbeat-evoked oscillations (HEO) analysis. \n', length(shortTrials));
% HEO = pop_rejepoch(HEO, shortTrials, 0);
%
% if params.clean_eeg
%     params.clean_eeg_step = 1;
%     [HEO, params] = clean_eeg(HEO,params);
% end
%
% if params.vis
%     if sum(strcmpi({HEP.chanlocs.labels}, 'cz')) > 0
%         elec = find(strcmpi({HEP.chanlocs.labels}, 'cz'));
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

