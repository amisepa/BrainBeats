function plot_features(Features,params)

disp('Plotting features...')

% set params if not in input (useful when users want to replot things
% without having to re-enter all the params)
if nargin==1 || ~isfield(params,'hrv') || ~isfield(params,'eeg') 

% HRV
    if any(contains(fieldnames(Features), {'HRV'}))
        params.hrv = true;
    else
        params.hrv = false;
    end
    if any(contains(fieldnames(Features.HRV), {'frequency'}))
        params.hrv_frequency = true;
        params.norm = true;     % use default (just for units)
    else
        params.hrv_frequency = false;
    end
    
    % EEG
    if any(contains(fieldnames(Features), {'EEG'}))
        params.eeg = true;
    else
        params.eeg = false;
    end
    if any(contains(fieldnames(Features.EEG), {'frequency'}))
        params.eeg_frequency = true;
        params.norm = true;     % use default (just for units)
    else
        params.eeg_frequency = false;
    end
    if any(contains(fieldnames(Features.EEG), {'nonlinear'}))
        params.eeg_nonlinear = true;
    else
        params.eeg_nonlinear = false;
    end
end

% EEG channel locations
if params.eeg && ~isfield(params,'chanlocs')
    errordlg('Sorry, you need to load your EEG channel locations into params.chanlocs to plot EEG features (see tutorial).')
end

% Pull features data
if params.hrv
    HRV = Features.HRV;
end
if params.eeg
    EEG = Features.EEG;
end

% abort if empty
if ~params.hrv && ~params.eeg
    fprintf('No features to plot. \n');
    return
end

%% PSD AND MFE PLOT


% PSD - HRV
if params.hrv && params.hrv_frequency

    figure('color','w');

    % Subplot mode if EEG PSD was also extracted
    if params.eeg && params.eeg_frequency
        % subplot(2,1,1); 
        nexttile([2 3])
    end
    
    hold on
    pwr = HRV.frequency.pwr; % power already averaged across windows
    freqs = HRV.frequency.pwr_freqs;
    bandNames = {'ULF' 'VLF' 'LF' 'HF'};
    bands = HRV.frequency.bands;
    baseval = min(pwr);

    idx = find(strcmp(bandNames,'ULF'));
    if isfield(HRV.frequency, 'ulf')
        x = freqs >= bands(idx,1) & freqs <= bands(idx,2);
        y = pwr(x);
        area(freqs(x),y,'BaseValue',baseval,'FaceColor',"#A2142F",'FaceAlpha',.7);
    else
        bandNames(idx) = [];
        bands(idx,:) = [];
    end

    idx = find(strcmp(bandNames,'VLF'));
    if isfield(HRV.frequency, 'vlf')
        x = freqs >= bands(idx,1) & freqs <= bands(idx,2);
        y = pwr(x);
        area(freqs(x),y,'BaseValue',baseval,'FaceColor',"#D95319",'FaceAlpha',.7)
    else
        bandNames(idx) = [];
        bands(idx,:) = [];
    end

    idx = find(strcmp(bandNames,'LF'));
    if isfield(HRV.frequency, 'lf')
        x = freqs >= bands(idx,1) & freqs <= bands(idx,2);
        y = pwr(x);
        area(freqs(x),y,'BaseValue',baseval,'FaceColor',"#EDB120",'FaceAlpha',.7)
    else
        bandNames(idx) = [];
        bands(idx,:) = [];
    end

    idx = find(strcmp(bandNames,'HF'));
    if isfield(HRV.frequency, 'hf')
        x = freqs >= bands(idx,1) & freqs <= bands(idx,2);
        y = pwr(x);
        area(freqs(x),y,'BaseValue',baseval,'FaceColor',"#0072BD",'FaceAlpha',.7)
    else
        bandNames(idx) = [];
        bands(idx,:) = [];
    end
    warning off
    legend(bandNames)
    warning on
    % xticklabels(unique(reshape(bands,1,[]))); xtickangle(45);
    axis tight; box on
    xlabel('Frequency (Hz)');
    if params.norm
        ylabel('Power (ms^2/Hz normalized)');
    else
        ylabel('Power (ms^2/Hz)');
    end
    title(sprintf('Power spectral density - HRV'))
    set(gcf,'Toolbar','none','Menu','none');  % remove toolbobar and menu
    set(gcf,'Name','Visualization of features','NumberTitle','Off')  % name

end

% % Multiscale fuzzy entropy (MFE) - HRV
% if params.hrv && params.hrv_nonlinear
%   nexttile      
%   hold on
%     mfe = HRV.nonlinear.MFE;
%     scales = HRV.nonlinear.MFE_scales;
%     area(scales,mfe,'FaceColor',"#A2142F",'FaceAlpha',.7);
%     title('Multiscale fuzzy entropy - HRV'); xlabel('Scale factors'); ylabel('Entropy')
%     axis tight; box on; grid on
% end

% PSD - EEG
if params.eeg && params.eeg_frequency

    if params.hrv && params.hrv_frequency    
        % subplot(2,1,2); 
        nexttile([2 3])
    else
        figure('color','w')
    end

    hold on
    pwr = trimmean(EEG.frequency.pwr_dB,20,1); % 20% trimmed mean across channels
    freqs = EEG.frequency.freqs(1,:);
    bands = [0 3; 3 7; 7 13; 13 30; 30 max(freqs)];
    baseval = min(pwr);

    % delta
    x = freqs >= bands(1,1) & freqs <= bands(1,2);
    y = pwr(x);
    area(freqs(x),y,'BaseValue',baseval,'FaceColor',"#A2142F",'FaceAlpha',.7)

    % theta
    x = freqs >= bands(2,1) & freqs <= bands(2,2);
    y = pwr(x);
    area(freqs(x),y,'BaseValue',baseval,'FaceColor',"#D95319",'FaceAlpha',.7)

    % alpha
    x = freqs >= bands(3,1) & freqs <= bands(3,2);
    y = pwr(x);
    area(freqs(x),y,'BaseValue',baseval,'FaceColor',"#EDB120",'FaceAlpha',.7)

    % beta
    x = freqs >= bands(4,1) & freqs <= bands(4,2);
    y = pwr(x);
    area(freqs(x),y,'BaseValue',baseval,'FaceColor',"#0072BD",'FaceAlpha',.7)

    % gamma
    x = freqs >= bands(5,1) & freqs <= bands(5,2);
    y = pwr(x);
    area(freqs(x),y,'BaseValue',baseval,'FaceColor',"#4DBEEE",'FaceAlpha',.7)

    % bandNames = {'Delta Theta Alpha Beta'};
    % xticks = freqs;
    % xticklabels(unique(reshape(bands,1,[])));
    % xtickangle(45); axis tight; box on
    title('Power spectral density - EEG'); 
    warning off
    legend({'delta' 'theta' 'alpha' 'beta' 'gamma'})
    warning on
    xlabel('Frequency (Hz)'); ylabel('Power (decibels)'); axis tight;

    set(gcf,'Toolbar','none','Menu','none');  % remove toolbobar and menu
    set(gcf,'Name','Visualization of features','NumberTitle','Off')  % name

end

% Multiscale fuzzy entropy (MFE) - EEG
% if params.eeg && params.eeg_nonlinear
% 
%     nexttile
%     hold on
%     scales = EEG.nonlinear.MFE_scales(1,:);
%     scaleBounds = EEG.nonlinear.MFE_scale_bounds(1,:);
%     % mfe = mean(EEG.nonlinear.MFE,1); % mean across channels
%     mfe = trimmean(EEG.nonlinear.MFE,20,1); % 20% trimmed mean across channels
%     % mfe = EEG.nonlinear.MFE(31,:); % Pz
%     area(scales,mfe(end:-1:1),'FaceColor',"#A2142F",'FaceAlpha',.7);
%     axis tight; box on; grid on
%     if ~isempty(scaleBounds)
%         xticks(scales);
%         xticklabels(scaleBounds(end:-1:1));
%         xtickangle(45) %xticklabels({f}); %
%     else
%         xlabel('Scale factors')
%     end
%     if max(mfe) < 1
%         ylim([0 1])
%     end
%     title('Multiscale fuzzy entropy - EEG'); ylabel('Entropy')
% 
%     set(findall(gcf,'type','axes'),'fontSize',11,'fontweight','bold');
%     % set(gca,'FontSize',11,'layer','top','fontweight','bold');
% end

%% Scalp topos 2D

if params.eeg && params.eeg_frequency

    mode = 1;  % 1 for 2D, 2 for 3D

    % Create figure
    figure('color','w')
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);    % enlarge
    set(gcf, 'Toolbar', 'none', 'Menu', 'none');                    % remove toolbobar and menu
    set(gcf, 'Name', 'EEG features', 'NumberTitle', 'Off')  % name

    % delta, theta, alpha, beta
    warning('off','all')
    nexttile
    plot_topo(gather(mean(EEG.frequency.delta,2)),params.chanlocs,mode,'psd');
    cb = colorbar; ylabel(cb,'Power (db)','Rotation',270,'fontSize',12,'fontweight','bold')
    % pos = get(cb,'position');  % move left & shrink to match head size
    % set(cb,'position',[pos(1) pos(2)+0.2 0.01 0.4],'fontSize',12,'fontweight','bold');
    title('Delta power'); %colorbar off;
    nexttile
    plot_topo(gather(mean(EEG.frequency.theta,2)),params.chanlocs,mode,'psd');
    cb = colorbar; ylabel(cb,'Power (db)','Rotation',270,'fontSize',12,'fontweight','bold')
    % pos = get(cb,'position');  % move left & shrink to match head size
    % set(cb,'position',[pos(1)-.005 pos(2)+0.05 pos(3)*0.7 pos(4)-0.1],'fontSize',12,'fontweight','bold');
    title('Theta power');
    nexttile
    plot_topo(gather(mean(EEG.frequency.alpha,2)),params.chanlocs,mode,'psd');
    cb = colorbar; ylabel(cb,'Power (db)','Rotation',270,'fontSize',12,'fontweight','bold')
    % pos = get(cb,'position');  % move left & shrink to match head size
    % set(cb,'position',[pos(1)-.005 pos(2)+0.05 pos(3)*0.7 pos(4)-0.1],'fontSize',12,'fontweight','bold');
    title('Alpha power');
    nexttile
    plot_topo(gather(mean(EEG.frequency.beta,2)),params.chanlocs,mode,'psd');
    cb = colorbar; ylabel(cb,'Power (db)','Rotation',270,'fontSize',12,'fontweight','bold')
    % pos = get(cb,'position');  % move left & shrink to match head size
    % set(cb,'position',[pos(1)-.005 pos(2)+0.05 pos(3)*0.7 pos(4)-0.1],'fontSize',12,'fontweight','bold');
    title('Beta power');
    % nexttile
    % plot_topo(gather(mean(EEG.frequency.low_gamma,2)),params.chanlocs,mode,'psd');
    % cb = colorbar; ylabel(cb,'Power (db)','Rotation',270,'fontSize',12,'fontweight','bold')
    % pos = get(cb,'position');  % move left & shrink to match head size
    % set(cb,'position',[pos(1)-.005 pos(2)+0.05 pos(3)*0.7 pos(4)-0.1],'fontSize',12,'fontweight','bold');
    % title('Mean gamma power');
    warning on

    % IAF
    try
        % subplot(3,3,5)
        nexttile
        plot_topo(gather(EEG.frequency.IAF),params.chanlocs,mode,'psd');
        cb = colorbar; ylabel(cb,'Power (db)','Rotation',270,'fontSize',12,'fontweight','bold')
        % pos = get(cb,'position');  % move left & shrink to match head size
        % set(cb,'position',[pos(1)-.005 pos(2)+0.05 pos(3)*0.7 pos(4)-0.1],'fontSize',12,'fontweight','bold');
        title('Individual alpha frequency (IAF)');
    catch
        warning('Could not plot the individual alpha frequency (IAF). IAF estimation may have failed (can happen if no clear alpha peak distribution is present)')
    end

    set(findall(gcf,'type','axes'),'fontSize',12,'fontweight','bold');

end

if params.eeg && params.eeg_nonlinear
    nexttile
    plot_topo(gather(EEG.nonlinear.SE),params.chanlocs,mode,'entropy');
    cb = colorbar; ylabel(cb,'Sample entropy','Rotation',270,'fontSize',12,'fontweight','bold')
    title('Sample entropy');

    % nexttile
    % plot_topo(gather(EEG.nonlinear.FE),params.chanlocs,mode,'entropy');
    % cb = colorbar; ylabel(cb,'Fuzzy entropy','Rotation',270,'fontSize',12,'fontweight','bold')
    % title('Fuzzy entropy');

    nexttile
    plot_topo(gather(EEG.nonlinear.FD),params.chanlocs,mode,'entropy');
    cb = colorbar; ylabel(cb,'Fractal dimension','Rotation',270,'fontSize',12,'fontweight','bold')
    title('Fractal dimension');
end

set(findall(gcf,'type','axes'),'fontSize',12,'fontweight','bold');
% set(findall(gca,'type','axes'),'fontSize',12,'fontweight','bold');
set(gcf,'Toolbar','none','Menu','none');  % remove toolbobar and menu
set(gcf,'Name','Visualization of features','NumberTitle','Off')  % name

% Asymmetry (has to be after ecause of colorbar issues
if params.eeg && params.eeg_frequency
    try
        warning('off','all')
        nexttile
        view = [-85 20];  % 'left'
        asy = EEG.frequency.asymmetry;
        pairNums = EEG.frequency.asymmetry_pairs_num;
        headplotparams = { 'meshfile','mheadnew.mat','transform',...
            [0.664455 -3.39403 -14.2521 -0.00241453 0.015519 -1.55584 11 10.1455 12],...
            'colormap',parula,'maplimits','absmax','cbar',1,...
            'electrodes','off','material','metal','verbose','off'};
        brainbeats_headplot('setup',params.chanlocs(pairNums(:,1)),...
            'tmp.spl',headplotparams{:}); % Generate temporary spline file
        brainbeats_headplot(asy,'tmp.spl','view',view,headplotparams{:});  % 3D headplot of asymmetry
        title('Alpha asymmetry')
        % cb = colorbar; ylabel(cb,'Asymmetry score','Rotation',270,'fontSize',12,'fontweight','bold')
        % nticks = length(cb.Ticks); % number of ticks for colorbar
        % increment = ceil(length(asy)/nticks);
        % increment = 1;
        % ticks = round(sort(asy),2);
        % cb.TickLabels = cellstr(string(ticks(1:increment:end))');
        set(findall(gcf,'type','axes'),'fontSize',12,'fontweight','bold');
        warning on

    catch
        warning('Failed to plot the 3D headplot of alpha asymmetry. This may happen if your EEG data are low-density (i.e., few EEG channels only)')
    end
end

%% Partial EEG coherence

% if params.eeg && params.eeg_frequency
% 
%     % f = EEG.frequency.eeg_coh_f;
%     chanlocs = params.chanlocs;
%     % my_connplot(EEG.frequency.eeg_pcoh_delta,'labels',{chanlocs.labels},'brainimg','off');
% 
%     % Partial coherence
%     figure('color','w');
% %     % subplot(2,2,1)
%       nexttile
%     % imagesc(EEG.frequency.eeg_pcoh_delta);
%     % % labels = {chanlocs.labels};
%     % % plot_corrmatrix(EEG.frequency.eeg_pcoh_delta,labels)
%     % title('Partial coherence - Delta'); colorbar
%     % xticks(1:length(chanlocs)); xticklabels({chanlocs.labels}); xtickangle(45)
%     % yticks(1:length(chanlocs)); yticklabels({chanlocs.labels});
% 
% %     % subplot(2,2,2)
%       nexttile
%     % imagesc(EEG.frequency.eeg_pcoh_theta);
%     % title('Partial coherence - Theta'); colorbar
%     % xticks(1:length(chanlocs)); xticklabels({chanlocs.labels}); xtickangle(45)
%     % yticks(1:length(chanlocs)); yticklabels({chanlocs.labels});
% 
% %     subplot(2,1,1)
%       nexttile
%     imagesc(EEG.frequency.eeg_pcoh_alpha);
%     title('Partial coherence - Alpha'); colorbar
%     xticks(1:length(chanlocs)); xticklabels({chanlocs.labels}); xtickangle(45)
%     yticks(1:length(chanlocs)); yticklabels({chanlocs.labels});
% 
% %     subplot(2,1,2)
%       nexttile
%     imagesc(EEG.frequency.eeg_pcoh_beta);
%     title('Partial coherence - Beta'); colorbar
%     xticks(1:length(chanlocs)); xticklabels({chanlocs.labels}); xtickangle(45)
%     yticks(1:length(chanlocs)); yticklabels({chanlocs.labels});
% 
% end
