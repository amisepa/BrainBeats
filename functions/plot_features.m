function plot_features(Features,params)

if params.hrv
    HRV = Features.HRV;
end
if params.eeg
    EEG = Features.EEG;
end
if ~params.hrv && ~params.eeg
    fprintf('No features to plot. \n');
    return
end

%% PSD AND MFE PLOT

figure('color','w');

% PSD - HRV
if params.hrv && params.hrv_frequency

    subplot(2,2,1); hold on
    pwr = HRV.frequency.pwr; % power already averaged across windows
    freqs = HRV.frequency.pwr_freqs;
    bandNames = {'ULF' 'VLF' 'LF' 'HF'};
    bands = HRV.frequency.bands;

    idx = find(strcmp(bandNames,'ULF'));
    if isfield(HRV.frequency, 'ulf')
        x = freqs >= bands(idx,1) & freqs <= bands(idx,2);
        y = pwr(x);
        area(freqs(x),y,'FaceColor',"#A2142F",'FaceAlpha',.7);
    else
        bandNames(idx) = [];
        bands(idx,:) = [];
    end

    idx = find(strcmp(bandNames,'VLF'));
    if isfield(HRV.frequency, 'vlf')
        x = freqs >= bands(idx,1) & freqs <= bands(idx,2);
        y = pwr(x);
        area(freqs(x),y,'FaceColor',"#D95319",'FaceAlpha',.7)
    else
        bandNames(idx) = [];
        bands(idx,:) = [];
    end

    idx = find(strcmp(bandNames,'LF'));
    if isfield(HRV.frequency, 'lf')
        x = freqs >= bands(idx,1) & freqs <= bands(idx,2);
        y = pwr(x);
        area(freqs(x),y,'FaceColor',"#EDB120",'FaceAlpha',.7)
    else
        bandNames(idx) = [];
        bands(idx,:) = [];
    end

    idx = find(strcmp(bandNames,'HF'));
    if isfield(HRV.frequency, 'hf')
        x = freqs >= bands(idx,1) & freqs <= bands(idx,2);
        y = pwr(x);
        area(freqs(x),y,'FaceColor',"#0072BD",'FaceAlpha',.7)
    else
        bandNames(idx) = [];
        bands(idx,:) = [];
    end
    legend(bandNames)
    % xticklabels(unique(reshape(bands,1,[]))); xtickangle(45);
    axis tight; box on
    xlabel('Frequency (Hz)');
    if params.norm
        ylabel('Power (ms^2 normalized)');
    else
        ylabel('Power (ms^2)');
    end
    title(sprintf('Power spectral density - HRV'))

end

% Multiscale fuzzy entropy (MFE) - HRV
if params.hrv && params.hrv_nonlinear

    subplot(2,2,3); hold on
    mfe = HRV.nonlinear.MFE;
    scales = HRV.nonlinear.MFE_scales;
    area(scales,mfe,'FaceColor',"#A2142F",'FaceAlpha',.7);
    title('Multiscale fuzzy entropy - HRV'); xlabel('Scale factors'); ylabel('Entropy')
    axis tight; box on; grid on
end

% PSD - EEG
if params.eeg && params.eeg_frequency

    subplot(2,2,2); hold on
    pwr = trimmean(EEG.frequency.pwr,20,1); % 20% trimmed mean across channels
    freqs = EEG.frequency.freqs(1,:);
    bands = [0 3; 3 7; 7 13; 13 30; 30 max(freqs)];
    % area(freqs,pwr,'FaceColor',[0 0.4470 0.7410],'FaceAlpha',.7);

    % delta
    x = freqs >= bands(1,1) & freqs <= bands(1,2);
    y = pwr(x);
    area(freqs(x),y,'FaceColor',"#A2142F",'FaceAlpha',.7)

    % theta
    x = freqs >= bands(2,1) & freqs <= bands(2,2);
    y = pwr(x);
    area(freqs(x),y,'FaceColor',"#D95319",'FaceAlpha',.7)

    % alpha
    x = freqs >= bands(3,1) & freqs <= bands(3,2);
    y = pwr(x);
    area(freqs(x),y,'FaceColor',"#EDB120",'FaceAlpha',.7)

    % beta
    x = freqs >= bands(4,1) & freqs <= bands(4,2);
    y = pwr(x);
    area(freqs(x),y,'FaceColor',"#0072BD",'FaceAlpha',.7)

    % gamma
    x = freqs >= bands(5,1) & freqs <= bands(5,2);
    y = pwr(x);
    area(freqs(x),y,'FaceColor',"#4DBEEE",'FaceAlpha',.7)

    % bandNames = {'Delta Theta Alpha Beta'};
    % xticks = freqs;
    % xticklabels(unique(reshape(bands,1,[])));
    % xtickangle(45); axis tight; box on
    title('Power spectral density - EEG'); legend({'delta' 'theta' 'alpha' 'beta' 'gamma'})
    xlabel('Frequency (Hz)'); ylabel('Power (ms^2)'); axis tight;
end

% Multiscale fuzzy entropy (MFE) - EEG
if params.eeg && params.eeg_nonlinear

    subplot(2,2,4); %hold on
    scales = EEG.nonlinear.MFE_scales(1,:);
    scaleBounds = EEG.nonlinear.MFE_scale_bounds(1,:);
    % mfe = mean(EEG.nonlinear.MFE,1); % mean across channels
    mfe = trimmean(EEG.nonlinear.MFE,20,1); % 20% trimmed mean across channels
    % mfe = EEG.nonlinear.MFE(31,:); % Pz
    area(scales,mfe(end:-1:1),'FaceColor',"#A2142F",'FaceAlpha',.7);
    axis tight; box on; grid on
    if ~isempty(scaleBounds)
        xticks(scales);
        xticklabels(scaleBounds(end:-1:1));
        xtickangle(45) %xticklabels({f}); %
    else
        xlabel('Scale factors')
    end
    if max(mfe) < 1
        ylim([0 1])
    end
    title('Multiscale fuzzy entropy - EEG'); ylabel('Entropy')

    set(findall(gcf,'type','axes'),'fontSize',11,'fontweight','bold');
    % set(gca,'FontSize',11,'layer','top','fontweight','bold');
end

%% Scalp topos 2D

if params.eeg && params.eeg_frequency
    mode = 1;
    figure('color','w')
    subplot(3,2,1)
    plot_topo(gather(mean(EEG.frequency.delta,2)),params.chanlocs,mode,'psd');
    title('Delta power'); %colorbar off;
    subplot(3,2,2)
    plot_topo(gather(mean(EEG.frequency.theta,2)),params.chanlocs,mode,'psd');
    title('Theta power');
    subplot(3,2,3)
    plot_topo(gather(mean(EEG.frequency.alpha,2)),params.chanlocs,mode,'psd');
    title('Alpha power');
    subplot(3,2,4)
    plot_topo(gather(mean(EEG.frequency.beta,2)),params.chanlocs,mode,'psd');
    title('Beta power');
    subplot(3,2,5)
    try
        plot_topo(gather(EEG.frequency.IAF),params.chanlocs,mode,'IAF');
        title('Individual alpha frequency (IAF)'); %colorbar off;
    catch
        warning('Could not plot IAF. IAF estimation may have failed (can happen when no clear alpha peak distribution is present)')
    end
    if params.eeg_nonlinear
        subplot(3,2,6)
        plot_topo(gather(EEG.nonlinear.FE),params.chanlocs,mode,'entropy');
        title('Fuzzy entropy'); %colorbar off;
    end
    
    set(findall(gcf,'type','axes'),'fontSize',10,'fontweight','bold');
end

%% Asymmetry (work in progress)



%% Coherence

if params.eeg && params.eeg_frequency

    try
        f = EEG.frequency.eeg_coh_f;
        chanlocs = params.chanlocs;
    
        % Coherence (measures coupling, symmetric)
        figure('color','w');
        COH = EEG.frequency.eeg_pcoh;
    
        subplot(2,2,1)
        title('Coherence - Delta');
        coh = mean(COH(:,:,f>0 & f<=3),3);
        % plotconnectivity(coh_delta,'labels',{chanlocs.labels},'brainimg','off');
        my_connplot(coh,'labels',{chanlocs.labels},'brainimg','off','threshold',0);
    
        subplot(2,2,2)
        title('Coherence - Theta')
        coh = mean(COH(:,:,f>=3 & f<=7),3);
        my_connplot(coh,'labels',{chanlocs.labels},'brainimg','off','threshold',0);
    
        subplot(2,2,3)
        title('Coherence - Alpha')
        coh = mean(COH(:,:,f>=8 & f<=13),3);
        my_connplot(coh,'labels',{chanlocs.labels},'brainimg','off','threshold',0);
    
        subplot(2,2,4)
        title('Coherence - Beta')
        coh = mean(COH(:,:,f>=14 & f<=30),3);
        my_connplot(coh,'labels',{chanlocs.labels},'brainimg','off','threshold',0);
    
        % Directed Coherence (DC; measures causality, directionality)
        DC = EEG.frequency.eeg_dc;
        figure('color','w');
        subplot(2,2,1)
        title('Directed coherence - Delta');
        dc = mean(DC(:,:,f>0 & f<=3),3);
        my_connplot(dc,'labels',{chanlocs.labels},'brainimg','off','threshold',0);
    
        subplot(2,2,2)
        title('Directed coherence - Theta')
        dc = mean(DC(:,:,f>=3 & f<=7),3);
        my_connplot(dc,'labels',{chanlocs.labels},'brainimg','off','threshold',0);
    
        subplot(2,2,3)  % alpha
        title('Directed coherence - Alpha')
        dc = mean(DC(:,:,f>=8 & f<=13),3);
        my_connplot(dc,'labels',{chanlocs.labels},'brainimg','off','threshold',0);
    
        subplot(2,2,4)  % beta
        title('Directed coherence - Beta')
        dc = mean(DC(:,:,f>=14 & f<=30),3);
        my_connplot(dc,'labels',{chanlocs.labels},'brainimg','off','threshold',0);
    catch
        warning('EEG coherence plots failed.')
    end

    % topoplot(coh, chanlocs, 'emarker2',{[1 2],'b','r'}); 

end
