    %%%% EEG coherence (magnitude squared coherence estimate) %%%%%
    % (only on pairs on electrodes that are not neighbors; see Nunez 2016:
    % https://escholarship.org/content/qt5fb0q5wx/qt5fb0q5wx.pdf)
    % FIXME: neighbors are incorrect with low-density montages, use actual
    % distance instead?.
    % nfft = fs*2;
    % noverlap = nfft/2;

    % neighbor electrodes
    % neighbors = get_channelneighbors(chanlocs);

    % all possible pairs
    % pairs = nchoosek({chanlocs.labels}, 2);

    % remove pairs that are neighbors
    % nPairs = size(pairs,1);
    % pairname = cell(nPairs,1);
    % cohr = nan(nPairs,nfft/2+1);
    % progressbar('Estimating EEG coherence on each (non-neighbor) channel pair')
    % fprintf('Estimating EEG coherence on all possible channel pairs... \n')
    % warning('Ignoring neighboring channels that are contaminated by volume conduction effects (see Nunez et al. 2016, Figure 7)');
    % for iPair = 1:nPairs
    %
    %     chan1 = pairs{iPair,1};
    %     chan2 = pairs{iPair,2};
    %     chan1_neighbors = neighbors(strcmp({neighbors.label},chan1)).neighblabel;
    %     pairname{iPair,:} = sprintf('%s - %s', pairs{iPair,1}, pairs{iPair,2});
    %
    %     % If chan2 is neighbor, skip to next pair, otherwise compute coherence
    %     if sum(contains(chan1_neighbors, chan2)) == 0
    %         fprintf('pair: %s \n', pairname{iPair,:})
    %         idx1 = strcmp({chanlocs.labels},chan1);
    %         idx2 = strcmp({chanlocs.labels},chan2);
    %         [cohr(iPair,:),f] = mscohere(signals(idx1,:),signals(idx2,:),hamming(nfft),noverlap,nfft,fs);
    %         % plot(f(f>=0 & f<45), squeeze(cohr(f>=0 & f<45)));
    %         % title(pairname); grid on; hold on;
    %     else
    %         fprintf('pair: %s (neighbors -> ignored) \n', pairname{iPair,:})
    %         continue
    %     end
    %
    %     progressbar(iPair/nPairs);
    % end
    %
    % % Remove empty rows
    % nans = isnan(cohr(:,1));
    % cohr(nans,:) = [];
    % pairname(nans,:) = [];
    %
    % % Export
    % eeg_features.frequency.eeg_coherence = cohr;
    % eeg_features.frequency.eeg_coherence_pair = pairname;
    % eeg_features.frequency.eeg_coherence_f = f;
    % fprintf('Coherence estimated on %g pairs after excluding neighbors. \n', length(pairname));

    % MULTIVARIATE ANALYSIS OF CAUSAL INTERACTIONS (see Faes and Nollo; 2011).
    % https://www.intechopen.com/chapters/12918
    % Here, we extract coherence (coh), partial coherence (pcoh), directed
    % coherence (DC) and partial directed coherence (PDC).
    % See below for descriptions.
    nfft = fs*2;  % 2-s windows
    fprintf('Computing multivariate analysis of causal interactions between EEG channels... \n')
    Su = eye(nChan,nChan);
    [dc,~,pdc,gpdc,~,coh,pcoh,~,~,~,~,f] = fdMVAR_5order(signals,Su,nfft,fs);

    % Deal with complex and negative values
    if ~isreal(coh)
        coh = abs(coh);
    end
    if ~isreal(pcoh)
        pcoh = abs(pcoh);
    end
    if ~isreal(dc)
        dc = abs(dc);  % abs to convert everything to positive, real to preserve polarity
    end
    if ~isreal(pdc)
        pdc = abs(pdc);
    end
    if ~isreal(gpdc)
        gpdc = abs(gpdc);
    end
    eeg_features.frequency.eeg_coh_f = f;
    eeg_features.frequency.eeg_coh = coh;
    eeg_features.frequency.eeg_pcoh = pcoh;
    eeg_features.frequency.eeg_dc = dc;
    eeg_features.frequency.eeg_pdc = pdc;

    % COHERENCE PER BAND
    % Coherence quantifies the linear correlation between two signals in the
    % frequency domain, i.e. the consistency of the phase difference and the
    % correlation of the amplitudes of the two signals across time windows.
    idx = f<=3;
    coh_delta = tril(mean(coh(:,:,idx),3,'omitnan'),-1);    % zero upper triangle
    coh_delta(logical(eye(size(coh_delta)))) = 1;           % one diagonal
    % imagesc(coh_delta); colorbar
    % labels = {chanlocs.labels};
    % plot_corrmatrix(coh_delta,labels)
    idx = f>=3 & f<8;
    coh_theta = tril(mean(coh(:,:,idx),3,'omitnan'),-1);    % zero upper triangle
    coh_theta(logical(eye(size(coh_theta)))) = 1;           % one diagonal
    % imagesc(coh_theta); colorbar
    idx = f>=8 & f<=13;
    coh_alpha = tril(mean(coh(:,:,idx),3,'omitnan'),-1);    % zero upper triangle
    coh_alpha(logical(eye(size(coh_alpha)))) = 1;           % one diagonal
    % imagesc(coh_alpha); colorbar
    idx = f>13 & f<=30;
    coh_beta = tril(mean(coh(:,:,idx),3,'omitnan'),-1);     % zero upper triangle
    coh_beta(logical(eye(size(coh_beta)))) = 1;             % one diagonal
    % imagesc(coh_beta); colorbar
    idx = f>30 & f<=50;
    coh_lowgamma = tril(mean(coh(:,:,idx),3,'omitnan'),-1); % zero upper triangle
    coh_lowgamma(logical(eye(size(coh_lowgamma)))) = 1;     % one diagonal
    % imagesc(coh_lowgamma); colorbar
    eeg_features.frequency.eeg_coh_delta = coh_delta;
    eeg_features.frequency.eeg_coh_theta = coh_theta;
    eeg_features.frequency.eeg_coh_alpha = coh_alpha;
    eeg_features.frequency.eeg_coh_beta = coh_beta;
    eeg_features.frequency.eeg_coh_lowgamma = coh_lowgamma;

    % PARTIAL COHERENCE PER BAND
    % Coherence might indicate that two brain areas are communicating or are
    % influenced by a common source. However, this doesn't necessarily mean
    % they're communicating directly. For instance, 2 regions might show
    % high coherence because they both receive input from a third region.
    % pcoh resolves this ambiguity. If two regions still show high partial
    % coherence after the influence of the third signal is accounted for,
    % it suggests a more direct relationship between them. i.e. it helps
    % determine if the observed coherence is direct or if it's mediated
    % through other signals.
    idx = f<=3;
    pcoh_delta = tril(mean(pcoh(:,:,idx),3,'omitnan'),-1);    % zero upper triangle
    pcoh_delta(logical(eye(size(pcoh_delta)))) = 1;           % one diagonal
    % imagesc(pcoh_delta); colorbar
    idx = f>=3 & f<8;
    pcoh_theta = tril(mean(pcoh(:,:,idx),3,'omitnan'),-1);    % zero upper triangle
    pcoh_theta(logical(eye(size(pcoh_theta)))) = 1;           % one diagonal
    % imagesc(pcoh_theta); colorbar
    idx = f>=8 & f<=13;
    pcoh_alpha = tril(mean(pcoh(:,:,idx),3,'omitnan'),-1);    % zero upper triangle
    pcoh_alpha(logical(eye(size(pcoh_alpha)))) = 1;           % one diagonal
    % figure; imagesc(pcoh_alpha); colorbar
    % labels = {chanlocs.labels};
    % figure; plot_corrmatrix(pcoh_alpha,labels)
    idx = f>13 & f<=30;
    pcoh_beta = tril(mean(pcoh(:,:,idx),3,'omitnan'),-1);     % zero upper triangle
    pcoh_beta(logical(eye(size(pcoh_beta)))) = 1;             % one diagonal
    % imagesc(pcoh_beta); colorbar
    idx = f>30 & f<=50;
    pcoh_lowgamma = tril(mean(pcoh(:,:,idx),3,'omitnan'),-1); % zero upper triangle
    pcoh_lowgamma(logical(eye(size(pcoh_lowgamma)))) = 1;     % one diagonal
    % imagesc(abs(pcoh_lowgamma)); colorbar
    eeg_features.frequency.eeg_pcoh_delta = pcoh_delta;
    eeg_features.frequency.eeg_pcoh_theta = pcoh_theta;
    eeg_features.frequency.eeg_pcoh_alpha = pcoh_alpha;
    eeg_features.frequency.eeg_pcoh_beta = pcoh_beta;
    eeg_features.frequency.eeg_pcoh_lowgamma = pcoh_lowgamma;

    % DIRECTED COHERENCE PER BAND
    % Provides information about directionality, indicating which signal has
    % a causal influence on the other, useful to identify the flow of
    % information in neural circuits. Based on Granger Causality.
    idx = f<=3;
    dc_delta = tril(mean(dc(:,:,idx),3,'omitnan'),-1);    % zero upper triangle
    dc_delta(logical(eye(size(dc_delta)))) = 1;           % one diagonal
    % imagesc(dc_delta); colorbar
    idx = f>=3 & f<8;
    dc_theta = tril(mean(dc(:,:,idx),3,'omitnan'),-1);    % zero upper triangle
    dc_theta(logical(eye(size(dc_theta)))) = 1;           % one diagonal
    % imagesc(dc_theta); colorbar
    idx = f>=8 & f<=13;
    dc_alpha = tril(mean(dc(:,:,idx),3,'omitnan'),-1);    % zero upper triangle
    dc_alpha(logical(eye(size(dc_alpha)))) = 1;           % one diagonal
    % figure;imagesc(dc_alpha); colorbar
    % labels = {chanlocs.labels};
    % plot_corrmatrix(dc_alpha,labels)
    idx = f>13 & f<=30;
    dc_beta = tril(mean(dc(:,:,idx),3,'omitnan'),-1);     % zero upper triangle
    dc_beta(logical(eye(size(dc_beta)))) = 1;             % one diagonal
    % imagesc(dc_beta); colorbar
    idx = f>30 & f<=50;
    dc_lowgamma = tril(mean(dc(:,:,idx),3,'omitnan'),-1); % zero upper triangle
    dc_lowgamma(logical(eye(size(dc_lowgamma)))) = 1;     % one diagonal
    % imagesc(dc_lowgamma); colorbar
    eeg_features.frequency.eeg_dc_delta = dc_delta;
    eeg_features.frequency.eeg_dc_theta = dc_theta;
    eeg_features.frequency.eeg_dc_alpha = dc_alpha;
    eeg_features.frequency.eeg_dc_beta = dc_beta;
    eeg_features.frequency.eeg_dc_lowgamma = dc_lowgamma;

    % PARTIAL DIRECTED COHERENCE PER BAND
    % By accounting for all other signals, it provides a more refined measure
    % of the direct influence of one signal over another. If there's a strong
    % PDC from signal X to signal Y, it suggests that X has a direct influence
    % on Y, not mediated through other signals.
    % DC = total influence (direct + mediated), PDC = direct influenc only.
    idx = f<=3;
    pdc_delta = tril(mean(pdc(:,:,idx),3,'omitnan'),-1);    % zero upper triangle
    pdc_delta(logical(eye(size(pdc_delta)))) = 1;           % one diagonal
    % imagesc(pdc_delta); colorbar
    idx = f>=3 & f<8;
    pdc_theta = tril(mean(pdc(:,:,idx),3,'omitnan'),-1);    % zero upper triangle
    pdc_theta(logical(eye(size(pdc_theta)))) = 1;           % one diagonal
    % imagesc(pdc_theta); colorbar
    idx = f>=8 & f<=13;
    pdc_alpha = tril(mean(pdc(:,:,idx),3,'omitnan'),-1);    % zero upper triangle
    pdc_alpha(logical(eye(size(pdc_alpha)))) = 1;           % one diagonal
    % figure; imagesc(pdc_alpha); colorbar
    % labels = {chanlocs.labels};
    % plot_corrmatrix(pdc_alpha,labels)
    idx = f>13 & f<=30;
    pdc_beta = tril(mean(pdc(:,:,idx),3,'omitnan'),-1);     % zero upper triangle
    pdc_beta(logical(eye(size(pdc_beta)))) = 1;             % one diagonal
    % imagesc(pdc_beta); colorbar
    idx = f>30 & f<=50;
    pdc_lowgamma = tril(mean(pdc(:,:,idx),3,'omitnan'),-1); % zero upper triangle
    pdc_lowgamma(logical(eye(size(pdc_lowgamma)))) = 1;     % one diagonal
    % imagesc(pdc_lowgamma); colorbar
    eeg_features.frequency.eeg_pdc_delta = pdc_delta;
    eeg_features.frequency.eeg_pdc_theta = pdc_theta;
    eeg_features.frequency.eeg_pdc_alpha = pdc_alpha;
    eeg_features.frequency.eeg_pdc_beta = pdc_beta;
    eeg_features.frequency.eeg_pdc_lowgamma = pdc_lowgamma;

    % GENERALIZED PARTIAL DIRECTED COHERENCE PER BAND
    % Extension of PDC that factors in the input noise covariance when
    % determining causality. This normalization is considering the relative
    % strength of the direct connections in the presence of background noise
    % or input variances. It provides a more detailed picture of the dynamics,
    % especially when the system's components (e.g., EEG channels or brain regions)
    % have different noise levels or variances. especially useful in scenarios
    % where there's variability in the noise levels or input variances across signals or brain regions.
    idx = f<=3;
    gpdc_delta = tril(mean(gpdc(:,:,idx),3,'omitnan'),-1);    % zero upper triangle
    gpdc_delta(logical(eye(size(gpdc_delta)))) = 1;           % one diagonal
    % figure; imagesc(gpdc_delta); colorbar
    idx = f>=3 & f<8;
    gpdc_theta = tril(mean(gpdc(:,:,idx),3,'omitnan'),-1);    % zero upper triangle
    gpdc_theta(logical(eye(size(gpdc_theta)))) = 1;           % one diagonal
    % figure; imagesc(gpdc_theta); colorbar
    idx = f>=8 & f<=13;
    gpdc_alpha = tril(mean(gpdc(:,:,idx),3,'omitnan'),-1);    % zero upper triangle
    gpdc_alpha(logical(eye(size(gpdc_alpha)))) = 1;           % one diagonal
    % figure; imagesc(gpdc_alpha); colorbar
    % labels = {chanlocs.labels};
    % plot_corrmatrix(gpdc_alpha,labels)
    idx = f>13 & f<=30;
    gpdc_beta = tril(mean(gpdc(:,:,idx),3,'omitnan'),-1);     % zero upper triangle
    gpdc_beta(logical(eye(size(gpdc_beta)))) = 1;             % one diagonal
    % figure; imagesc(gpdc_beta); colorbar
    idx = f>30 & f<=50;
    gpdc_lowgamma = tril(mean(gpdc(:,:,idx),3,'omitnan'),-1); % zero upper triangle
    gpdc_lowgamma(logical(eye(size(gpdc_lowgamma)))) = 1;     % one diagonal
    % figure; imagesc(gpdc_lowgamma); colorbar
    eeg_features.frequency.eeg_gpdc_delta = gpdc_delta;
    eeg_features.frequency.eeg_gpdc_theta = gpdc_theta;
    eeg_features.frequency.eeg_gpdc_alpha = gpdc_alpha;
    eeg_features.frequency.eeg_gpdc_beta = gpdc_beta;
    eeg_features.frequency.eeg_gpdc_lowgamma = gpdc_lowgamma;

    warning('Causal relationships between neighboring channels should be interpreted with caution due to volume conduction effects. See for example https://escholarship.org/content/qt5fb0q5wx/qt5fb0q5wx.pdf')

