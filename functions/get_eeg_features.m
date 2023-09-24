%% Extract EEG features in time, fequency, and nonlinear domains.
%
% Copyright (C) - Cedric Cannard, 2023

function eeg_features = get_eeg_features(signals,params)

disp('----------------------------------------------------')
disp('               Extracting EEG features ')
disp('----------------------------------------------------')

tstart = tic;

fs = params.fs;
chanlocs = params.chanlocs;

% Parallel computing
if params.parpool
    % ps = parallel.Settings; 
    % ps.Pool.AutoCreate = params.parpool; 

    % delete(gcp('nocreate')) %shut down opened parpool
    p = gcp('nocreate');
    if isempty(p) % if not already on, launch it
        disp('Initiating parrallel computing (all cores and threads -1)...')
        c = parcluster; % cluster profile
        % N = feature('numcores');        % physical number of cores
        N = getenv('NUMBER_OF_PROCESSORS'); % all processors (including threads)
        if ischar(N), N = str2double(N); end
        c.NumWorkers = N-1;  % update cluster profile to include all workers
        c.parpool();
    end
end

% Use parallel GPUs computing (if multiple GPUS are available)
if params.parpool && params.gpu
    availableGPUs = gpuDeviceCount("available");
    if availableGPUs > 1
        parpool('Processes',availableGPUs);
        fprintf('%g GPUs detected. Using them in parallel pool. \n',availableGPUs)
    else
        fprintf('Only one GPU detected. Using normal GPU and parallel pool computing. \n')
    end
end

%% Pre-allocate matrices: for memory and especially important if bad 
% channels were not interpolated)



%% Time domain

if params.eeg_time
    disp('Calculating time-domain EEG features...')
    eeg_features.time.rms = rms(signals,2);
    eeg_features.time.mode = mode(signals,2);
    eeg_features.time.var = var(signals,0,2);
    eeg_features.time.skewness = skewness(signals,0,2);
    eeg_features.time.kurtosis = kurtosis(signals,0,2);
    eeg_features.time.iqr = iqr(signals,2);
end

%% Frequency domain

if params.eeg_frequency

    nChan = size(signals,1);
    fRange = [1 40];    % FIXME: lowpass/highpass should be in params to make sure these are within filtered signal
    winSize = 4;        % window size (in s). Default = 2 (at least 2 s recommended by Smith et al, 2017 for asymmetry)
    winType = 'hamming';
    overlap = 50;       % 50% default (Smith et al. 2017)

    % progressbar (only when not in parpool)
    if ~params.parpool
        progressbar('Estimating EEG power spectral density on each channel')
    end
    disp('Calculating band-power on each EEG channel')
    for iChan = 1:nChan

        fprintf('EEG channel %g \n', iChan)

        % Compute PSD using pwelch
        [pwr, pwr_dB, f] = compute_psd(signals(iChan,:),fs*winSize,winType,overlap,[],fs,fRange,'psd',params.gpu);
        eeg_features.frequency.pwr(iChan,:) = pwr;
        eeg_features.frequency.pwr_dB(iChan,:) = pwr_dB;
        eeg_features.frequency.freqs(iChan,:) = f;

        % Delta
        eeg_features.frequency.delta(iChan,:) = mean( pwr_dB(f >= f(1) & f <= 3) );
        eeg_features.frequency.delta_norm(iChan,:) = eeg_features.frequency.delta(iChan,:) ./ sum(pwr_dB); % normalized by total power of same channel
        try
            bounds = get_freqBounds(pwr, f, fs, [f(1) 3], winSize,0.25);  % individualized frequency bounds
            eeg_features.frequency.delta_indiv(iChan,:) = mean( pwr_dB(f >= bounds(1) & f <= bounds(2)) );
        catch
            fprintf('Could not detect individualized frequency bounds for Delta \n')
        end

        % Theta
        eeg_features.frequency.theta(iChan,:) = mean( pwr_dB(f >= 3 & f <= 7) );
        eeg_features.frequency.theta_norm(iChan,:) = eeg_features.frequency.theta(iChan,:) ./ sum(pwr_dB); % normalized by total power of same channel
        try
            bounds = get_freqBounds(pwr, f, fs, [2 8], winSize, 0.25);  % individualized lower/upper bounds
            eeg_features.frequency.theta_indiv(iChan,:) = mean( pwr_dB(f >= bounds(1) & f <= bounds(2)) );
        catch
            fprintf('Could not detect individualized frequency bounds for Theta \n')
        end

        % Alpha
        eeg_features.frequency.alpha(iChan,:) = mean( pwr_dB(f >= 7.5 & f <= 13) );
        eeg_features.frequency.alpha_norm(iChan,:) = eeg_features.frequency.alpha(iChan,:) ./ sum(pwr_dB); % normalized by total power of same channel
        try
            bounds = get_freqBounds(pwr, f, fs, [7 14], winSize, 1);  % individualized lower/upper bounds
            eeg_features.frequency.alpha_indiv(iChan,:) = mean( pwr_dB(f >= bounds(1) & f <= bounds(2)) );
        catch
            fprintf('Could not detect individualized frequency bounds for Alpha \n')
        end

        % Beta
        eeg_features.frequency.beta(iChan,:) = mean( pwr_dB(f >= 13.5 & f <= 30) );
        eeg_features.frequency.beta_norm(iChan,:) = eeg_features.frequency.beta(iChan,:) ./ sum(pwr_dB); % normalized by total power of same channel
        try
            bounds = get_freqBounds(pwr, f, fs, [13 30], winSize, 0.25);  % individualized lower/upper bounds
            eeg_features.frequency.beta_indiv(iChan,:) = mean( pwr_dB(f >= bounds(1) & f <= bounds(2)) );
        catch
            fprintf('Could not detect individualized frequency bounds for Beta \n')
        end

        % Low gamma
        eeg_features.frequency.low_gamma(iChan,:) = mean( pwr_dB(f >= 31 & f <= fRange(2)) );
        eeg_features.frequency.low_gamma_norm(iChan,:) = eeg_features.frequency.low_gamma(iChan,:) ./ sum(pwr_dB); % normalized by total power of same channel
        try
            bounds = get_freqBounds(pwr, f, fs, [30 45], winSize, 0.25);  % individualized lower/upper bounds
            eeg_features.frequency.low_gamma_indiv(iChan,:) = mean( pwr_dB(f >= bounds(1) & f <= bounds(2)) );
        catch
            fprintf('Could not detect individualized frequency bounds for Low gamma \n')
        end

        % Individual alpha frequency (IAF) (my code, not working)
        % iaf = detect_iaf(pwr(iChan,:), freqs, winSize, params)

        if ~params.parpool
            progressbar(iChan/nChan)
        end

    end

    % Individual alpha frequency (IAF; only export CoG since it's the best)
    disp('Detecting individual alpha frequency (IAF) for each EEG channel...')
    [pSum, pChans, f] = restingIAF(signals, size(signals,1), 1, [1 30], fs, [7 13], 11, 5);
    eeg_features.frequency.IAF_mean = pSum.cog;
    eeg_features.frequency.IAF = [pChans.gravs]';

    % Asymmetry (use log(pwr) no pwr_dB) - on all pairs
    nChan = size(chanlocs,2);
    % if rem(nChan,2) == 0 %even
    nPairs = floor(nChan/2); 
    % % else
    % %     nPairs = floor(nChan/2+1);
    % end 
    % pairs = nan(nPairs,2);
    % pairLabels = cell(nPairs,1);
    % asy = nan(nPairs,1);
    fprintf('Extracting (z-normalized) EEG asymmetry on %g channel pairs... \n',nPairs)
    for iChan = 1:nChan
        
        % if not left electrode, skip
        if chanlocs(iChan).theta >= 0 
            % pairs(iChan,:) = [NaN NaN];
            % pairLabels(iChan,:) = {' '  ' '};
            continue; 
        end

        % find matching pair from right side using theta distance
        for iChan2 = 1:nChan
            if iChan2 == iChan
                elec_dist(iChan2,:) = NaN;
            else
                elec_dist(iChan2,:) = diff([abs(chanlocs(iChan).theta) abs(chanlocs(iChan2).theta) ]);
            end
        end
        [~, match] = min(abs(elec_dist));
        pairs(iChan,:) = [iChan match];
        pairLabels(iChan,:) = { sprintf('%s %s', chanlocs(iChan).labels, chanlocs(match).labels) };
        
        % Ignore midline electrodes
        % if contains(pairLabels{iChan}, 'z'), continue, end

        % Ensure pairs are always left-right to remove duplicates
        % if rem(str2double(pairLabels{iChan}(end)),2) ~= 0 % second elec should be even number
        %     pairs(iChan,:) = [match iChan];
        %     pairLabels(iChan,:) = { sprintf('%s %s', chanlocs(match).labels, chanlocs(iChan).labels) };
        % end
    end
    
    % Remove duplicate or empty pairs
    % pairLabels(pairs(:,1) == 0) = [];
    % pairs(pairs(:,1) == 0,:) = [];
    % [pairLabels, idx] = unique(pairLabels);
    % pairs = pairs(idx,:);

    % asy = nan(length(pairs),length(f));
    for iPair = 1:size(pairs,1)
        % Z-normalize by correcting for overall alpha power (see Allen et al. 2004 and Smith et al. 2017)
        alpha_left = mean(eeg_features.frequency.alpha(pairs(iPair,1),:));
        alpha_right = mean(eeg_features.frequency.alpha(pairs(iPair,2),:));
        alpha_left = alpha_left / sum(mean(eeg_features.frequency.alpha,2));
        alpha_right = alpha_right / sum(mean(eeg_features.frequency.alpha,2));

        % compute log asymmetry and export
        asy(iPair,:) = log(alpha_left) - log(alpha_right);
    end
    eeg_features.frequency.asymmetry = asy(~isnan(asy));
    eeg_features.frequency.asymmetry_pairs = pairLabels(~isnan(asy));

    % if params.vis
        % topoplot(asy, chanlocs, 'emarker2',{[3 17],'b','r'}); 
    % end

    %%%% EEG coherence (magnitude squared coherence estimate) %%%%%
    % (only on pairs on electrodes that are not neighbors; see Nunez 2016:
    % https://escholarship.org/content/qt5fb0q5wx/qt5fb0q5wx.pdf)
    % FIXME: neighbors are incorrect with low-density montages, should use
    % actual distance instead?. 
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
    
    % MULTIVARIATE ANALYSIS OF CAUSAL INTERACTIONS (see Faes and Nollo (2011). 
    % https://www.intechopen.com/chapters/12918
    % Here, we extract coherence (coh), partial coherence (pcoh), directed
    % coherence (DC) and partial directed coherence (PDC). 
    % See below for descriptions. 
    nfft = fs*2;  % 2-s windows
    fprintf('Computing multivariate analysis of causal interactions between all EEG channels... \n')
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
    % imagesc(abs(coh_delta)); colorbar
    % labels = {chanlocs.labels};
    % plot_corrmatrix(coh_delta,labels)
    idx = f>=3 & f<8;
    coh_theta = tril(mean(coh(:,:,idx),3,'omitnan'),-1);    % zero upper triangle
    coh_theta(logical(eye(size(coh_theta)))) = 1;           % one diagonal
    % imagesc(abs(coh_theta)); colorbar
    idx = f>=8 & f<=13;
    coh_alpha = tril(mean(coh(:,:,idx),3,'omitnan'),-1);    % zero upper triangle
    coh_alpha(logical(eye(size(coh_alpha)))) = 1;           % one diagonal
    % imagesc(abs(coh_alpha)); colorbar
    idx = f>13 & f<=30;
    coh_beta = tril(mean(coh(:,:,idx),3,'omitnan'),-1);     % zero upper triangle
    coh_beta(logical(eye(size(coh_beta)))) = 1;             % one diagonal
    % imagesc(abs(coh_beta)); colorbar
    idx = f>30 & f<=50;
    coh_lowgamma = tril(mean(coh(:,:,idx),3,'omitnan'),-1); % zero upper triangle
    coh_lowgamma(logical(eye(size(coh_lowgamma)))) = 1;     % one diagonal
    % imagesc(abs(coh_lowgamma)); colorbar
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
    % imagesc(abs(pcoh_delta)); colorbar
    idx = f>=3 & f<8;
    pcoh_theta = tril(mean(pcoh(:,:,idx),3,'omitnan'),-1);    % zero upper triangle
    pcoh_theta(logical(eye(size(pcoh_theta)))) = 1;           % one diagonal
    % imagesc(abs(pcoh_theta)); colorbar
    idx = f>=8 & f<=13;
    pcoh_alpha = tril(mean(pcoh(:,:,idx),3,'omitnan'),-1);    % zero upper triangle
    pcoh_alpha(logical(eye(size(pcoh_alpha)))) = 1;           % one diagonal
    figure; imagesc(abs(pcoh_alpha)); colorbar
    labels = {chanlocs.labels};
    figure; plot_corrmatrix(pcoh_alpha,labels)
    idx = f>13 & f<=30;
    pcoh_beta = tril(mean(pcoh(:,:,idx),3,'omitnan'),-1);     % zero upper triangle
    pcoh_beta(logical(eye(size(pcoh_beta)))) = 1;             % one diagonal
    % imagesc(abs(pcoh_beta)); colorbar
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
    % imagesc(abs(dc_delta)); colorbar
    idx = f>=3 & f<8;
    dc_theta = tril(mean(dc(:,:,idx),3,'omitnan'),-1);    % zero upper triangle
    dc_theta(logical(eye(size(dc_theta)))) = 1;           % one diagonal
    % imagesc(abs(dc_theta)); colorbar
    idx = f>=8 & f<=13;
    dc_alpha = tril(mean(dc(:,:,idx),3,'omitnan'),-1);    % zero upper triangle
    dc_alpha(logical(eye(size(dc_alpha)))) = 1;           % one diagonal
    % figure;imagesc(abs(dc_alpha)); colorbar
    labels = {chanlocs.labels};
    plot_corrmatrix(dc_alpha,labels)
    idx = f>13 & f<=30;
    dc_beta = tril(mean(dc(:,:,idx),3,'omitnan'),-1);     % zero upper triangle
    dc_beta(logical(eye(size(dc_beta)))) = 1;             % one diagonal
    % imagesc(abs(dc_beta)); colorbar
    idx = f>30 & f<=50;
    dc_lowgamma = tril(mean(dc(:,:,idx),3,'omitnan'),-1); % zero upper triangle
    dc_lowgamma(logical(eye(size(dc_lowgamma)))) = 1;     % one diagonal
    % imagesc(abs(dc_lowgamma)); colorbar
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
    % imagesc(abs(pdc_delta)); colorbar
    idx = f>=3 & f<8;
    pdc_theta = tril(mean(pdc(:,:,idx),3,'omitnan'),-1);    % zero upper triangle
    pdc_theta(logical(eye(size(pdc_theta)))) = 1;           % one diagonal
    % imagesc(abs(pdc_theta)); colorbar
    idx = f>=8 & f<=13;
    pdc_alpha = tril(mean(pdc(:,:,idx),3,'omitnan'),-1);    % zero upper triangle
    pdc_alpha(logical(eye(size(pdc_alpha)))) = 1;           % one diagonal
    % figure; imagesc(abs(pdc_alpha)); colorbar
    labels = {chanlocs.labels};
    plot_corrmatrix(pdc_alpha,labels)
    idx = f>13 & f<=30;
    pdc_beta = tril(mean(pdc(:,:,idx),3,'omitnan'),-1);     % zero upper triangle
    pdc_beta(logical(eye(size(pdc_beta)))) = 1;             % one diagonal
    % imagesc(abs(pdc_beta)); colorbar
    idx = f>30 & f<=50;
    pdc_lowgamma = tril(mean(pdc(:,:,idx),3,'omitnan'),-1); % zero upper triangle
    pdc_lowgamma(logical(eye(size(pdc_lowgamma)))) = 1;     % one diagonal
    % imagesc(abs(pdc_lowgamma)); colorbar
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
    % figure; imagesc(abs(gpdc_delta)); colorbar
    idx = f>=3 & f<8;
    gpdc_theta = tril(mean(gpdc(:,:,idx),3,'omitnan'),-1);    % zero upper triangle
    gpdc_theta(logical(eye(size(gpdc_theta)))) = 1;           % one diagonal
    % figure; imagesc(abs(gpdc_theta)); colorbar
    idx = f>=8 & f<=13;
    gpdc_alpha = tril(mean(gpdc(:,:,idx),3,'omitnan'),-1);    % zero upper triangle
    gpdc_alpha(logical(eye(size(gpdc_alpha)))) = 1;           % one diagonal
    % figure; imagesc(abs(gpdc_alpha)); colorbar
    % labels = {chanlocs.labels};
    % plot_corrmatrix(gpdc_alpha,labels)
    idx = f>13 & f<=30;
    gpdc_beta = tril(mean(gpdc(:,:,idx),3,'omitnan'),-1);     % zero upper triangle
    gpdc_beta(logical(eye(size(gpdc_beta)))) = 1;             % one diagonal
    figure; imagesc(abs(gpdc_beta)); colorbar
    idx = f>30 & f<=50;
    gpdc_lowgamma = tril(mean(gpdc(:,:,idx),3,'omitnan'),-1); % zero upper triangle
    gpdc_lowgamma(logical(eye(size(gpdc_lowgamma)))) = 1;     % one diagonal
    % figure; imagesc(abs(gpdc_lowgamma)); colorbar
    eeg_features.frequency.eeg_gpdc_delta = gpdc_delta;
    eeg_features.frequency.eeg_gpdc_theta = gpdc_theta;
    eeg_features.frequency.eeg_gpdc_alpha = gpdc_alpha;
    eeg_features.frequency.eeg_gpdc_beta = gpdc_beta;
    eeg_features.frequency.eeg_gpdc_lowgamma = gpdc_lowgamma;

    warning('Causal relationships between neighboring channels should be interpreted with caution due to volume conduction effects. See for example https://escholarship.org/content/qt5fb0q5wx/qt5fb0q5wx.pdf')

end

%% Entropy

if params.eeg_nonlinear

    % default parameters
    m = 2;
    r = .15;
    n = 2;
    tau = 1;
    coarseType = 'Standard deviation';	% coarse graining method
    nScales = 30;						% number of scale factors to compute
    filtData = true;  					% bandpass filter each scale factor (see Kosciessa et al. 2020)

    % Initiate progressbar (only when not in parpool)
    disp('Calculating nonlinear-domain EEG features...')
    if ~params.parpool
        progressbar('Extracting entropy on each EEG channel')
    end
    for iChan = 1:nChan

        fprintf('Processing channel %g \n', iChan);

        % Downsample to accelerate on data >5 min with sample rate of 256 hz
        if size(signals,2) > 76800 
            new_fs = 90;  % for Nyquist freq = lowpass cutoff (i.e. 45 Hz)
            fac = fs / new_fs; % downsample factor

            % downsample if integer, otherwise decimate to round factor
            if fac ~= floor(fac)
                fac = round(fac);
                signals_res = decimate(signals(iChan,:), fac);
                fprintf('Decimating EEG data to %g Hz sample rate to compute entropy on these large datasets... \n',new_fs)
            else
                signals_res = resample(signals(iChan,:), 1, fac);
                fprintf('Downsampling EEG data to %g Hz sample rate to compute entropy on these large datasets... \n',new_fs)
            end
            % Plot to check
            % times_res = (0:1/new_fs:(length(signals_res(iChan,:))-1)/new_fs)*1000;
            % figure; plot(times(1:fs*5), signals(iChan,1:fs*5)); % plot 5 s of data
            % hold on; plot(times_res(1:new_fs*5), signals_res(iChan,1:new_fs*5));

            % Fuzzy entropy
            fprintf('Computing fuzzy entropy...')
            fe = compute_fe(signals_res, m, r, n, tau, params.gpu);

            % Multiscale fuzzy entropy
            fprintf('Computing multiscale fuzzy entropy...')
            [mfe, scales, scale_bounds] = compute_mfe(signals_res, m, r, tau, coarseType, nScales, filtData, new_fs, n, params.gpu);
            % plot(scales(end:-1:1),mfe(end:-1:1));  hold on;
            % title('downsampled'); axis tight; box on; grid on
            % xticks(scales); xticklabels(scale_bounds(end:-1:1)); xtickangle(45)

            % Refined composite multiscale fuzzy entropy (without filtering)
            % [rcmfe, scales] = compute_rcmfe(signals_res, m, r, tau, coarseType, nScales, new_fs, n, params.gpu);
			
			if ~params.parpool
				progressbar(iChan / nChan);
			end

        else

            % Fuzzy entropy
            disp('Computing fuzzy entropy...')
            fe = compute_fe(signals(iChan,:), m, r, n, tau, params.gpu);

            % Multiscale fuzzy entropy
            disp('Computing multiscale fuzzy entropy...')
            [mfe, scales, scale_bounds] = compute_mfe(signals(iChan,:), m, r, tau, coarseType, nScales, filtData, fs, n, params.gpu);
            plot(scales(end:-1:1),mfe(end:-1:1)); hold on; axis tight; box on; grid on
            xticks(scales); xticklabels(scale_bounds(end:-1:1)); xtickangle(45)

            % Refined composite multiscale fuzzy entropy (without filtering)
            % disp('Computing refined composite multiscale fuzzy entropy...')
            % [rcmfe, scales] = compute_rcmfe(signals(iChan,:), m, r, tau, coarseType, nScales, fs, n, params.gpu);
            % plot(scales(end:-1:1),rcmfe(end:-1:1)); hold on; axis tight; box on; grid on
            % xticks(scales); xticklabels(scale_bounds(end:-1:1)); xtickangle(45)
			
			if ~params.parpool
				progressbar(iChan / nChan);
			end
        end

        % Outputs
        eeg_features.nonlinear.FE(iChan,:) = fe;
        eeg_features.nonlinear.MFE_scales(iChan,:) = scales;
        eeg_features.nonlinear.MFE_scale_bounds(iChan,:) = scale_bounds;
        eeg_features.nonlinear.MFE(iChan,:) = mfe;
        eeg_features.nonlinear.MFE_mean(iChan,:) = mean(mfe);
        eeg_features.nonlinear.MFE_sd(iChan,:) = std(mfe);
        [~,eeg_features.nonlinear.MFE_peak(iChan,:)] = max(mfe);
        eeg_features.nonlinear.MFE_area(iChan,:) = trapz(mfe);

    end
end

fprintf('Time to extract EEG features: %g min \n', round(toc(tstart)/60,1))
