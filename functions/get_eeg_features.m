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
    eeg_features.frequency.IAF = [pChans.gravs];

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

        % compute asymmetry and export
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
    
    % Use this method instead? 
    % MVAR coefficients
    % Faes and Nollo (2011). Multivariate Frequency Domain Analysis of 
    % Causal Interactions in Physiological Time Series. https://www.intechopen.com/chapters/12918
    nfft = fs*2;  % 2-s windows
    fprintf('Computing multivariate frequency domain analysis of causal interactions between all EEG channels...')
    Su = eye(nChan,nChan);
    [DC,DTF,PDC,~,~,COH,PCOH,~,~,~,~,f] = fdMVAR_5order(signals,Su,nfft,fs);

    % Deal with complex and negative values
    if ~isreal(DC)
        DC = real(DC);
        % DC(DC < 0) = 0;
    end
    if ~isreal(DTF)
        DTF = real(DTF);
        % DTF(DTF < 0) = 0;
    end
    if ~isreal(PDC)
        PDC = real(PDC);
        % PDC(PDC < 0) = 0;
    end
    if ~isreal(COH)
        COH = real(COH);
        % COH(COH < 0) = 0;
    end
    if ~isreal(PCOH)
        PCOH = real(PCOH);
        % PCOH(PCOH < 0) = 0;
    end
    
    eeg_features.frequency.eeg_coh = COH;
    eeg_features.frequency.eeg_coh_f = f;
    eeg_features.frequency.eeg_pcoh = PCOH;
    eeg_features.frequency.eeg_dc = DC;
    eeg_features.frequency.eeg_pdc = PDC;
    eeg_features.frequency.eeg_dtf = DTF;

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
