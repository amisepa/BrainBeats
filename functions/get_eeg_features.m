%% Extract EEG features in time, fequency, and nonlinear domains.
%
% Copyright (C) - Cedric Cannard, 2023

function eeg_features = get_eeg_features(signals,params)

tstart = tic;

% Default parameters: frequency domain
if params.eeg_frequency
    fRange = [1 40];    % FIXME: lowpass/highpass should be in params to make sure these are within filtered signal
    winSize = 2;        % window size (in s). Default = 2 (at least 2 s recommended by Smith et al, 2017 for asymmetry)
    winType = 'hamming';
    overlap = 50;       % 50% default (Smith et al. 2017)
end

% Default parameters: nonlinear domain
if params.eeg_nonlinear
    m = 2;
    r = .15;
    n = 2;
    tau = 1;
    % coarseType = 'Standard deviation';	% coarse graining method
    % nScales = 30;						    % number of scale factors to compute
    % filtData = true;  					% bandpass filter each scale factor (see Kosciessa et al. 2020)
end

% General parameters
fs = params.fs;
chanlocs = params.chanlocs;
usegpu = params.gpu;
useparpool = params.parpool;


disp('----------------------------------------------------')
disp('               Extracting EEG features ')
disp('----------------------------------------------------')


% Parallel computing
if useparpool
    p = gcp('nocreate');
    % delete(gcp('nocreate')) % shut down opened parpool
    if isempty(p) % if not already on, launch it
        disp('Initiating parrallel computing (all cores and threads -1)...')
        c = parcluster; % cluster profile
        % N = feature('numcores');        % physical number of cores
        N = getenv('NUMBER_OF_PROCESSORS'); % all processors (including threads)
        if ischar(N), N = str2double(N); end
        c.NumWorkers = N-2;  % update cluster profile to include all workers
        c.parpool();
    end
end

% Use parallel GPUs computing (if multiple GPUS are available)
% if useparpool && usegpu
%     availableGPUs = gpuDeviceCount("available");
%     if availableGPUs > 1
%         parpool('Processes',availableGPUs);
%         fprintf('%g GPUs detected. Using them in parallel pool. \n',availableGPUs)
%     else
%         fprintf('Only one GPU detected. Using normal GPU and parallel pool computing. \n')
%     end
% end


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

    % progressbar (only when not in parpool)
    if ~useparpool
        progressbar('Estimating EEG power spectral density on each channel')
    end

    %%%%% Band power %%%%%

    % number of frequency bins to preallocate memory
    samplesPerWindow = fs * winSize;
    nfft = 2^nextpow2(samplesPerWindow);
    freqResolution = fs / nfft;
    numFrequencyBins = floor((fRange(2) - fRange(1)) / freqResolution) + 1;
    
    % preallocate memory
    pwr = nan(nChan,numFrequencyBins);
    pwr_dB = nan(nChan,numFrequencyBins);
    delta = nan(nChan,1);
    delta_norm = nan(nChan,1);
    delta_indiv = nan(nChan,1);
    theta = nan(nChan,1);
    theta_norm = nan(nChan,1);
    theta_indiv = nan(nChan,1);
    alpha = nan(nChan,1);
    alpha_norm = nan(nChan,1);
    alpha_indiv = nan(nChan,1);
    beta = nan(nChan,1);
    beta_norm = nan(nChan,1);
    beta_indiv = nan(nChan,1);
    low_gamma = nan(nChan,1);
    low_gamma_norm = nan(nChan,1);
    low_gamma_indiv = nan(nChan,1);

    disp('Calculating band-power on each EEG channel:')
    for iChan = 1:nChan

        fprintf('  - channel %g \n', iChan)

        % Compute PSD using pwelch
        [pwr(iChan,:), pwr_dB(iChan,:), f] = compute_psd(signals(iChan,:),fs*winSize,winType,overlap,[],fs,fRange,'psd',usegpu);

        % Delta
        delta(iChan,:) = mean(pwr_dB(iChan,f >= f(1) & f <= 3));
        delta_norm(iChan,:) = delta(iChan,:) ./ sum(pwr_dB(iChan,:)); % normalized by total power of same channel
        try
            bounds = get_freqBounds(pwr(iChan,:), f, fs, [f(1) 3], winSize,0.25);  % individualized frequency bounds
            delta_indiv(iChan,:) = mean(pwr_dB(iChan,f >= bounds(1) & f <= bounds(2)));
        catch
            % fprintf('Could not detect individualized frequency bounds for Delta \n')
            delta_indiv(iChan,:) = NaN;
        end

        % Theta
        theta(iChan,:) = mean(pwr_dB(iChan,f >= 3 & f <= 7));
        theta_norm(iChan,:) = theta(iChan,:) ./ sum(pwr_dB(iChan,:)); % normalized by total power of same channel
        try
            bounds = get_freqBounds(pwr(iChan,:), f, fs, [2 8], winSize, 0.25);  % individualized lower/upper bounds
            theta_indiv(iChan,:) = mean(pwr_dB(iChan,f >= bounds(1) & f <= bounds(2)));
        catch
            % fprintf('Could not detect individualized frequency bounds for Theta \n')
            theta_indiv(iChan,:) = NaN;
        end

        % Alpha
        alpha(iChan,:) = mean(pwr_dB(iChan,f >= 7.5 & f <= 13));
        alpha_norm(iChan,:) = alpha(iChan,:) ./ sum(pwr_dB(iChan,:)); % normalized by total power of same channel
        try
            bounds = get_freqBounds(pwr(iChan,:), f, fs, [7 14], winSize, 1);  % individualized lower/upper bounds
            alpha_indiv(iChan,:) = mean(pwr_dB(iChan,f >= bounds(1) & f <= bounds(2)));
        catch
            % fprintf('Could not detect individualized frequency bounds for Alpha \n')
            alpha_indiv(iChan,:) = NaN;
        end

        % Beta
        beta(iChan,:) = mean(pwr_dB(iChan,f >= 13.5 & f <= 30));
        beta_norm(iChan,:) = beta(iChan,:) ./ sum(pwr_dB(iChan,:)); % normalized by total power of same channel
        try
            bounds = get_freqBounds(pwr(iChan,:), f, fs, [13 30], winSize, 0.25);  % individualized lower/upper bounds
            beta_indiv(iChan,:) = mean(pwr_dB(iChan,f >= bounds(1) & f <= bounds(2)));
        catch
            % fprintf('Could not detect individualized frequency bounds for Beta \n')
            beta_indiv(iChan,:) = NaN;
        end

        % Low gamma
        low_gamma(iChan,:) = mean(pwr_dB(iChan,f >= 31 & f <= fRange(2)));
        low_gamma_norm(iChan,:) = low_gamma(iChan,:) ./ sum(pwr_dB(iChan,:)); % normalized by total power of same channel
        try
            bounds = get_freqBounds(pwr(iChan,:), f, fs, [30 45], winSize, 0.25);  % individualized lower/upper bounds
            low_gamma_indiv(iChan,:) = mean(pwr_dB(iChan,f >= bounds(1) & f <= bounds(2)));
        catch
            % fprintf('Could not detect individualized frequency bounds for Low gamma \n')
            low_gamma_indiv(iChan,:) = NaN;
        end

        if ~useparpool
            progressbar(iChan/nChan)
        end
        
    end

    % Outputs
    eeg_features.frequency.freqs = f;
    eeg_features.frequency.pwr = pwr;
    eeg_features.frequency.pwr_dB = pwr_dB;
    eeg_features.frequency.delta = delta;
    eeg_features.frequency.delta_norm = delta_norm;
    eeg_features.frequency.delta_indiv = delta_indiv;
    eeg_features.frequency.theta = theta;
    eeg_features.frequency.theta_norm = theta_norm;
    eeg_features.frequency.theta_indiv = theta_indiv;
    eeg_features.frequency.alpha = alpha;
    eeg_features.frequency.alpha_norm = alpha_norm;
    eeg_features.frequency.alpha_indiv = alpha_indiv;
    eeg_features.frequency.beta = beta;
    eeg_features.frequency.beta_norm = beta_norm;
    eeg_features.frequency.beta_indiv = beta_indiv;
    eeg_features.frequency.low_gamma = low_gamma;
    eeg_features.frequency.low_gamma_norm = low_gamma_norm;
    eeg_features.frequency.low_gamma_indiv = low_gamma_indiv;

    %%%%% Individual alpha frequency (IAF) %%%%%
    % Use alpha center of gravity (CoG) since it's the best
    disp('Attempting to find the individual alpha frequency (IAF) for each EEG channel...')
    [pSum, pChans, ~] = restingIAF(signals, size(signals,1), 1, [1 30], fs, [7 13], 11, 5);
    eeg_features.frequency.IAF_mean = pSum.cog;
    eeg_features.frequency.IAF = [pChans.gravs]';
    if ~isnan(eeg_features.frequency.IAF_mean)
        fprintf('Mean IAF across all channels: %g \n', eeg_features.frequency.IAF_mean)
    elseif sum(isnan(eeg_features.frequency.IAF)) == length(chanlocs)
        warning("Failed to find the IAF on all EEG channels. This can be due to improperly preprocessed data or lack of alpha peak in the power spectral distribution.")
    end

    %%%%% Alpha asymmetry %%%%%
    % on log(pwr) no pwr_dB - on all possible symmetric pairs of electrodes
    norm = false;  % normalize by dividing by total alpha power
    alpha_pwr = mean(pwr(:,f >= 7.5 & f <= 13),2);  % IMPORTANT: use power in μV2/Hz NOT in decibels
    [asy, pairLabels, pairNums] = compute_asymmetry(alpha_pwr, norm, chanlocs, false);
    eeg_features.frequency.asymmetry = asy;
    eeg_features.frequency.asymmetry_pairs_labels = pairLabels;
    eeg_features.frequency.asymmetry_pairs_num = pairNums;

end

%% Entropy

if params.eeg_nonlinear


    % Initiate progressbar (only when not in parpool)
    disp('Computing EEG features in the nonlinear-domain (this may take a while)...')
    if ~useparpool
        progressbar('Computing nonlinear features on all EEG channels')
    end
    
    % Downsample/decimate if data are >2 min long and > 100 Hz sample rate
    if fs>100 && size(signals,2)/fs/60 > 2
        
        new_fs = 90;        % for Nyquist freq = default lowpass cutoff (i.e. 45 Hz)
        fac = fs / new_fs;  % downsample factor
        if fac ~= floor(fac)
            fac = round(fac);
            fprintf('Decimating EEG data to a sample rate of %g Hz sample rate to avoid memory issues and increase speed... \n',new_fs)
        else
            fprintf('Downsampling EEG data to a sample rate of %g Hz sample rate to avoid memory issues and increase speed... \n',new_fs)
        end

        % Downsample if integer, otherwise decimate to round factor
        signals_res = nan(nChan,ceil(size(signals,2)/ceil(fac)));
        for iChan = 1:nChan
            if fac ~= floor(fac)
                fac = round(fac);
                signals_res(iChan,:) = decimate(signals(iChan,:), fac);
            else
                signals_res(iChan,:) = resample(signals(iChan,:), 1, fac);
            end
            % Plot to check
            % times_res = (0:1/new_fs:(length(signals(iChan,:))-1)/new_fs)*1000;
            % figure; plot(times(1:fs*5), signals(iChan,1:fs*5)); % plot 5 s of data
            % hold on; plot(times_res(1:new_fs*5), signals_res(iChan,1:new_fs*5));
        end
        signals = signals_res;
        fs = new_fs;
    end
    
    % tic
    parfor iChan = 1:nChan
        
        if usegpu
            sig = gpuArray(signals(iChan,:));
        else
            sig = signals(iChan,:);
        end

        fprintf(' channel %g... \n', iChan);

        % Sample entropy (fast method)
        se(iChan,:) = compute_se_fast(sig,m,r);
        % se(iChan,:) = compute_se(sig,m,r,tau);

        % Fractal dimension
        fd(iChan,:) = fractal_volatility(sig);

        % Fuzzy entropy
        % fe(iChan,:) = compute_fe(sig, m, r, n, tau, usegpu);

        % Multiscale fuzzy entropy
        % disp('Computing multiscale fuzzy entropy...')
        % [mfe, scales, scale_bounds] = compute_mfe(sig, m, r, tau, coarseType, nScales, filtData, fs, n, usegpu);
        % plot(scales(end:-1:1),mfe(end:-1:1)); hold on; axis tight; box on; grid on
        % xticks(scales); xticklabels(scale_bounds(end:-1:1)); xtickangle(45)

        % Refined composite multiscale fuzzy entropy (without filtering)
        % disp('Computing refined composite multiscale fuzzy entropy...')
        % [rcmfe, scales] = compute_rcmfe(sig, m, r, tau, coarseType, nScales, fs, n, usegpu);
        % plot(scales(end:-1:1),rcmfe(end:-1:1)); hold on; axis tight; box on; grid on
        % xticks(scales); xticklabels(scale_bounds(end:-1:1)); xtickangle(45)

        if ~useparpool
            progressbar(iChan / nChan);
        end
    end
    % toc

    % Outputs
    eeg_features.nonlinear.SE = se;     % sample entropy
    eeg_features.nonlinear.FD = fd;     % fractal dimension
    % eeg_features.nonlinear.FE = fe;     % fuzzy entropy
    % eeg_features.nonlinear.MFE_scales(iChan,:) = scales;
    % eeg_features.nonlinear.MFE_scale_bounds(iChan,:) = scale_bounds;
    % eeg_features.nonlinear.MFE(iChan,:) = mfe;
    % eeg_features.nonlinear.MFE_mean(iChan,:) = mean(mfe);
    % eeg_features.nonlinear.MFE_sd(iChan,:) = std(mfe);
    % [~,eeg_features.nonlinear.MFE_peak(iChan,:)] = max(mfe);
    % eeg_features.nonlinear.MFE_area(iChan,:) = trapz(mfe);

    % figure('color','w')
    % subplot(2,2,1)
    % plot_topo(se,params.chanlocs,1,'entropy');
    % title('sample entropy')
    % subplot(2,2,2)
    % plot_topo(fd,params.chanlocs,1,'entropy');
    % title('fractal dimension')
    % subplot(2,2,3)
    % plot_topo(fe,params.chanlocs,1,'entropy');
    % title('fuzzy entropy')

end

fprintf('Time to extract EEG features: %g min \n', round(toc(tstart)/60,1))
