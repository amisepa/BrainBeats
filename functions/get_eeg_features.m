 %% Extract EEG features in time, fequency, and nonlinear domains.
%
% Example:
%   eeg_features = get_eeg_features(EEG.data,params)
% 
% Copyright (C) - Cedric Cannard, 2023

function [eeg_features, params] = get_eeg_features(signals,params)

tstart = tic;

% General parameters
if isfield(params,'fs') && ~isempty(params.fs)
    fs = params.fs;
else
    errordlg("The 'params' structure must contain EEG sample frequency. Example to add it: params.fs = 256;  % in Hz","error in get_eeg_features.m")
end
if isfield(params,'chanlocs') && ~isempty(params.chanlocs)
    chanlocs = params.chanlocs;
else
    errordlg("The 'params' structure must contain EEG channel locations. Example to add it: params.chanlocs = EEG.chanlocs;","error in get_eeg_features.m")
end
if isfield(params,'gpu') && ~isempty(params.gpu)
    usegpu = params.gpu;
else
    usegpu = false;
end
if isfield(params,'parpool') && ~isempty(params.parpool)
    useparpool = params.parpool;
else
    useparpool = false;
end

% Frequency domain parameters
if params.eeg_frequency
    if isfield(params,'eeg_frange') && ~isempty(params.eeg_frange)
        fRange = params.eeg_frange;
    else
        fRange = [1 40];        % overall frequency range to compute PSD (in Hz)
        params.fRange = fRange; % to export
    end
    if isfield(params,'eeg_wintype') && ~isempty(params.eeg_wintype)
        wintype = params.eeg_wintype;
    else
        wintype = 'hamming';    % window type. Default = 'hamming' (see Smith et al, 2017 for asymmetry)
        params.wintype = wintype; % to export
    end
    if isfield(params,'eeg_winlen') && ~isempty(params.eeg_winlen)
        winlen = params.eeg_winlen;
    else
        winlen = 2;            % window size (in s). Default = 2 (see Smith et al, 2017 for asymmetry)
        params.winlen = winlen; % to export
    end
    if isfield(params,'eeg_winoverlap') && ~isempty(params.eeg_winoverlap)
        overlap = params.eeg_winoverlap;
    else
        overlap = 50;           % window overlap. Default = 50% (see Smith et al, 2017 for asymmetry)
        params.overlap = overlap; % to export
    end
    if isfield(params,'eeg_freqbounds') && ~isempty(params.eeg_freqbounds)
        freqbounds = params.eeg_freqbounds;
    else
        freqbounds = 'conventional';      % freq bounds for band-power: 'conventional' (default) or 'individualized' (see Corcoran et al. 2018)
        params.freqbounds = freqbounds; % to export
    end
    if isfield(params,'eeg_norm') && ~isempty(params.eeg_norm)
        eeg_norm = params.eeg_norm;
    else
        eeg_norm = 1;      % none (0), normalize to decibels (1), normalize to decibels + divide by total power (2)
        params.eeg_norm = eeg_norm; % to export
    end
    if isfield(params,'asy_norm') && ~isempty(params.asy_norm)
        asy_norm = params.asy_norm;
    else
        asy_norm = false;      % normalization by dividing electrode's alpha power by total power (true) or not (false). see Smith et al. (2017)
        params.asy_norm = asy_norm; % to export
    end
    
end

% Nonlinear domain parameters
if params.eeg_nonlinear
    m = 2;
    r = .15;
    n = 2;
    tau = 1;
    % coarseType = 'Standard deviation';	% coarse graining method
    % nScales = 30;						    % number of scale factors to compute
    % filtData = true;  					% bandpass filter each scale factor (see Kosciessa et al. 2020)
end

disp('----------------------------------------------------')
disp('               Extracting EEG features ')
disp('----------------------------------------------------')


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

    % Get freqs once here to avoid parloop issue
    [~, ~, f] = compute_psd(signals(1,:),fs*winlen,wintype,overlap,[],fs,fRange,'psd',usegpu);

    % % number of frequency bins to preallocate memory
    % samplesPerWindow = fs * winlen;
    % nfft = 2^nextpow2(samplesPerWindow);
    % freqRes = fs / nfft;
    % % nFreqBins = ceil((fRange(2) - fRange(1)) / freqRes) + 1;
    % nFreqBins = floor((fRange(2) - fRange(1)) / freqRes);
    nFreqBins = length(f);

    % preallocate memory
    PWR = nan(nChan,nFreqBins);
    PWR_DB = nan(nChan,nFreqBins);
    DELTA = nan(nChan,1);
    THETA = nan(nChan,1);
    ALPHA = nan(nChan,1);
    BETA = nan(nChan,1);
    GAMMA = nan(nChan,1);

    % progressbar (only when not in parpool)
    if ~useparpool
        progressbar('Computing power spectral density for each EEG channel')
    end
    % tic
    disp('Calculating band-power on each EEG channel:')
    parfor iChan = 1:nChan

        delta = [];
        theta = [];
        alpha = [];
        beta = [];
        gamma = [];

        fprintf('  - channel %g \n', iChan)

        if usegpu
            sig = gpuArray(signals(iChan,:));
        else
            sig = signals(iChan,:);
        end

        % Compute PSD using pwelch
        [pwr, pwr_dB, ~] = compute_psd(sig,fs*winlen,wintype,overlap,[],fs,fRange,'psd',usegpu);

        % Delta
        if strcmp(freqbounds, 'conventional')
            if eeg_norm == 0
                delta = mean(pwr(f >= f(1) & f <= 3));       % no normalization (uV^2/Hz)
            elseif eeg_norm == 1
                delta = mean(pwr_dB(f >= f(1) & f <= 3));    % db
            elseif eeg_norm == 2
                delta = mean(pwr_dB(f >= f(1) & f <= 3)) ./ sum(pwr_dB);   % normalized by total power of same channel
            end
        elseif strcmp(freqbounds, 'individualized')
            try
                bounds = get_freqBounds(pwr, f, fs, [f(1) 3.5], winlen, 0.25);  % individualized frequency bounds
                delta = mean(pwr_dB(:,f >= bounds(1) & f <= bounds(2)));
            catch
                delta = NaN;
            end
        end

        % Theta
        if strcmp(freqbounds, 'conventional')
            if eeg_norm == 0
                theta = mean(pwr(f >= f(3) & f <= 7));       % no normalization (uV^2/Hz)
            elseif eeg_norm == 1
                theta= mean(pwr_dB(f >= f(3) & f <= 7));    % db
            elseif eeg_norm == 2
                theta = mean(pwr_dB(f >= f(3) & f <= 7)) ./ sum(pwr_dB);   % normalized by total power of same channel
            end
        elseif strcmp(freqbounds, 'individualized')
            try
                bounds = get_freqBounds(pwr, f, fs, [3 7], winlen, 0.25);  % individualized frequency bounds
                theta = mean(pwr_dB(iChan,f >= bounds(1) & f <= bounds(2)));
            catch
                theta = NaN;
            end
        end

        % Alpha
        if strcmp(freqbounds, 'conventional')
            if eeg_norm == 0
                alpha = mean(pwr(f >= f(8) & f <= 13));       % no normalization (uV^2/Hz)
            elseif eeg_norm == 1
                alpha = mean(pwr_dB(f >= f(8) & f <= 13));    % db
            elseif eeg_norm == 2
                alpha = mean(pwr_dB(f >= f(8) & f <= 13)) ./ sum(pwr_dB);   % normalized by total power of same channel
            end
        elseif strcmp(freqbounds, 'individualized')
            try
                bounds = get_freqBounds(pwr, f, fs, [7 14], winlen, 1);  % individualized frequency bounds
                alpha = mean(pwr_dB(f >= bounds(1) & f <= bounds(2)));
            catch
                alpha = NaN;
            end
        end

        % Beta
        if strcmp(freqbounds, 'conventional')
            if eeg_norm == 0
                beta = mean(pwr(:,f >= f(13) & f <= 30));       % no normalization (uV^2/Hz)
            elseif eeg_norm == 1
                beta = mean(pwr_dB(:,f >= f(13) & f <= 30));    % db
            elseif eeg_norm == 2
                beta = mean(pwr_dB(:,f >= f(13) & f <= 30)) ./ sum(pwr_dB);   % normalized by total power of same channel
            end
        elseif strcmp(freqbounds, 'individualized')
            try
                bounds = get_freqBounds(pwr, f, fs, [13 30], winlen, 0.25);  % individualized frequency bounds
                beta = mean(pwr_dB(iChan,f >= bounds(1) & f <= bounds(2)));
            catch
                beta = NaN;
            end
        end

        % Low gamma
        if strcmp(freqbounds, 'conventional')
            if eeg_norm == 0
                gamma = mean(pwr(:,f >= f(30) & f <= fRange(2)));       % no normalization (uV^2/Hz)
            elseif eeg_norm == 1
                gamma = mean(pwr_dB(:,f >= f(30) & f <= fRange(2)));    % db
            elseif eeg_norm == 2
                gamma = mean(pwr_dB(:,f >= f(30) & f <= fRange(2))) ./ sum(pwr_dB);   % normalized by total power of same channel
            end
        elseif strcmp(freqbounds, 'individualized')
            try
                bounds = get_freqBounds(pwr, f, fs, [30 fRange(2)], winlen, 0.25);  % individualized frequency bounds
                gamma = mean(pwr_dB(:,f >= bounds(1) & f <= bounds(2)));
            catch
                gamma = NaN;
            end
        end

        if ~useparpool
            progressbar(iChan/nChan)
        end
        
        PWR(iChan,:) = pwr;
        PWR_DB(iChan,:) = pwr_dB;
        DELTA(iChan,:) = delta;
        THETA(iChan,:) = theta;
        ALPHA(iChan,:) = alpha;
        BETA(iChan,:) = beta;
        GAMMA(iChan,:) = gamma;
        
    end
    % toc

    % Outputs
    eeg_features.frequency.freqs = f;
    if eeg_norm == 0
        eeg_features.frequency.pwr = PWR;
    elseif eeg_norm == 1 || eeg_norm == 2
        eeg_features.frequency.pwr = PWR_DB;
    end
    eeg_features.frequency.delta = round(DELTA,3);
    eeg_features.frequency.theta = round(THETA,3);
    eeg_features.frequency.alpha = round(ALPHA,3);
    eeg_features.frequency.beta = round(BETA,3);
    eeg_features.frequency.gamma = round(GAMMA,3);

    %%%%% Individual alpha frequency (IAF) %%%%%
    % Use alpha center of gravity (CoG) since it's the best
    disp('Attempting to find the individual alpha frequency (IAF) for each EEG channel...')
    [pSum, pChans, ~] = restingIAF(signals, size(signals,1), 1, [1 30], fs, [7 14], 11, 5);
    eeg_features.frequency.IAF_mean = round(pSum.cog,3);
    eeg_features.frequency.IAF = round([pChans.gravs]',3);
    if ~isnan(eeg_features.frequency.IAF_mean)
        fprintf('Mean IAF across all channels: %g \n', eeg_features.frequency.IAF_mean)
    elseif sum(isnan(eeg_features.frequency.IAF)) == length(chanlocs)
        warning("Failed to find the IAF on all EEG channels. This can be due to improperly preprocessed data or lack of alpha peak in the power spectral distribution.")
    end

    %%%%% Alpha asymmetry %%%%%
    alpha_pwr = mean(PWR(:,f >= 8 & f <= 13),2,'omitnan');  % IMPORTANT: use power in μV^2/Hz here, NOT in log or decibels
    [asy, pairLabels, pairNums] = compute_asymmetry(alpha_pwr, asy_norm, chanlocs, false);
    eeg_features.frequency.asymmetry = round(asy,3);
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
            fprintf('Decimating EEG data to a sample rate of %g Hz to avoid memory issues and increase speed... \n',new_fs)
        else
            fprintf('Downsampling EEG data to a sample rate of %g Hz to avoid memory issues and increase speed... \n',new_fs)
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
        % se(iChan,:) = compute_se_fast(sig,m,r);
        % se(iChan,:) = compute_se(sig,m,r,tau);

        % Fractal dimension
        fd(iChan,:) = fractal_volatility(sig);

        % Fuzzy entropy
        fe(iChan,:) = compute_fe(sig, m, r, n, tau, usegpu);

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
    % eeg_features.nonlinear.SE = se;     % sample entropy
    eeg_features.nonlinear.FD = fd;     % fractal dimension
    eeg_features.nonlinear.FE = fe;     % fuzzy entropy
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
