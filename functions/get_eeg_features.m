% GET_EEG_FEATURES - Extract EEG features in the time, frequency and nonlinear domains.
%
% Usage:
%   [eeg_features, params] = get_eeg_features(signals, params)
%
% Inputs:
%   signals - EEG data (channels x samples, uV)
%   params  - structure with fields:
%     fs, chanlocs   - sample rate (Hz) and EEGLAB channel locations (required)
%     eeg_time, eeg_frequency, eeg_nonlinear - compute each domain (true/false, required)
%     eeg_frange     - PSD frequency range in Hz (default [1 40])
%     eeg_wintype    - pwelch taper (default 'hamming')
%     eeg_winlen     - pwelch window length in s (default 2)
%     eeg_winoverlap - window overlap in % (default 50)
%     eeg_freqbounds - 'conventional' (default) or 'individualized': alpha
%                      bounds from the alpha peak of the spectra (median
%                      across channels with a detectable peak, see
%                      get_freqBounds); theta then ends and beta starts at
%                      these bounds, the other limits are conventional.
%                      Falls back to 'conventional' if no alpha peak is found.
%     eeg_norm       - band power units: 0 = uV^2/Hz, 1 = dB (default),
%                      2 = uV^2/Hz divided by the channel's total power
%     asy_norm       - normalize alpha asymmetry (default false, see compute_asymmetry)
%     gpu, parpool   - use GPU / parallel pool (default false). Without
%                      parpool, parfor loops run serially.
%
% Outputs:
%   eeg_features - structure (one row per channel):
%     time      - rms, mode, var, skewness, kurtosis, iqr
%     frequency - freqs, pwr (PSD, channels x freqs; in dB if eeg_norm > 0),
%                 mean band power: delta (f(1)-3 Hz), theta (4-7), alpha (8-13),
%                 beta (13-30), gamma (30-fRange(2)) with conventional bounds;
%                 bands (band limits used, Hz, rows delta to gamma); IAF and
%                 IAF_mean (alpha center of gravity, Hz); alpha asymmetry
%                 and electrode pairs
%     qeeg      - absolute (uV^2) and relative band power, total power,
%                 alpha/theta, theta/beta, alpha/beta, alpha/(theta+beta),
%                 peak alpha frequency (7-13 Hz), median and 90% spectral edge frequency (Hz)
%     nonlinear - FD (fractal dimension), FE (fuzzy entropy, m=2, r=0.15, n=2,
%                 tau=1), computed on data resampled to ~90 Hz when fs > 100 Hz
%   params       - input params plus the defaults used (eeg_frange, eeg_wintype,
%                  eeg_winlen, eeg_winoverlap, eeg_freqbounds, eeg_norm, asy_norm)
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
% parfor(..., 0) runs serially in the client without starting a pool
if useparpool, maxWorkers = Inf; else, maxWorkers = 0; end
nChan = size(signals,1);

% Frequency domain parameters
if params.eeg_frequency
    if isfield(params,'eeg_frange') && ~isempty(params.eeg_frange)
        fRange = params.eeg_frange;
    else
        fRange = [1 40];        % overall frequency range to compute PSD (in Hz)
        params.eeg_frange = fRange; % to export
    end
    if isfield(params,'eeg_wintype') && ~isempty(params.eeg_wintype)
        wintype = params.eeg_wintype;
    else
        wintype = 'hamming';    % window type. Default = 'hamming' (see Smith et al, 2017 for asymmetry)
        params.eeg_wintype = wintype; % to export
    end
    if isfield(params,'eeg_winlen') && ~isempty(params.eeg_winlen)
        winlen = params.eeg_winlen;
    else
        winlen = 2;            % window size (in s). Default = 2 (see Smith et al, 2017 for asymmetry)
        params.eeg_winlen = winlen; % to export
    end
    if isfield(params,'eeg_winoverlap') && ~isempty(params.eeg_winoverlap)
        overlap = params.eeg_winoverlap;
    else
        overlap = 50;           % window overlap. Default = 50% (see Smith et al, 2017 for asymmetry)
        params.eeg_winoverlap = overlap; % to export
    end
    if isfield(params,'eeg_freqbounds') && ~isempty(params.eeg_freqbounds)
        freqbounds = lower(params.eeg_freqbounds);
    else
        freqbounds = 'conventional';    % band limits for band power (see header)
        params.eeg_freqbounds = freqbounds; % to export
    end
    if isfield(params,'eeg_norm') && ~isempty(params.eeg_norm)
        eeg_norm = params.eeg_norm;
    else
        eeg_norm = 1;      % band power in uV^2/Hz (0), dB (1), or divided by the channel's total power (2)
        params.eeg_norm = eeg_norm; % to export
    end
    if isfield(params,'asy_norm') && ~isempty(params.asy_norm)
        asy_norm = params.asy_norm;
    else
        asy_norm = false;      % normalize alpha asymmetry (true) or not (false), see compute_asymmetry and Smith et al. (2017)
        params.asy_norm = asy_norm; % to export
    end

end

% Nonlinear domain parameters (fuzzy entropy)
if params.eeg_nonlinear
    m = 2;      % embedding dimension
    r = .15;    % tolerance (fraction of the signal's SD)
    n = 2;      % fuzzy power
    tau = 1;    % time lag
end

disp('----------------------------------------------------')
disp('               Extracting EEG features ')
disp('----------------------------------------------------')


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

    % Frequency vector, computed once outside the parfor loop so that all
    % iterations (and workers) share the same one
    [~, ~, f] = compute_psd(signals(1,:),fs*winlen,wintype,overlap,[],fs,fRange,'psd',usegpu);
    f = gather(f);

    % Band limits (Hz), rows: delta, theta, alpha, beta, gamma
    bands = [f(1) 3; 4 7; 8 13; 13 30; 30 fRange(2)];

    % Individualized bounds: alpha bounds (minima on either side of the alpha
    % peak, 7-14 Hz search window) of each channel with a detectable peak,
    % median across these channels (similar to the individual alpha window
    % of restingIAF). Theta ends and beta starts at these bounds. The other
    % limits stay conventional, as peaks are rarely detectable in the
    % delta, theta, beta and gamma bands of resting-state spectra.
    if strcmp(freqbounds,'individualized')
        disp('Estimating individualized alpha band bounds...')
        pwr_all = gather(compute_psd(signals,fs*winlen,wintype,overlap,[],fs,fRange,'psd',usegpu));
        alphaBounds = nan(nChan,2);
        for iChan = 1:nChan
            try
                alphaBounds(iChan,:) = get_freqBounds(pwr_all(iChan,:), f, fs, [7 14], fs*winlen, 1);
            catch
            end
        end
        alphaBounds(any(isnan(alphaBounds),2),:) = [];
        alphaBounds = median(alphaBounds,1);
        if ~isempty(alphaBounds) && all(~isnan(alphaBounds)) && alphaBounds(1) > bands(2,1) && alphaBounds(2) < bands(4,2)
            bands(2,2) = alphaBounds(1);
            bands(3,:) = alphaBounds;
            bands(4,1) = alphaBounds(2);
            fprintf('Individualized alpha band: %.2f-%.2f Hz \n', alphaBounds)
        else
            warning("No alpha peak detected to individualize the frequency bands. Using conventional bands.")
        end
    end

    % Preallocate qEEG outputs
    ABS_DELTA = nan(nChan,1); ABS_THETA = nan(nChan,1); ABS_ALPHA = nan(nChan,1); ABS_BETA = nan(nChan,1); ABS_GAMMA = nan(nChan,1);
    REL_DELTA = nan(nChan,1); REL_THETA = nan(nChan,1); REL_ALPHA = nan(nChan,1); REL_BETA = nan(nChan,1); REL_GAMMA = nan(nChan,1);
    TOT_PWR   = nan(nChan,1);
    
    R_AT = nan(nChan,1); R_TB = nan(nChan,1); R_AB = nan(nChan,1); R_A_TB = nan(nChan,1);
    
    IAF_Hz = nan(nChan,1); MF_Hz = nan(nChan,1); SEF90_Hz = nan(nChan,1);
    
    % Frequency resolution (Hz), to integrate the PSD into power
    df = mean(diff(f));

    % progressbar (only when not in parpool)
    if ~useparpool
        progressbar('Computing power spectral density for each EEG channel')
    end
    disp('Calculating band-power on each EEG channel:')
    parfor (iChan = 1:nChan, maxWorkers)

        fprintf('  - channel %g \n', iChan)

        if usegpu
            sig = gpuArray(signals(iChan,:));
        else
            sig = signals(iChan,:);
        end

        % Compute PSD using pwelch
        [pwr, pwr_dB, ~] = compute_psd(sig,fs*winlen,wintype,overlap,[],fs,fRange,'psd',usegpu);

        % Band masks (limits set above)
        idxD = f >= bands(1,1) & f <= bands(1,2);
        idxT = f >= bands(2,1) & f <= bands(2,2);
        idxA = f >= bands(3,1) & f <= bands(3,2);
        idxB = f >= bands(4,1) & f <= bands(4,2);
        idxG = f >= bands(5,1) & f <= bands(5,2);

        % Mean band power in the units set by eeg_norm
        bp = nan(1,5);
        if eeg_norm == 0       % no normalization (uV^2/Hz)
            bp = [mean(pwr(idxD)) mean(pwr(idxT)) mean(pwr(idxA)) mean(pwr(idxB)) mean(pwr(idxG))];
        elseif eeg_norm == 1    % dB
            bp = [mean(pwr_dB(idxD)) mean(pwr_dB(idxT)) mean(pwr_dB(idxA)) mean(pwr_dB(idxB)) mean(pwr_dB(idxG))];
        elseif eeg_norm == 2    % normalized by total power of same channel
            bp = [mean(pwr(idxD)) mean(pwr(idxT)) mean(pwr(idxA)) mean(pwr(idxB)) mean(pwr(idxG))] ./ sum(pwr);
        end

        PWR(iChan,:) = pwr;
        PWR_DB(iChan,:) = pwr_dB;
        DELTA(iChan,:) = bp(1);
        THETA(iChan,:) = bp(2);
        ALPHA(iChan,:) = bp(3);
        BETA(iChan,:) = bp(4);
        GAMMA(iChan,:) = bp(5);

        % qEEG features, all computed on linear power (PSD integrated over frequency)

        % Total power within fRange for this channel (uV^2)
        TOT_PWR(iChan) = sum(pwr) * df;

        % Absolute band power in linear units (µV^2)
        absD = sum(pwr(idxD)) * df;
        absT = sum(pwr(idxT)) * df;
        absA = sum(pwr(idxA)) * df;
        absB = sum(pwr(idxB)) * df;
        absG = sum(pwr(idxG)) * df;
        
        ABS_DELTA(iChan) = absD; ABS_THETA(iChan) = absT; ABS_ALPHA(iChan) = absA; ABS_BETA(iChan) = absB; ABS_GAMMA(iChan) = absG;
        
        % Relative band power (to total in fRange)
        den = max(TOT_PWR(iChan), eps);
        REL_DELTA(iChan) = absD / den;
        REL_THETA(iChan) = absT / den;
        REL_ALPHA(iChan) = absA / den;
        REL_BETA(iChan)  = absB / den;
        REL_GAMMA(iChan) = absG / den;
        
        % Ratios on linear power
        R_AT(iChan)   = absA / max(absT, eps);              % alpha/theta
        R_TB(iChan)   = absT / max(absB, eps);              % theta/beta
        R_AB(iChan)   = absA / max(absB, eps);              % alpha/beta
        R_A_TB(iChan) = absA / max(absT + absB, eps);       % alpha/(theta+beta)
        
        % Peak alpha frequency: frequency of maximum power between 7 and 13 Hz
        aSearch = (f >= 7 & f <= 13);
        if any(aSearch)
            fA = f(aSearch);
            pA = pwr(aSearch);
            [~, ix] = max(pA);
            IAF_Hz(iChan) = fA(ix);
        else
            IAF_Hz(iChan) = NaN;
        end
        
        % Median frequency (MF), the frequency that splits the power 
        % spectrum into two equal halves.
        % Spectral edge frequency 90 (SEF90), the frequency below which 90% 
        % of the total spectral power is contained.
        cs = cumsum(pwr) * df;   % df = mean(diff(f)) defined before the parfor
        if cs(end) > 0
            MF_Hz(iChan)    = interp1(cs, f, 0.5*cs(end), 'linear', 'extrap');
            SEF90_Hz(iChan) = interp1(cs, f, 0.9*cs(end), 'linear', 'extrap');
        else
            MF_Hz(iChan) = NaN;
            SEF90_Hz(iChan) = NaN;
        end


        if ~useparpool
            progressbar(iChan/nChan)
        end

    end

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
    eeg_features.frequency.bands = bands;

    % qEEG features (not rounded)
    eeg_features.qeeg.delta_abs   = ABS_DELTA;
    eeg_features.qeeg.theta_abs   = ABS_THETA;
    eeg_features.qeeg.alpha_abs   = ABS_ALPHA;
    eeg_features.qeeg.beta_abs    = ABS_BETA;
    eeg_features.qeeg.gamma_abs   = ABS_GAMMA;
    eeg_features.qeeg.total_abs   = TOT_PWR;
    
    eeg_features.qeeg.delta_rel   = REL_DELTA;
    eeg_features.qeeg.theta_rel   = REL_THETA;
    eeg_features.qeeg.alpha_rel   = REL_ALPHA;
    eeg_features.qeeg.beta_rel    = REL_BETA;
    eeg_features.qeeg.gamma_rel   = REL_GAMMA;
    
    eeg_features.qeeg.alpha_theta = R_AT;
    eeg_features.qeeg.theta_beta  = R_TB;
    eeg_features.qeeg.alpha_beta  = R_AB;
    eeg_features.qeeg.alpha_tplusb = R_A_TB;
    
    eeg_features.qeeg.iaf    = IAF_Hz;
    eeg_features.qeeg.median = MF_Hz;
    eeg_features.qeeg.sef90  = SEF90_Hz;


    %%%%% Individual alpha frequency (IAF) %%%%%
    % Alpha center of gravity (CoG) from restingIAF (Corcoran et al. 2018):
    % 1-30 Hz, alpha search window 7-14 Hz, Savitzky-Golay frame width 11 and
    % order 5, mean CoG requires at least 1 channel (restingIAF requires an
    % integer sample rate)
    disp('Attempting to find the individual alpha frequency (IAF) for each EEG channel...')
    [pSum, pChans, ~] = restingIAF(signals, size(signals,1), 1, [1 30], round(fs), [7 14], 11, 5);
    eeg_features.frequency.IAF_mean = round(pSum.cog,3);
    eeg_features.frequency.IAF = round([pChans.gravs]',3);
    if ~isnan(eeg_features.frequency.IAF_mean)
        fprintf('Mean IAF across all channels: %g \n', eeg_features.frequency.IAF_mean)
    elseif sum(isnan(eeg_features.frequency.IAF)) == length(chanlocs)
        warning("Failed to find the IAF on all EEG channels. This can be due to improperly preprocessed data or lack of alpha peak in the power spectral distribution.")
    end

    %%%%% Alpha asymmetry %%%%%
    if length(chanlocs)>1
        alpha_pwr = mean(PWR(:,f >= bands(3,1) & f <= bands(3,2)),2,'omitnan');  % IMPORTANT: use power in μV^2/Hz here, NOT in log or decibels
        tot_pwr = mean(PWR,2,'omitnan');    % mean PSD over fRange (same units), for asy_norm
        [asy, pairLabels, pairNums] = compute_asymmetry(alpha_pwr, asy_norm, chanlocs, false, tot_pwr);
        eeg_features.frequency.asymmetry = round(asy,3);
        eeg_features.frequency.asymmetry_pairs_labels = pairLabels;
        eeg_features.frequency.asymmetry_pairs_num = pairNums;
    else
        warning("Only one EEG channel detected. Cannot compute alpha asymmetry.")
    end
end

%% Entropy

if params.eeg_nonlinear

    % Initiate progressbar (only when not in parpool)
    disp('Computing EEG features in the nonlinear-domain (this may take a while)...')
    if ~useparpool
        progressbar('Computing nonlinear features on all EEG channels')
    end
    
    % Resample to ~90 Hz when fs > 100 Hz, to limit memory use and computation
    % time of the entropy measures. Done regardless of the recording length
    % so that entropy values are always computed at the same sample rate
    % (they depend on it) and remain comparable across files.
    if fs > 100
        new_fs = 90;        % Nyquist freq = default lowpass cutoff (i.e. 45 Hz)
        [p, q] = rat(new_fs/fs, 1e-4);  % resampling factors (resample applies an anti-aliasing filter)
        fs = fs*p/q;        % actual new sample rate
        fprintf('Resampling EEG data to %g Hz to avoid memory issues and increase speed... \n', round(fs,2))
        signals_res = nan(nChan, ceil(size(signals,2)*p/q));
        for iChan = 1:nChan
            signals_res(iChan,:) = resample(double(signals(iChan,:)), p, q);
        end
        signals = signals_res;
    end
    
    parfor (iChan = 1:nChan, maxWorkers)

        if usegpu
            sig = gpuArray(signals(iChan,:));
        else
            sig = signals(iChan,:);
        end

        fprintf(' channel %g... \n', iChan);

        % Fractal dimension (box counting)
        fd(iChan,:) = fractal_volatility(sig);

        % Fuzzy entropy
        fe(iChan,:) = compute_fe(sig, m, r, n, tau);

        if ~useparpool
            progressbar(iChan / nChan);
        end
    end

    % Outputs
    eeg_features.nonlinear.FD = fd;     % fractal dimension
    eeg_features.nonlinear.FE = fe;     % fuzzy entropy

end

fprintf('Time to extract EEG features: %g min \n', round(toc(tstart)/60,1))
