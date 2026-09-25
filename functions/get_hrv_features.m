% GET_HRV_FEATURES - Heart-rate variability features from an NN interval series.
%
% Usage:
%   [HRV, params] = get_hrv_features(NN, NN_times, params)
%
% Inputs:
%   NN       - NN intervals (s), e.g. from clean_rr
%   NN_times - time of each NN interval (s)
%   params   - struct. Required: .hrv_time, .hrv_frequency, .hrv_nonlinear
%              (true/false for each domain). Optional:
%              .hrv_spec    - PSD method: 'LombScargle_norm' (default),
%                             'LombScargle', 'welch' or 'fft' (the last two
%                             on NN resampled at 7 Hz, mean removed)
%              .hrv_norm    - divide band powers by total power (default false)
%              .hrv_overlap - sliding-window overlap, fraction (default 0.25)
%
% Outputs:
%   HRV    - struct with fields present for the domains computed:
%     .time      - SDNN, RMSSD (ms), pNN50 (%)
%     .frequency - ulf, vlf, lf, hf: band power (integral of the PSD over
%                  the band) averaged over windows. Units: ms^2 for
%                  'LombScargle', 'welch' and 'fft' (PSD of NN in s^2/Hz,
%                  x1e6); unitless for 'LombScargle_norm' (periodogram
%                  normalized by the variance, comparable across bands and
%                  recordings, not with ms^2); fraction of ttlpwr if
%                  hrv_norm. Bands the recording is too short for are
%                  absent; all four are NaN if none can be estimated (< 34 s).
%                  Also lfhf, ttlpwr (hrv_norm only), pwr/pwr_freqs
%                  (spectrum of the last window of each band, concatenated;
%                  s^2/Hz or unitless), bands (Hz)
%     .nonlinear - Poincare SD1, SD2 (ms) and SD1SD2, FE (fuzzy entropy),
%                  FD (fractal volatility), PRSA_AC, PRSA_DC (ms)
%   params - input params plus the settings used, for export (hrv_spec,
%            hrv_norm, hrv_overlap, hrv_band_freqs, hrv_band_names, ...)
%
% Notes:
%   Bands: ULF 0-0.003, VLF 0.003-0.04, LF 0.04-0.15, HF 0.15-0.40 Hz. Each
%   band is estimated in sliding windows of 5 cycles of its lowest
%   frequency (Task Force, 1996): 24 h (ULF), ~28 min (VLF), 125 s (LF),
%   34 s (HF). A band is skipped, with a warning, if the recording is
%   shorter. Lomb-Scargle is the default because it handles the uneven
%   sampling and gaps of NN series without resampling. Time and nonlinear
%   features use the whole series.
%
% Copyright (C) - Cedric Cannard, 2023

function [HRV, params] = get_hrv_features(NN, NN_times, params)

%% Time domain

if params.hrv_time

    disp('Extracting HRV features in the time domain...')

    % SDNN (standard deviation of the NN intervals)
    HRV.time.SDNN = round(std(NN.*1000),1);   % in ms

    % RMSSD (root mean square of successive NN differences)
    HRV.time.RMSSD = round(sqrt(mean(diff(NN.*1000).^2)),1);  % in ms

    % pNN50 (percentage of successive differences >= alpha = 50 ms);
    % requires at least 2 min of data (Shaffer & Ginsberg, 2017)
    alpha = 50;
    HRV.time.pNN50 = round(100 * sum( abs(diff(NN)) >= alpha/1000 )/length(diff(NN)),1);  % in %

end

%% Frequency domain

if params.hrv_frequency

    disp('Extracting HRV features in the frequency domain...')

    % Parameters (the settings used are copied into params for export)
    if isfield(params,'hrv_spec') && ~isempty(params.hrv_spec)
        hrv_spec = params.hrv_spec;
    else
        hrv_spec = 'LombScargle_norm';
    end
    if ~any(strcmp(hrv_spec, {'LombScargle_norm' 'LombScargle' 'welch' 'fft'}))
        error("get_hrv_features: unknown hrv_spec '%s'. Use 'LombScargle_norm', 'LombScargle', 'welch' or 'fft'.", hrv_spec)
    end
    if isfield(params,'hrv_norm') && ~isempty(params.hrv_norm)
        norm = params.hrv_norm;
    else
        norm = false;
    end
    if isfield(params,'hrv_overlap') && ~isempty(params.hrv_overlap)
        overlap = params.hrv_overlap;
    else
        overlap = .25;  % window overlap (default = 25 %)
    end
    params.hrv_spec = hrv_spec;      % for exportation for users
    params.hrv_norm = norm;
    params.hrv_overlap = overlap;

    % Band power units: the normalized Lomb-Scargle periodogram is unitless;
    % the other PSDs are in s^2/Hz (NN in s), so band powers x1e6 are in ms^2
    if strcmp(hrv_spec, 'LombScargle_norm')
        pwrScale = 1;
    else
        pwrScale = 1e6;
    end
    HRV.frequency = struct();
    PWR = {}; PWR_freqs = {};

    % HRV frequency bands (ULF; VLF; LF; HF)
    bands = [ 0 .003; 0.003 .04; .04 .15; 0.15 0.40 ];
    bandNames = {'ULF' 'VLF' 'LF' 'HF'};
    params.hrv_band_freqs = bands;  % for exportation for users
    params.hrv_band_names = bandNames;  % for exportation for users

    % Minimum data length (s) for each band, also used as window length:
    % 24 h for ULF, then 5 cycles of the lowest frequency of VLF, LF and HF
    minLength = ceil([ 86400  5/0.003 5/0.04 5/0.15 ]);

    for iBand = 1:size(bands,1)

        if NN_times(end) < minLength(iBand)
            warning('File length is too short for estimating %s power reliably. At least %1.1f minutes are required. Cannot export this variable.', bandNames{iBand},minLength(iBand)/60)
        end

        % Sliding windows (none if the recording is shorter than winLength)
        winLength = minLength(iBand);  % in sec
        stepSize = floor(winLength * (1 - overlap));
        nWindows = floor((NN_times(end) - winLength) / stepSize) + 1;

        fprintf('HRV frequency band: %s \n', bandNames{iBand})

        % Compute PSD on each sliding window
        for iWin = 1:nWindows
            fprintf(' - window %g \n', iWin)

            % window bounds in s
            start_idx = (iWin - 1) * stepSize + 1;
            end_idx = start_idx + winLength - 1;
            win_idx = NN_times>=start_idx & NN_times<=end_idx;

            % warn if the NN intervals in this window cover less than 85% of it
            % (a gap left by the cleaning of RR artifacts by clean_rr)
            if sum(NN(win_idx)) < .85*winLength
                warning("This window contains a gap greater than 15% of the minimum window required for this band (likely from cleaning of RR artifacts by clean_rr).")
            end

            % Frequency vector for Lomb-Scargle: step 1/nfft Hz, nfft = next
            % power of 2 of the number of NN intervals in the window
            nfft = 2^nextpow2(length(NN(win_idx)));
            fvec = bands(iBand,1):1/nfft:bands(iBand,2);

            % Lomb-Scargle periodogram (no resampling required)
            if strcmp(hrv_spec, 'LombScargle_norm')
                [pwr,freqs] = plomb(NN(win_idx),NN_times(win_idx),fvec,'normalized');
                fprintf('Computing normalized Lomb-Scargle periodogram on the NN series... \n')
            elseif strcmp(hrv_spec, 'LombScargle')
                [pwr,freqs] = plomb(NN(win_idx),NN_times(win_idx),fvec,'psd');
                fprintf('Computing standard Lomb-Scargle periodogram on the NN series... \n')

            % Welch or FFT (require resampling)
            else
                % Resample at 7 Hz (cubic spline) and remove the mean, whose
                % leakage would otherwise inflate the band powers
                resamp_freq = 7;
                NN_resamp = resample_NN(NN_times(win_idx),NN(win_idx),resamp_freq,'cub');
                NN_resamp = NN_resamp(:) - mean(NN_resamp);

                % Pwelch (one segment spanning the window), one-sided PSD in s^2/Hz
                if strcmp(hrv_spec, 'welch')
                    welchWin = min(length(NN_resamp), round(minLength(iBand)*resamp_freq));  % window in samples, not s
                    [pwr,freqs] = pwelch(NN_resamp,welchWin,[],[],resamp_freq);
                    fprintf('Computing pwelch on the NN series... \n')

                % FFT (periodogram), one-sided PSD in s^2/Hz
                elseif strcmp(hrv_spec, 'fft')
                    nSamp = length(NN_resamp);
                    nHalf = floor(nSamp/2) + 1;             % 0 Hz to Nyquist
                    pwr = abs(fft(NN_resamp)).^2 / (resamp_freq*nSamp);
                    pwr = pwr(1:nHalf);
                    pwr(2:end-1+mod(nSamp,2)) = 2*pwr(2:end-1+mod(nSamp,2));  % fold negative frequencies (not 0 Hz or Nyquist)
                    freqs = (0:nHalf-1)' * resamp_freq / nSamp;
                    fprintf('Computing FFT on the NN series... \n')
                end
            end

            % Freq index
            freq_idx = bands(iBand,1) <= freqs & freqs <= bands(iBand,2);
            freq_res = freqs(2)-freqs(1); % resolution

            % Band power: integral of the PSD over the band (ms^2, or
            % unitless for the normalized Lomb-Scargle)
            bandPwr = sum(pwr(freq_idx)*freq_res) * pwrScale;
            if iBand == 1
                HRV.frequency.ulf(iWin,:) = bandPwr;      % ULF
            elseif iBand == 2
                HRV.frequency.vlf(iWin,:) = bandPwr;      % VLF
            elseif iBand == 3
                HRV.frequency.lf(iWin,:) = bandPwr;       % LF
            elseif iBand == 4
                HRV.frequency.hf(iWin,:) = bandPwr;       % HF
            end

            % Spectrum for export (overwritten each window: the last one is kept)
            PWR{iBand,1} = pwr;
            PWR_freqs{iBand,1} = freqs;

        end
    end

    % Recording shorter than the HF window: no band could be estimated
    if ~any(isfield(HRV.frequency, {'ulf' 'vlf' 'lf' 'hf'}))
        warning('Recording too short (%1.1f s) for HRV frequency features (at least %g s are needed for HF). Band powers set to NaN.', NN_times(end), minLength(end))
        [HRV.frequency.ulf, HRV.frequency.vlf, HRV.frequency.lf, HRV.frequency.hf] = deal(NaN);
    end

    % LF/HF ratio (of the power averaged across windows)
    if isfield(HRV.frequency,'lf') && isfield(HRV.frequency,'hf')
        HRV.frequency.lfhf = round(mean(HRV.frequency.lf) / mean(HRV.frequency.hf) * 100)/100;
    end

    % Total power (sum of the available band means, from ULF down to LF+HF)
    % and normalization of each band by it
    if norm
        try
            HRV.frequency.ttlpwr = sum([mean(HRV.frequency.ulf) mean(HRV.frequency.vlf) ...
                mean(HRV.frequency.lf) mean(HRV.frequency.hf)]);
        catch
            warndlg('HRV total power does not include ULF power (likely due to the short length of the cardiovascular time series).')
            warning('HRV total power does not include ULF power (likely due to the short length of the cardiovascular time series).')
            try
                HRV.frequency.ttlpwr = sum([mean(HRV.frequency.vlf) mean(HRV.frequency.lf) ...
                    mean(HRV.frequency.hf)]);
            catch
                warndlg('HRV total power does not include VLF power (likely due to the short length of the cardiovascular time series).')
                warning('HRV total power does not include VLF power (likely due to the short length of the cardiovascular time series).')
                try
                    HRV.frequency.ttlpwr = sum([mean(HRV.frequency.lf) mean(HRV.frequency.hf)]);
                catch
                    warndlg('Sorry, LF-HRV and HF-HRV power could not be normalized to total power.')
                    warning('Sorry, LF-HRV and HF-HRV power could not be normalized to total power.')
                end
            end
        end
        
        % Relative power: contribution of each band to the total
        if isfield(HRV.frequency,'ttlpwr')
            disp('Normalizing HRV power to overall power')
            if isfield(HRV.frequency,'ulf')
                HRV.frequency.ulf = mean(HRV.frequency.ulf) / HRV.frequency.ttlpwr;
            end
            if isfield(HRV.frequency,'vlf')
                HRV.frequency.vlf = mean(HRV.frequency.vlf) / HRV.frequency.ttlpwr;
            end
            if isfield(HRV.frequency,'lf')
                HRV.frequency.lf = mean(HRV.frequency.lf) / HRV.frequency.ttlpwr;
            end
            if isfield(HRV.frequency,'hf')
                HRV.frequency.hf = mean(HRV.frequency.hf) / HRV.frequency.ttlpwr;
            end
            % (the LF/HF ratio is already scale-free: not normalized)
        end
    end

    % remove empty cells (bands without any window)
    PWR(cellfun(@isempty,PWR)) = [];
    PWR_freqs(cellfun(@isempty,PWR_freqs)) = [];

    % Merge spectra from each band and export
    HRV.frequency.pwr_freqs = [cat(1, PWR_freqs{:})];
    HRV.frequency.pwr = [cat(1, PWR{:})];
    HRV.frequency.bands = bands;

    % Average across time windows (already scalars if normalized). ms^2 are
    % rounded to 2 decimals, unitless values to 4 significant digits
    if strcmp(hrv_spec, 'LombScargle_norm') || norm
        rnd = @(x) round(x, 4, 'significant');
    else
        rnd = @(x) round(x, 2);
    end
    if isfield(HRV.frequency,'ulf')
        HRV.frequency.ulf = rnd(mean(HRV.frequency.ulf,'omitnan'));
    end
    if isfield(HRV.frequency,'vlf')
        HRV.frequency.vlf = rnd(mean(HRV.frequency.vlf,'omitnan'));
    end
    if isfield(HRV.frequency,'lf')
        HRV.frequency.lf = rnd(mean(HRV.frequency.lf,'omitnan'));
    end
    if isfield(HRV.frequency,'hf')
        HRV.frequency.hf = rnd(mean(HRV.frequency.hf,'omitnan'));
    end
end

%% Nonlinear domain
if params.hrv_nonlinear

    disp('Extracting HRV features in the nonlinear domain (Poincare, fuzzy entropy, fractal dimension, PRSA)...')

    % Poincare plot descriptors
    SDSD = std(diff(NN));
    SDRR = std(NN);
    SD1 = (1 / sqrt(2)) * SDSD;     % measures the width of poincare cloud
    SD2 = sqrt((2 * SDRR^2) - (0.5 * SDSD^2));      % measures the length of the poincare cloud
    HRV.nonlinear.Poincare.SD1 = round(SD1*1000,3);      % in ms
    HRV.nonlinear.Poincare.SD2 = round(SD2*1000,3);      % in ms
    HRV.nonlinear.Poincare.SD1SD2 = round(SD1/SD2,3);

    % Fuzzy entropy parameters
    m = 2;      % embedding dimension
    r = .15;    % similarity bound (x SD)
    tau = 1;    % time lag
    n = 2;      % fuzzy power
    params.entropy_m = m;  % for exportation for users
    params.entropy_r = r;  % for exportation for users
    params.entropy_tau = tau;  % for exportation for users
    params.entropy_n = n;  % for exportation for users

    % Run Fuzzy entropy
    HRV.nonlinear.FE = compute_fe(NN, m, r, n, tau);

    % Fractal dimension (fractal volatility)
    HRV.nonlinear.FD = fractal_volatility(NN);

    % Phase rectified signal averaging (PRSA). Anchors: beats whose NN
    % shortens (AC) or lengthens (DC) by less than thresh % (larger changes
    % are treated as artifacts)
    fprintf('Computing phase rectified signal averaging (PRSA)... \n')
    thresh = 20;
    params.prsa_thresh = 20;  % for exportation for users
    lowAnchor = 1-thresh/100-0.0001; % lower limit for the AC anchor selection
    highAnchor = 1+thresh/100;      % The upper limit for the DC anchor selection
    NN = NN(:);
    drr_per = NN(2:end)./NN(1:end-1);   % drr_per(k) compares NN(k+1) with NN(k)
    ac_anchor = find((drr_per > lowAnchor) & (drr_per <= .9999)) + 1;  % shortening beats (anchor = NN(k+1))
    dc_anchor = find((drr_per > 1) & (drr_per <= highAnchor)) + 1;    % lengthening beats
    % PRSA (Bauer et al. 2006): average the NN segment around each anchor,
    % X(-2..1), then AC/DC = [X(0) + X(1) - X(-1) - X(-2)] / 4
    HRV.nonlinear.PRSA_AC = prsa_capacity(NN, ac_anchor);  % acceleration capacity (in ms, negative)
    HRV.nonlinear.PRSA_DC = prsa_capacity(NN, dc_anchor);  % deceleration capacity (in ms)

end

%% Subfunction
function cap = prsa_capacity(NN, anchors)
% Acceleration/deceleration capacity (ms) from phase-rectified signal
% averaging; anchors without two beats before and one after are dropped
anchors = anchors(anchors >= 3 & anchors <= length(NN)-1);
if isempty(anchors)
    cap = NaN; return
end
X = mean(NN(anchors + (-2:1)), 1);    % [X(-2) X(-1) X(0) X(1)]
cap = round(1000 * (X(3) + X(4) - X(2) - X(1)) / 4, 2);
