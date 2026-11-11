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
%              .hrv_spec    - PSD method: 'LombScargle' (default),
%                             'LombScargle_norm', 'welch' or 'fft' (the last
%                             two on NN resampled at 7 Hz, mean removed)
%              .hrv_norm    - divide band powers by total power (default false)
%              .hrv_overlap - minimum overlap of the sliding windows,
%                             fraction from 0 to < 1 (default 0.25; see Notes)
%              .prsa_thresh - PRSA anchors: NN changes larger than this
%                             (%) are left out (default 20; Bauer et al.,
%                             2006: 5)
%
% Outputs:
%   HRV    - struct with fields present for the domains computed:
%     .time      - SDNN, RMSSD (ms), pNN50 (% of successive differences > 50 ms)
%     .frequency - ulf, vlf, lf, hf: band power (integral of the PSD over
%                  the band) averaged over windows. Units: ms^2 for
%                  'LombScargle', 'welch' and 'fft' (PSD of NN in s^2/Hz,
%                  x1e6); unitless for 'LombScargle_norm' (periodogram
%                  divided by twice the variance of each window). LF and HF
%                  come from the same windows, so the normalized LF/HF is
%                  close to the ms^2 one, but normalized values grow with the
%                  number of beats per window, i.e. the heart rate: same
%                  modulation in ms, LF +41% at 90 vs 60 bpm, ms^2 -3%. With
%                  hrv_norm, fraction of ttlpwr, the sum of the bands
%                  estimated (LF + HF below ~28 min, i.e. normalized units;
%                  NaN without LF, below 120 s). Bands the recording is too
%                  short for are absent; all four are NaN if none can be
%                  estimated (< 60 s). Also lfhf, ttlpwr (hrv_norm only),
%                  pwr/pwr_freqs (spectrum of the last window of each window
%                  length over its bands, concatenated: ULF, VLF, then LF and
%                  HF; ms^2/Hz, or unitless), bands (Hz)
%     .nonlinear - Poincare SD1, SD2 (ms) and SD1SD2, FE (fuzzy entropy) and
%                  FD (box-counting fractal dimension) from the ASCENT plugin
%                  (see BB_FUZZY_FRACTAL; NaN without it), PRSA_AC, PRSA_DC (ms)
%   params - input params plus the settings used, for export (hrv_spec,
%            hrv_norm, hrv_overlap, hrv_band_freqs, hrv_band_names, ...)
%
% Notes:
%   Bands: ULF 0-0.003, VLF 0.003-0.04, LF 0.04-0.15, HF 0.15-0.40 Hz. LF
%   and HF come from the same sliding windows of 120 s (Task Force, 1996:
%   about 1 min is needed for HF and 2 min for LF), one spectrum per window.
%   VLF uses windows of 5 cycles of its lowest frequency (~28 min), ULF 24 h.
%   The windows are spread evenly from the first to the last beat, as many
%   as needed for an overlap of at least 'hrv_overlap', so that the whole
%   recording is used (a fixed 90-s step left out 3-35% of 6 recordings of
%   3-5 min). A band is skipped, with a warning, if the recording is shorter
%   than its window, except HF, which is taken from the whole recording when
%   it lasts 60-120 s (LF is then skipped).
%   Window length (tests/compare_hrv_windows.m, synthetic series with a
%   1/f-like spectrum): shorter windows let VLF power into LF (Lomb-Scargle
%   LF +24%, +13%, +10% with 60, 90 and 120 s; welch +28%, +7%, +2%);
%   180 s did little better (+8%, 0%) but needs 3 min for LF and leaves one
%   window on 3-min recordings. Until v1.6, LF (125-s windows) and HF (34-s
%   windows, 0.03-Hz resolution) came from different windows: LF power
%   leaked into HF (Lomb-Scargle HF +14% on the same series, +5% now; on 6
%   recordings the old HF was 2-24% higher), and normalized LF and HF had
%   different variances (normalized LF/HF 1.02 vs 1.67 in ms^2 on the same
%   series; now 1.79 vs 1.80).
%   Lomb-Scargle is the default because it handles the uneven sampling and
%   gaps of NN series without resampling. It has no taper: VLF power adds
%   to LF (+7% with a 0.025-Hz, 50-ms component; welch -1%). Time and
%   nonlinear features use the whole series.
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

    % pNN50 (percentage of successive differences > alpha = 50 ms; Task
    % Force, 1996); requires at least 2 min of data (Shaffer & Ginsberg,
    % 2017). Differences rounded to 1 us first: one of exactly 50 ms (common
    % at 500 or 1000 Hz) was counted or not by floating-point chance
    alpha = 50;
    dNN = round(abs(diff(NN))*1e6)/1e3;     % in ms
    HRV.time.pNN50 = round(100 * sum(dNN > alpha)/length(dNN),1);  % in %

end

%% Frequency domain

if params.hrv_frequency

    disp('Extracting HRV features in the frequency domain...')

    % Parameters (the settings used are copied into params for export)
    if isfield(params,'hrv_spec') && ~isempty(params.hrv_spec)
        hrv_spec = params.hrv_spec;
    else
        hrv_spec = 'LombScargle';
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
        overlap = .25;  % minimum window overlap (default = 25 %)
    end
    if ~isscalar(overlap) || ~isnumeric(overlap) || overlap < 0 || overlap >= 1
        error('get_hrv_features: hrv_overlap must be a fraction from 0 to less than 1 (e.g. 0.25 for 25%%), not %s.', mat2str(overlap))
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
    bandFields = lower(bandNames);
    params.hrv_band_freqs = bands;  % for exportation for users
    params.hrv_band_names = bandNames;  % for exportation for users

    % Window length (s) of each band: 24 h for ULF, 5 cycles of the lowest
    % VLF frequency, and the same 120-s windows for LF and HF (one spectrum
    % per window gives both). minLength: shortest recording for each band
    % (HF: 60 s, then from the whole recording)
    lfhfWin = 120;
    winLength = [ 86400 ceil(5/0.003) lfhfWin lfhfWin ];
    minLength = [ 86400 ceil(5/0.003) lfhfWin 60 ];

    % Recording span covered by the NN intervals: from the R-peak that
    % starts the first interval to the last R-peak
    NN = NN(:); NN_times = NN_times(:);
    t0 = NN_times(1) - NN(1);
    span = NN_times(end) - t0;
    canDo = span >= minLength;
    shortHF = canDo & span < winLength;
    winLength(shortHF) = span;
    for iBand = find(~canDo)
        warning('File length is too short for estimating %s power reliably. At least %1.1f minutes are required. Cannot export this variable.', bandNames{iBand},minLength(iBand)/60)
    end
    if any(shortHF)
        fprintf('HF power from the whole recording (%1.1f s, shorter than the %g-s windows). \n', span, lfhfWin)
    end

    % Bands with the same window length are estimated from the same windows
    winLengths = unique(winLength(canDo), 'stable');
    for iLen = 1:numel(winLengths)

        winLen = winLengths(iLen);  % in s
        iBands = find(canDo & winLength == winLen);
        fRange = [ bands(iBands(1),1) bands(iBands(end),2) ];
        fprintf('HRV frequency band(s): %s (%g-s windows) \n', strjoin(bandNames(iBands),', '), round(winLen,1))

        % Sliding windows spread evenly from the first to the last beat: as
        % many as needed for an overlap of at least 'hrv_overlap'
        nWindows = ceil((span - winLen) / (winLen*(1 - overlap)) - 1e-9) + 1;
        if nWindows > 1
            winStart = t0 + (0:nWindows-1) * (span - winLen) / (nWindows - 1);
        else
            winStart = t0;
        end

        % Compute PSD on each sliding window
        for iWin = 1:nWindows
            fprintf(' - window %g \n', iWin)

            % NN intervals ending in the window
            win_idx = NN_times > winStart(iWin) & NN_times <= winStart(iWin) + winLen;

            % warn if the NN intervals in this window cover less than 85% of it
            % (a gap left by the cleaning of RR artifacts by clean_rr)
            if sum(NN(win_idx)) < .85*winLen
                warning("This window contains a gap greater than 15% of its length (likely from cleaning of RR artifacts by clean_rr).")
            end

            % One spectrum per window (s^2/Hz, or unitless)
            [pwr, freqs, freq_res] = nn_spectrum(NN(win_idx), NN_times(win_idx), winLen, fRange, hrv_spec);

            % Band power: integral of the PSD over the band (ms^2, or
            % unitless for the normalized Lomb-Scargle). Bands sharing the
            % spectrum are half-open, so that a frequency at the border is
            % counted once; the highest one includes its upper limit
            for iBand = iBands
                if iBand == iBands(end)
                    freq_idx = bands(iBand,1) <= freqs & freqs <= bands(iBand,2);
                else
                    freq_idx = bands(iBand,1) <= freqs & freqs < bands(iBand,2);
                end
                HRV.frequency.(bandFields{iBand})(iWin,1) = sum(pwr(freq_idx)*freq_res) * pwrScale;
            end
        end

        % Spectrum for export: last window, over the bands of this window
        % length, in ms^2/Hz or unitless
        keep = fRange(1) <= freqs & freqs <= fRange(2);
        PWR{end+1,1} = pwr(keep) * pwrScale; %#ok<AGROW>
        PWR_freqs{end+1,1} = freqs(keep); %#ok<AGROW>
    end

    % Recording shorter than the shortest HF window: no band could be estimated
    if ~any(isfield(HRV.frequency, {'ulf' 'vlf' 'lf' 'hf'}))
        warning('Recording too short (%1.1f s) for HRV frequency features (at least %g s are needed for HF). Band powers set to NaN.', span, minLength(end))
        [HRV.frequency.ulf, HRV.frequency.vlf, HRV.frequency.lf, HRV.frequency.hf] = deal(NaN);
    end

    % LF/HF ratio (of the power averaged across windows; the same windows
    % for both)
    if isfield(HRV.frequency,'lf') && isfield(HRV.frequency,'hf')
        HRV.frequency.lfhf = round(mean(HRV.frequency.lf) / mean(HRV.frequency.hf) * 100)/100;
    end

    % Total power (sum of the available band means: LF + HF, plus VLF and
    % ULF when the recording is long enough for them) and normalization of
    % each band by it. ULF needs 24 h and VLF ~28 min, so on most recordings
    % the total is LF + HF (LF and HF in normalized units). Without LF (60 to
    % 120 s), HF cannot be normalized: NaN, not ms^2 under a relative field.
    if norm
        have = bandFields(isfield(HRV.frequency, bandFields));
        if all(ismember({'lf' 'hf'}, have))
            HRV.frequency.ttlpwr = sum(cellfun(@(b) mean(HRV.frequency.(b)), have));
            fprintf('HRV band powers normalized by their total (%s). \n', strjoin(upper(have), ' + '))
            for b = have
                HRV.frequency.(b{1}) = mean(HRV.frequency.(b{1})) / HRV.frequency.ttlpwr;
            end
            % (the LF/HF ratio is already scale-free: not normalized)
        else
            warning('HRV band powers cannot be normalized without LF power (recording shorter than %g s): set to NaN.', lfhfWin)
            for b = have
                HRV.frequency.(b{1}) = NaN;
            end
            HRV.frequency.ttlpwr = NaN;
        end
    end

    % Merge spectra from each window length and export
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

    % Fuzzy entropy and fractal dimension (box counting), from the ASCENT plugin
    [HRV.nonlinear.FE, HRV.nonlinear.FD, params.ascent_version] = bb_fuzzy_fractal(NN(:)', m, r, n, tau);

    % Phase rectified signal averaging (PRSA). Anchors: beats whose NN
    % shortens (AC) or lengthens (DC) by less than thresh % (larger changes
    % are treated as artifacts). 20% by default: Bauer et al. (2006) used 5%
    % on 24-h recordings, but on 6 ds003838 resting recordings (NN already
    % cleaned by clean_rr) 5% left out 18-61% of the anchors and halved AC
    % and DC, so values are not comparable with their risk cut-offs.
    fprintf('Computing phase rectified signal averaging (PRSA)... \n')
    if isfield(params,'prsa_thresh') && ~isempty(params.prsa_thresh)
        thresh = params.prsa_thresh;
    else
        thresh = 20;
    end
    params.prsa_thresh = thresh;  % for exportation for users
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

%% Subfunctions
function [pwr, freqs, freq_res] = nn_spectrum(NN, NN_times, winLen, fRange, hrv_spec)
% One-sided PSD of the NN intervals of one window (s^2/Hz, or unitless for
% the normalized Lomb-Scargle), its frequencies (Hz) and their step (Hz)

% Lomb-Scargle periodogram (no resampling required). Frequency vector: at
% least 4 points per resolution step (1/window length), and at most 1/512 Hz
% apart, so that the band integral and the plotted spectrum follow the peaks
% (a step of 1/(number of beats) gave ~2 points per HF peak)
if strcmp(hrv_spec, 'LombScargle_norm') || strcmp(hrv_spec, 'LombScargle')
    df = min(1/512, 1/(4*winLen));
    fvec = fRange(1):df:fRange(2);
    if strcmp(hrv_spec, 'LombScargle_norm')
        [pwr,freqs] = plomb(NN,NN_times,fvec,'normalized');
        fprintf('Computing normalized Lomb-Scargle periodogram on the NN series... \n')
    else
        [pwr,freqs] = plomb(NN,NN_times,fvec,'psd');
        fprintf('Computing standard Lomb-Scargle periodogram on the NN series... \n')
    end

% Welch or FFT (require resampling)
else
    % Resample at 7 Hz (cubic spline) and remove the mean, whose
    % leakage would otherwise inflate the band powers
    resamp_freq = 7;
    NN_resamp = resample_NN(NN_times,NN,resamp_freq,'cub');
    NN_resamp = NN_resamp(:) - mean(NN_resamp);

    % Pwelch (one Hamming-tapered segment spanning the window; the
    % averaging is across the sliding windows), one-sided PSD in s^2/Hz
    if strcmp(hrv_spec, 'welch')
        [pwr,freqs] = pwelch(NN_resamp,length(NN_resamp),[],[],resamp_freq);
        fprintf('Computing pwelch on the NN series... \n')

    % FFT (periodogram), one-sided PSD in s^2/Hz
    else
        nSamp = length(NN_resamp);
        nHalf = floor(nSamp/2) + 1;             % 0 Hz to Nyquist
        pwr = abs(fft(NN_resamp)).^2 / (resamp_freq*nSamp);
        pwr = pwr(1:nHalf);
        pwr(2:end-1+mod(nSamp,2)) = 2*pwr(2:end-1+mod(nSamp,2));  % fold negative frequencies (not 0 Hz or Nyquist)
        freqs = (0:nHalf-1)' * resamp_freq / nSamp;
        fprintf('Computing FFT on the NN series... \n')
    end
end
freq_res = freqs(2) - freqs(1);

function cap = prsa_capacity(NN, anchors)
% Acceleration/deceleration capacity (ms) from phase-rectified signal
% averaging; anchors without two beats before and one after are dropped
anchors = anchors(anchors >= 3 & anchors <= length(NN)-1);
if isempty(anchors)
    cap = NaN; return
end
X = mean(NN(anchors + (-2:1)), 1);    % [X(-2) X(-1) X(0) X(1)]
cap = round(1000 * (X(3) + X(4) - X(2) - X(1)) / 4, 2);
