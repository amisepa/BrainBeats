function [mask, badSegments] = detect_ecg_artifacts(ECG, varargin)
% DETECT_ECG_ARTIFACTS - Flag bad stretches of a continuous ECG recording.
%
% Flags high-frequency bursts (muscle, electrode rubbing, cable movement),
% optionally large slow excursions (body movement), and non-finite samples,
% so they can be excluded before R-peak detection (get_RR).
%
% Usage:
%   [mask, badSegments] = detect_ecg_artifacts(ECG, 'key', val, ...)
%
% Inputs:
%   ECG - EEGLAB structure (.data, .srate; only the first channel is used)
%         or a numeric vector (then 'SampleRate' is required)
% Optional name-value pairs:
%   'SampleRate' - sampling rate (Hz), for vector input only (default [])
%   'HFcutoff'   - highpass cutoff (Hz) of the noise band (default 30)
%   'Method'     - 'moving' (default): sample by sample, as up to v1.6;
%                  'window': sliding windows of WinSec every StepSec, as
%                  clean_rawdata does (see Method below; recommended for long
%                  or sleep recordings)
%   'HFthresh'   - robust z threshold on HF power; lower flags more
%                  (default 5; ~15 to catch only gross artifacts;
%                  'window': 8, on log power)
%   'AmpThresh'  - 'moving' only: robust z threshold on the < 5 Hz amplitude
%                  (default Inf = off)
%   'AmpRatio'   - 'window' only: a window is flagged when the range of its
%                  < 5 Hz signal exceeds AmpRatio x the median range of the
%                  windows within AmpLocalSec (default 2.5; Inf = off)
%   'AmpLocalSec'- 'window' only: half-width of that local reference (s,
%                  default 60): a stretch of taller R waves raises its own
%                  reference, a movement of a few seconds does not
%   'HighPass'   - the ECG is first high-passed at this edge (Hz) as
%                  pop_eegfiltnew(EEG,'locutoff',HighPass) does (default 0 =
%                  off; 'window': 1). Baseline drift is not an artifact: on an
%                  unfiltered ECG the < 5 Hz rule flags it.
%   'WinSec'     - window of the HF power (s, default 0.5; 'window': 1):
%                  short enough to isolate a burst, long enough that one QRS
%                  does not register as one
%   'StepSec'    - 'window' only: step between windows (s, default 0.3, i.e.
%                  70% overlap with 1-s windows)
%   'ErodeWin'   - 'moving' only: majority-filter length, fraction of WinSec
%                  (default 0.4)
%   'MinArtSec'  - flagged runs shorter than this are dropped (s, default 0.5)
%   'MinGapSec'  - runs this close or closer are merged (s, default 0.1;
%                  'window': 1)
%   'Plot'       - diagnostic figure (default false)
%
% Outputs:
%   mask        - 1 x nSamples logical, true where the ECG is bad
%   badSegments - nSeg x 2 [start stop] sample indices of the flagged runs,
%                 usable with pop_select(..., 'nopoint', badSegments)
%
% Method: HF power = moving average of the squared > HFcutoff signal (where
% the QRS contributes little), as a robust z score (median, 1.4826 x MAD,
% so artifacts do not inflate their own baseline); then a 60% majority
% filter, the optional slow-amplitude mask, the MinArtSec and MinGapSec rules.
% 'window': on the ECG high-passed at HighPass Hz, per window of WinSec,
% every StepSec, the log mean > HFcutoff power as a robust z score across the
% recording's windows, and the range (max - min) of the < 5 Hz signal as a
% ratio to the median range of the windows within AmpLocalSec; a sample is
% flagged when any window holding it is flagged, then the MinArtSec and
% MinGapSec rules. The range of a window catches a slow movement excursion
% whole, where the sample-by-sample < 5 Hz amplitude flags only its extremes.
% Limits: a movement lasting about AmpLocalSec or longer raises its own local
% reference and is not flagged by the range rule (only by HFthresh, if it
% carries high-frequency power); the last StepSec or less after the last
% full window is not assessed. On an overnight sleep ECG (200 Hz),
% scored against each heartbeat's own ECG shape (correlation < 0.8 with the
% median heartbeat, or peak amplitude outside 0.5-2x; 180 distorted
% heartbeats in REM and N2 sleep, all of abnormal shape), 'window' with its
% defaults caught 137 and dropped 0.98% of the clean REM heartbeats (those
% within 2 s of a flag). Against the night's median range instead of the
% local one (a robust z of 6: 1.6 x the median) it caught 158 for 1.7%, and
% flagged stretches of clean, taller R waves; without the amplitude rule, 61
% (bursts and saturation only); without the high-pass, baseline drift;
% 'moving' with HFthresh 5 caught 168 with AmpThresh 6 (7.6% dropped) and 9
% with AmpThresh Inf. On a second night (2 distorted heartbeats), 'window'
% caught both and dropped no clean REM heartbeat, 'moving' neither.
% Thresholds from two recordings: check yours.
%
% Notes:
%   - Non-finite samples are set to 0 for filtering and added to the mask
%     before the MinArtSec rule, so NaN runs shorter than MinArtSec are
%     reported but not returned.
%   - The < 5 Hz filter reuses the highpass FIR order, which is short for
%     5 Hz: a coarse movement detector.
%   - Requires the Signal Processing (fir1, filtfilt) and Statistics (mad)
%     toolboxes.
%
% See also: get_RR, clean_rr, get_sqi_ecg
%
% Copyright (C) - Cedric Cannard, 2026

ip = inputParser;
addParameter(ip, 'SampleRate', []);
addParameter(ip, 'Method',    'moving', @(x) any(strcmpi(x, {'moving','window'})));
addParameter(ip, 'StepSec',   0.3);
addParameter(ip, 'HighPass',  0);
addParameter(ip, 'AmpRatio',  2.5);
addParameter(ip, 'AmpLocalSec', 60);
addParameter(ip, 'HFcutoff',  30);
addParameter(ip, 'HFthresh',  5);
addParameter(ip, 'AmpThresh', Inf);
addParameter(ip, 'WinSec',    0.5);
addParameter(ip, 'MinGapSec', 0.1);
addParameter(ip, 'MinArtSec', 0.5);
addParameter(ip, 'ErodeWin',  0.4);
addParameter(ip, 'Plot',      false);
parse(ip, varargin{:});

hfCut     = ip.Results.HFcutoff;
hfThr     = ip.Results.HFthresh;
ampThr    = ip.Results.AmpThresh;
winSec    = ip.Results.WinSec;
minGap    = ip.Results.MinGapSec;
minArtSec = ip.Results.MinArtSec;
erodeWin  = ip.Results.ErodeWin;
doPlot    = ip.Results.Plot;
method    = lower(char(ip.Results.Method));
stepSec   = ip.Results.StepSec;
hpEdge    = ip.Results.HighPass;
ampRatio  = ip.Results.AmpRatio;
ampLocal  = ip.Results.AmpLocalSec;
if strcmp(method, 'window')
    % the window method's own defaults, for the values not given
    dflt = ip.UsingDefaults;
    if ismember('WinSec', dflt),    winSec = 1;  end
    if ismember('HFthresh', dflt),  hfThr  = 8;  end
    if ismember('HighPass', dflt),  hpEdge = 1;  end
    if ismember('MinGapSec', dflt), minGap = 1;  end
end

% Accept an EEGLAB structure or a bare vector. A structure carries its own
% srate, so SampleRate is only consulted for the vector form.
setname = '';
if isstruct(ECG)
    if ~isfield(ECG,'data') || ~isfield(ECG,'srate')
        error('detect_ecg_artifacts:badInput', ...
            'ECG structure must have .data and .srate fields.');
    end
    sr = ECG.srate;
    if size(ECG.data,1) > 1
        warning('detect_ecg_artifacts:multiChannel', ...
            'ECG has %d channels; using the first one.', size(ECG.data,1));
    end
    sig = double(ECG.data(1,:));
    if isfield(ECG,'setname'), setname = char(ECG.setname); end
else
    sr = ip.Results.SampleRate;
    if isempty(sr)
        error('detect_ecg_artifacts:noSampleRate', ...
            'Pass ''SampleRate'' when ECG is a numeric vector.');
    end
    sig = double(ECG(:)');
end
n   = numel(sig);
nyq = sr / 2;

if hfCut >= nyq
    error('detect_ecg_artifacts:cutoffAboveNyquist', ...
        'HFcutoff (%g Hz) must be below Nyquist (%g Hz).', hfCut, nyq);
end

% Flag and replace non-finite samples before filtering
nonfinite_mask = ~isfinite(sig);
if any(nonfinite_mask)
    fprintf('    Warning: %d non-finite samples in ECG (%.2f%%) -- replacing with 0 for filtering.\n', ...
        sum(nonfinite_mask), 100*mean(nonfinite_mask));
    sig(nonfinite_mask) = 0;
end
if hpEdge > 0
    sig = highpass_eegfiltnew(sig, sr, hpEdge);
end

% High-frequency (noise) band: FIR highpass, even order, zero-phase
order  = 3 * round(sr / hfCut);
order  = order + mod(order, 2);
if order < 4
    error('detect_ecg_artifacts:filterTooShort', ...
        ['HFcutoff (%g Hz) is too close to the sampling rate (%g Hz) to design ' ...
         'a usable high-pass. Lower HFcutoff or resample.'], hfCut, sr);
end
b_hf   = fir1(order, hfCut/nyq, 'high');
hf_sig = filtfilt(b_hf, 1, sig);

fprintf('  Artifact detection (%s):\n', method);
if strcmp(method, 'window')
    w  = max(1, round(winSec * sr));
    st = max(1, round(stepSec * sr));
    if n < w
        error('detect_ecg_artifacts:tooShort', ...
            'The ECG (%d samples) is shorter than one window (%d).', n, w);
    end
    starts = 1:st:(n - w + 1);
    c    = [0 cumsum(hf_sig.^2)];
    pw   = (c(starts + w) - c(starts)) / w;           % mean HF power per window
    hf_z = robust_z(log(max(pw, 1e-12 * max(pw))));   % floor relative to the data: any units
    bad  = false(size(starts));
    if isempty(hf_z)
        warning('detect_ecg_artifacts:degenerateScale', ...
            'High-frequency power has zero spread; no HF artifact can be detected.');
        hf_z = zeros(size(starts));
    else
        bad = hf_z > hfThr;
    end
    fprintf('    HF burst windows (>%g Hz):  %.2f%%\n', hfCut, 100*mean(bad));
    if isfinite(ampRatio)
        b_lf   = fir1(order, 5/nyq, 'low');
        lf_sig = filtfilt(b_lf, 1, sig);
        rng_lf = movmax(lf_sig, [0 w-1]) - movmin(lf_sig, [0 w-1]);
        rng_lf = rng_lf(starts);
        k      = round(ampLocal * sr / st);          % windows on each side
        loc    = movmedian(rng_lf, [k k]);           % ends: shrinking window
        badA   = rng_lf > ampRatio * loc;
        fprintf('    Slow movement windows (<5Hz): %.2f%%\n', 100*mean(badA));
        bad = bad | badA;
    end
    % a sample is flagged when any flagged window covers it
    d = accumarray(starts(bad)', 1, [n+1 1])' - accumarray((starts(bad) + w)', 1, [n+1 1])';
    mask_hf = cumsum(d(1:n)) > 0;
    mask_lf = false(1, n);
    hf_t    = (starts - 1 + w/2) / sr / 60;           % window centres (min), for the plot
else
winSamp = max(1, round(winSec * sr));
hf_pow  = movmean(hf_sig.^2, winSamp);
% Robust scale. MAD is zero for a flat or constant-power signal, and dividing
% by it would turn every sample into Inf or NaN and then flag all of them.
hf_scale = 1.4826 * mad(hf_pow, 1);
if hf_scale <= 0 || ~isfinite(hf_scale)
    warning('detect_ecg_artifacts:degenerateScale', ...
        'High-frequency power has zero spread; no HF artifact can be detected.');
    hf_z    = zeros(1, n);
    mask_hf = false(1, n);
else
    hf_z    = (hf_pow - median(hf_pow)) / hf_scale;
    mask_hf = hf_z > hfThr;
end

fprintf('    HF burst (>%g Hz, raw):  %.2f%%\n', hfCut, 100*mean(mask_hf));

% Majority filter (> 60% flagged within ErodeWin*WinSec): drops isolated
% samples and closes pinholes inside bursts
erosionSamp = max(1, round(winSamp * erodeWin));
mask_hf = logical(conv(double(mask_hf), ones(1,erosionSamp)/erosionSamp, 'same') > 0.6);

% Optional slow-movement mask: robust z of the < 5 Hz amplitude
if isfinite(ampThr)
    b_lf   = fir1(order, 5/nyq, 'low');
    lf_sig = filtfilt(b_lf, 1, sig);
    lf_scale = 1.4826 * mad(lf_sig, 1);
    if lf_scale <= 0 || ~isfinite(lf_scale)
        warning('detect_ecg_artifacts:degenerateScale', ...
            'Low-frequency amplitude has zero spread; no movement artifact can be detected.');
        mask_lf = false(1, n);
    else
        lf_z    = (lf_sig - median(lf_sig)) / lf_scale;
        mask_lf = abs(lf_z) > ampThr;
    end
    fprintf('    Slow movement (<5Hz): %.2f%%\n', 100*mean(mask_lf));
else
    mask_lf = false(1, n);
end
hf_t = (0:n-1) / sr / 60;
end

% Combine, convert to [start stop] runs, drop short runs, merge close ones
mask = mask_hf | mask_lf | nonfinite_mask;
minArtSamp = round(minArtSec * sr);
badSegments      = reshape(find(diff([false mask false])), 2, [])';
badSegments(:,2) = badSegments(:,2) - 1;
if ~isempty(badSegments)
    tooShort = (badSegments(:,2) - badSegments(:,1) + 1) < minArtSamp;
    badSegments(tooShort, :) = [];
end

if ~isempty(badSegments)
    badSegments = merge_segments(badSegments, round(minGap * sr));
end

% Rebuild the mask from the final segments
mask = false(1, n);
for i = 1:size(badSegments, 1)
    mask(badSegments(i,1):badSegments(i,2)) = true;
end

fprintf('    Combined (after erosion + min-dur filter): %.2f%%\n', 100*mean(mask));
fprintf('    Segments after merging: %d\n', size(badSegments, 1));

if doPlot
    t = (0:n-1) / sr / 60;
    % The figure is left open, so its Name identifies the recording (setname
    % may be empty). sprintf runs once on it, so '%%' gives one literal '%'.
    figure('Color','w','Position',[100 100 1400 500], 'NumberTitle','off', ...
        'Name', sprintf('ECG artifacts | %s | %.1f%% flagged | %d seg', ...
                        setname, 100*mean(mask), size(badSegments,1)));
    subplot(2,1,1);
    plot(t, sig, 'Color',[0.2 0.2 0.6], 'LineWidth',0.5); hold on;
    sig_bad = sig; sig_bad(~mask) = NaN;
    plot(t, sig_bad, 'r', 'LineWidth',0.8);
    ylabel('ECG');
    title(sprintf('ECG artifact detection  (%.2f%% flagged)', 100*mean(mask)));
    legend({'Clean','Artifact'}, 'Box','off', 'Location','best');
    set(gca,'FontSize',11,'TickDir','out','Box','off');
    subplot(2,1,2);
    plot(hf_t, hf_z, 'Color',[0.2 0.6 0.2], 'LineWidth',0.5); hold on;
    yline(hfThr, 'r--', sprintf('thresh=%g', hfThr));
    ylabel('HF power (z)');
    title(sprintf('High-frequency power (>%g Hz)', hfCut));
    set(gca,'FontSize',11,'TickDir','out','Box','off');
    xlabel('Time (min)');
    finish_figure(gcf)
end

end

%% Helpers
function y = highpass_eegfiltnew(x, sr, edge)
% High-pass as pop_eegfiltnew(EEG,'locutoff',edge) without EEGLAB: firws
% windowed sinc (Hamming), transition band min(max(0.25*edge, 2), edge) Hz,
% -6 dB at edge minus half of it, order 3.3/(band/sr) rounded up to even, one
% pass made zero-phase, the signal padded with its first and last values.
df    = min(max(0.25 * edge, 2), edge);
m     = ceil(3.3 / (df / sr) / 2) * 2;
fc    = (edge - df/2) / sr;                % cycles per sample
k     = -m/2:m/2;
b     = sin(2*pi*fc*k) ./ k;
b(k == 0) = 2*pi*fc;
b     = b .* hamming(m + 1)';
b     = -b / sum(b);                       % unity gain at DC, then spectral inversion
b(m/2 + 1) = b(m/2 + 1) + 1;
h     = m / 2;
y     = conv([repmat(x(1), 1, h) x repmat(x(end), 1, h)], b, 'valid');
end

function z = robust_z(x)
% Robust z score (median, 1.4826 x MAD); empty when the spread is zero.
s = 1.4826 * mad(x, 1);
if s <= 0 || ~isfinite(s), z = []; else, z = (x - median(x)) / s; end
end

function segs = merge_segments(segs, minGap)
% Merge consecutive segments whose start-to-previous-stop distance is at
% most minGap samples.
if isempty(segs), return; end
i = 1;
while i < size(segs, 1)
    if segs(i+1,1) - segs(i,2) <= minGap
        segs(i,2)   = segs(i+1,2);
        segs(i+1,:) = [];
    else
        i = i + 1;
    end
end
end
