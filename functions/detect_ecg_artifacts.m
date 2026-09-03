function [mask, badSegments] = detect_ecg_artifacts(ECG, varargin)
% detect_ecg_artifacts - Flag bad portions of a continuous ECG recording.
%
% Marks stretches of ECG that should not be trusted for R-peak detection:
% high-frequency bursts (muscle activity, electrode rubbing, cable movement),
% optionally large slow excursions (body movement, electrode pull), and any
% non-finite samples. It is meant to run BEFORE get_RR / get_rwave3, so that
% peak detection and RR cleaning never see the corrupted stretches at all.
%
% Strategy:
%   1. Replace non-finite samples with 0 so the filters stay well posed, and
%      remember where they were.
%   2. High-pass the raw signal above HFcutoff (default 30 Hz). At that cutoff
%      the QRS complex contributes little, so what remains is dominated by
%      noise rather than by the beats themselves.
%   3. Track the moving-average power of that high-passed signal over a WinSec
%      window and convert it to a robust z score (median / 1.4826*MAD, so the
%      artifacts being measured do not inflate their own baseline). Samples
%      above HFthresh are candidate artifacts.
%   4. Smooth the candidate mask with a moving-average majority filter over
%      ErodeWin*WinSec seconds, which drops isolated flagged samples and closes
%      pinholes inside genuine bursts.
%   5. Optionally add a slow-drift mask: low-pass below 5 Hz, robust z score of
%      the amplitude itself, flag |z| > AmpThresh. Off by default (AmpThresh Inf).
%   6. Drop flagged runs shorter than MinArtSec, then merge runs separated by
%      less than MinGapSec, so the output is a small number of contiguous
%      segments rather than a speckled mask.
%
% Usage:
%   mask = detect_ecg_artifacts(ECG);
%   mask = detect_ecg_artifacts(ECG, 'HFthresh', 15, 'AmpThresh', 6);
%   [mask, seg] = detect_ecg_artifacts(ecg_vector, 'SampleRate', 500);
%
% Inputs:
%   ECG - EEGLAB structure holding one ECG channel (uses .data and .srate; if
%         .data has several channels only the first is used), OR a plain
%         numeric vector, in which case 'SampleRate' is required.
%
% Optional name-value inputs:
%   'SampleRate' (default: [])    - sampling rate (Hz). Required, and only
%                                   used, when ECG is a numeric vector.
%   'HFcutoff'   (default: 30)    - high-pass cutoff (Hz) defining the band
%                                   scanned for noise bursts.
%   'HFthresh'   (default: 5)     - robust z threshold on high-frequency power.
%                                   Lower flags more. 15 is a typical value when
%                                   only gross artifacts should be removed.
%   'AmpThresh'  (default: Inf)   - robust z threshold on the <5 Hz amplitude,
%                                   for slow movement artifacts. Inf disables
%                                   this mask entirely.
%   'WinSec'     (default: 0.5)   - moving-average window (s) for the
%                                   high-frequency power estimate. Sets the
%                                   time resolution of the detector: short
%                                   enough to isolate a burst, long enough that
%                                   a single QRS does not register as one.
%   'ErodeWin'   (default: 0.4)   - majority-filter length as a fraction of
%                                   WinSec.
%   'MinArtSec'  (default: 0.5)   - flagged runs shorter than this are dropped.
%   'MinGapSec'  (default: 0.1)   - flagged runs separated by less than this
%                                   are merged into one.
%   'Plot'       (default: false) - draw a diagnostic figure: the signal with
%                                   flagged stretches in red, and the
%                                   high-frequency z trace against its threshold.
%
% Outputs:
%   mask        - 1 x nSamples logical, true where the ECG is flagged as bad.
%   badSegments - nSeg x 2 sample indices [start stop] of the flagged runs.
%
% Notes:
%   - The mask is sample-resolved and time-aligned with the input, so it can be
%     handed straight to pop_select(..., 'nopoint', badSegments), or used to
%     drop R-peaks that fall inside a flagged stretch.
%   - Non-finite samples are folded into the mask before the MinArtSec filter,
%     so a run of NaNs shorter than MinArtSec is reported in the console but
%     does not survive into the returned mask. Handle isolated NaNs separately.
%   - AmpThresh reuses the FIR order designed for the high-pass, which is short
%     for a 5 Hz cutoff. It is a coarse movement detector, not a precise filter.
%   - Requires the Signal Processing Toolbox (fir1, filtfilt) and the Statistics
%     Toolbox (mad).
%
% See also: get_RR, clean_rr, get_sqi_ecg
%
% Changelog:
%   v1.1 - September 2026 (Cedric Cannard)
%     - Moved into BrainBeats from an external preprocessing script.
%     - WinSec is now honoured as written. It was previously passed through
%       min(WinSec, 0.5), so any value above 0.5 s silently became 0.5 s. The
%       default is 0.5 s, which reproduces the old effective behaviour; callers
%       that used to pass a larger value were already getting 0.5 s.
%     - Accepts a plain numeric vector plus 'SampleRate', not only an EEGLAB
%       structure.
%     - Returns the segment boundaries as a second output.
%     - Guards a degenerate MAD of zero (flat or constant signal) instead of
%       dividing by it and flagging on the resulting Inf/NaN z scores.
%     - Validates that HFcutoff leaves a usable filter order below Nyquist.
%     - Uses only the first channel when handed multi-channel data, instead of
%       flattening every channel into one series of the wrong length.
%
% Cedric Cannard, 2026

ip = inputParser;
addParameter(ip, 'SampleRate', []);
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

order  = 3 * round(sr / hfCut);
order  = order + mod(order, 2);
if order < 4
    error('detect_ecg_artifacts:filterTooShort', ...
        ['HFcutoff (%g Hz) is too close to the sampling rate (%g Hz) to design ' ...
         'a usable high-pass. Lower HFcutoff or resample.'], hfCut, sr);
end
b_hf   = fir1(order, hfCut/nyq, 'high');
hf_sig = filtfilt(b_hf, 1, sig);

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

fprintf('  Artifact detection:\n');
fprintf('    HF burst (>%g Hz, raw):  %.2f%%\n', hfCut, 100*mean(mask_hf));

erosionSamp = max(1, round(winSamp * erodeWin));
mask_hf = logical(conv(double(mask_hf), ones(1,erosionSamp)/erosionSamp, 'same') > 0.6);

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

mask = false(1, n);
for i = 1:size(badSegments, 1)
    mask(badSegments(i,1):badSegments(i,2)) = true;
end

fprintf('    Combined (after erosion + min-dur filter): %.2f%%\n', 100*mean(mask));
fprintf('    Segments after merging: %d\n', size(badSegments, 1));

if doPlot
    t = (0:n-1) / sr / 60;
    % Per-recording diagnostic, never closed here, so it needs the recording in
    % its Name. ECG often arrives via pop_select, so setname may be inherited or
    % missing entirely, hence the guard where it was read above.
    % '%%' (not '%%%%') -- the figure Name is a plain string, so sprintf runs
    % once and '%%' is what yields a single literal per-cent sign.
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
    plot(t, hf_z, 'Color',[0.2 0.6 0.2], 'LineWidth',0.5); hold on;
    yline(hfThr, 'r--', sprintf('thresh=%g', hfThr));
    ylabel('HF power (z)');
    title(sprintf('High-frequency power (>%g Hz)', hfCut));
    set(gca,'FontSize',11,'TickDir','out','Box','off');
    xlabel('Time (min)');
end

end

%% Helper
function segs = merge_segments(segs, minGap)
% Merge segments separated by minGap samples or fewer.
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
