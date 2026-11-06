% FIND_BADTRIALS - Detect bad epochs from their amplitude and high-frequency content.
%
% Usage:
%   badTrials = find_badTrials(EEG, method, vis)
%   badTrials = find_badTrials(EEG, method, vis, opts)
%
% For each epoch, two metrics are computed over the analysis window
% (opts.win): amplitude (RMS across channels of each channel's RMS) and
% high-frequency noise (90th percentile across channels of the RMS of the
% data minus their 45-Hz low-pass: muscle is usually on a few channels).
% Epochs that are outliers on either metric (isoutlier) are flagged.
% The high-frequency metric needs data above 45 Hz: for HEP it is computed
% on opts.hf, the epochs of the broadband copy the ICA is fitted on (the
% HEP data are low-passed at 30 Hz by default); without it, on the data if
% opts.lowpass is above 50 Hz, else it is skipped. On 6 recordings with 15
% added muscle-like bursts each (tests/compare_badepoch_metrics.m: 300 ms
% of 30-100 Hz noise, 15 uV RMS, 10 temporal channels), this flagged 15/15
% on 5 recordings (0/15 on one with much muscle activity), with 0-3 other
% epochs on the clean ones; the former metrics (whole padded epoch, median
% absolute deviation over time, on the 0.5-30 Hz data) flagged 0-2/15.
%
% Inputs:
%   EEG       - epoched EEGLAB EEG structure
%   method    - isoutlier method: 'median' (more aggressive), 'quartiles',
%               'grubbs' (moderate, clean_eeg default), 'mean' (more lax)
%   vis       - true to plot the flagged epochs (eegplot)
%   opts      - optional struct: win ([start end] ms, default: the whole
%               epoch), hf (channels x times x epochs, the same epochs and
%               latencies as EEG.data, for the high-frequency metric),
%               lowpass (Hz, of EEG.data; default Inf)
% Outputs:
%   badTrials - indices of bad epochs (row vector)
%
% Copyright (C) - Cedric Cannard, 2022

function badTrials = find_badTrials(EEG,method,vis,opts)

if nargin < 4, opts = struct; end
disp('Detecting bad trials...')
inWin = true(1, EEG.pnts);
if isfield(opts,'win') && ~isempty(opts.win)
    inWin = EEG.times >= opts.win(1) & EEG.times < opts.win(2);
end
hfData = [];
if isfield(opts,'hf') && ~isempty(opts.hf)
    hfData = opts.hf;
elseif ~isfield(opts,'lowpass') || isempty(opts.lowpass) || opts.lowpass > 50
    hfData = EEG.data;
else
    fprintf('High-frequency metric skipped: the data are low-passed at %g Hz. \n', opts.lowpass)
end
% Low-pass FIR (order 100, pass 0-45 Hz, stop from 50 Hz)
b = design_fir(100,[2*[0 45 50]/EEG.srate 1],[1 1 0 0]);
nE = size(EEG.data,3);
sig_amp = nan(1,nE);
sig_hf = nan(1,nE);
for iEpoch = 1:nE
    sig_amp(iEpoch) = rms(rms(double(EEG.data(:,inWin,iEpoch)),2));
    if ~isempty(hfData)
        x = double(hfData(:,:,iEpoch));
        r = x - filtfilt_fast(b,1,x')';          % whole epoch filtered, then the window
        sig_hf(iEpoch) = prctile(rms(r(:,inWin),2), 90);
    end
end
badRMS = isoutlier(sig_amp,method);
badHF = false(1,nE);
if ~isempty(hfData), badHF = isoutlier(sig_hf,method); end
badTrials = unique([find(badRMS) find(badHF)]);
fprintf('Bad epochs: %d by amplitude, %d by high-frequency noise. \n', sum(badRMS), sum(badHF));

if ~isempty(badTrials)
    if vis
        % bad epochs only, with their own events; eegplot sets the channel
        % spacing from the data (bad epochs are larger than usual)
        BAD = pop_select(EEG, 'trial', badTrials);
        eegplot(BAD.data,'srate',BAD.srate,'events',BAD.event,'eloc_file',BAD.chanlocs, ...
            'winlength',min(5,numel(badTrials)),'title','Epochs removed','plottitle','Bad epochs');
        set(gcf,'Menu', 'none','Name','Bad epochs removed','NumberTitle','Off')
        finish_figure(gcf)
    end
else
    disp("No bad trials detected")
end

fprintf('Trials detected: %g \n', length(badTrials));
