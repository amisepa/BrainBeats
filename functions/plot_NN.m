% PLOT_NN - Plot the cardiac signal with its peaks, and the RR and NN intervals.
%
% Figure panels:
%   1. Whole signal with the detected peaks (red) and the peaks kept after
%      cleaning (dark red).
%   2. RR intervals before cleaning (red dashed) and NN intervals (blue).
%   3. Scrollable 8-s window of the signal with both sets of peaks (only if
%      bb_scrollplot works on this system; otherwise panels 1-2 only).
%
% Usage:
%   plot_NN(sig_t, sig, RR_t, RR, Rpeaks, NN_t, NN, Npeaks, sigtype)
%   plot_NN(..., sigtype, unit)
%
% Inputs:
%   sig_t   - signal time vector (s), length N
%   sig     - cardiac signal, e.g. as returned by get_RR (filtered, in the
%             units of the recording)
%   RR_t    - time of each detected peak / RR interval (s), length M
%   RR      - RR intervals before cleaning (s), length M
%   Rpeaks  - sample indices of the detected peaks, length M
%   NN_t    - time of each NN interval (s), length K
%   NN      - NN intervals after cleaning (s), length K
%   Npeaks  - sample indices of the NN peaks, length K (NaN for synthetic
%             beats, whose amplitude is then interpolated at NN_t); if
%             empty or not of length K, all amplitudes are interpolated
%   sigtype - 'ecg' or 'ppg' (titles)
%   unit    - y-axis unit of sig (default 'a.u.': EEGLAB stores no channel
%             units, and recordings hold uV, mV or raw amplifier units)
%
% Copyright (C) - Cedric Cannard, 2023

function plot_NN(sig_t, sig, RR_t, RR, Rpeaks, NN_t, NN, Npeaks, sigtype, unit)

if nargin < 10 || isempty(unit), unit = 'a.u.'; end

figure("color","w",'name','Time series & corresponding RR/NN intervals', ...
    'numberTitle','off', 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
try icadefs; set(gcf, 'color', BACKCOLOR); catch; end

% Npeaks must pair with NN_t (Rpeaks does not once beats were removed or
% inserted): otherwise interpolate the signal at every NN_t
if numel(Npeaks) ~= numel(NN_t)
    if ~isempty(Npeaks)
        warning('plot_NN: Npeaks (%d) and NN_t (%d) differ in length; NN peak amplitudes are interpolated at NN_t.', numel(Npeaks), numel(NN_t))
    end
    Npeaks = nan(size(NN_t));
end

% Amplitude at NN peak positions:
%   - real beats:      exact sample sig(Npeaks), no timestamp rounding error
%   - synthetic beats: Npeaks is NaN, so interpolate the signal at NN_t
nn_amp        = nan(size(NN_t));
valid         = ~isnan(Npeaks);
nn_amp(valid) = sig(Npeaks(valid));
nn_amp(~valid)= interp1(sig_t, sig, NN_t(~valid), 'linear', NaN);

%% Subplot 3: scrollable signal
try
    subplot(3,1,3)
    win_len = 8;  % visible window (s)
    bb_scrollplot( ...
        {sig_t, sig,         'color', '#0072BD'}, {'X'}, win_len, ...
        {RR_t,  sig(Rpeaks), '+', 'MarkerSize', 10, 'LineWidth', 2, 'color', 'r',}, ...
        {NN_t,  nn_amp,      '+', 'MarkerSize', 11, 'LineWidth', 2, 'color', 'k'});
    scroll = true;
    if strcmp(sigtype, 'ecg')
        title('ECG signal + R peaks (press -> arrow to scroll)');
    else
        title('PPG signal + Pulse wave peaks (press -> arrow to scroll)');
    end
    ylabel(unit);
    xlabel('Time (s)')
    legend('signal', 'before correction', 'after correction');
catch
    warning('Scroll plot failed. Please submit an issue at: https://github.com/amisepa/BrainBeats/issues')
    scroll = false;
end

%% Subplot 1: full signal (y-limits at +/-10 SD to hide large artifacts)
if scroll, subplot(3,1,1); else, subplot(2,1,1); end
plot(sig_t, sig,         'color', '#0072BD'); hold on;
plot(RR_t,  sig(Rpeaks), '.', 'MarkerSize', 7, 'color', 'r');
plot(NN_t,  nn_amp,      '.', 'MarkerSize', 7, 'color', [0.6350 0.0780 0.1840]);
axis tight
ylim([-std(sig,'omitnan')*10  std(sig,'omitnan')*10])
if strcmp(sigtype, 'ecg')
    title('Entire ECG time series + R-peaks');
else
    title('Entire PPG time series + peaks');
end
ylabel(unit);

%% Subplot 2: RR and NN interval time series
if scroll, subplot(3,1,2); else, subplot(2,1,2); end
plot(RR_t, RR, '--', 'color', 'r',        'linewidth', 0.5); hold on;
plot(NN_t, NN, '-',  'color', '#0072BD',  'LineWidth',  1);
title('NN intervals (blue) & RR artifacts before correction (red)');
ylabel('NN intervals (s)'); xlabel('Time (s)');
legend('RR (before)', 'NN (after)');
axis tight; box on

set(findall(gcf,'type','axes'), 'fontSize', 11, 'fontweight', 'bold');
finish_figure(gcf)
