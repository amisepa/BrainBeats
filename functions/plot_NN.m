% plot_NN - Plot cardiac signal with R-peaks, RR intervals, and NN intervals.
%
% Produces a 3-panel figure:
%   Panel 1: Full signal time series with original R-peaks (red) and
%            corrected NN peaks (dark red) overlaid.
%   Panel 2: RR interval time series (red dashed) vs. cleaned NN intervals
%            (blue), allowing visual inspection of artifact correction.
%   Panel 3: Scrollable signal window (10 s) for detailed peak inspection.
%            Only available on systems where scrollplot is functional.
%
% NN peak amplitudes are estimated via linear interpolation of the signal
% at NN timestamps, which correctly handles synthetic beats (whose sample
% indices are NaN) without crashing.
%
% Usage:
%   plot_NN(sig_t, sig, RR_t, RR, Rpeaks, NN_t, NN, Npeaks, sigtype)
%
% Inputs:
%   sig_t   - signal timestamps (s), vector (1 x N)
%   sig     - raw cardiac signal (ECG in uV or PPG in a.u.), vector (1 x N)
%   RR_t    - timestamps of original RR intervals (s), vector (1 x M)
%   RR      - original RR interval durations (s), vector (1 x M)
%   Rpeaks  - sample indices of original detected R-peaks, vector (1 x M)
%   NN_t    - timestamps of cleaned NN intervals (s), vector (1 x K)
%   NN      - cleaned NN interval durations (s), vector (1 x K)
%   Npeaks  - sample indices of cleaned peaks (NaN for synthetic beats),
%             vector (1 x K); if empty, Rpeaks is used instead
%   sigtype - signal type string: 'ecg' or 'ppg'
%
% Copyright (C) - Cedric Cannard, 2023

function plot_NN(sig_t, sig, RR_t, RR, Rpeaks, NN_t, NN, Npeaks, sigtype)

figure("color","w",'name','Time series & corresponding RR/NN intervals', ...
    'numberTitle','off', 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
try icadefs; set(gcf, 'color', BACKCOLOR); catch; end

% Use RR peaks if no NN peaks provided
if isempty(Npeaks)
    Npeaks = Rpeaks;
end

% Amplitude at NN peak positions:
%   - real beats:      use exact sample index sig(nPeaks) - avoids timestamp rounding error
%   - synthetic beats: fall back to interp1 since nPeaks is NaN
nn_amp        = nan(size(NN_t));
valid         = ~isnan(Npeaks);
nn_amp(valid) = sig(Npeaks(valid));
nn_amp(~valid)= interp1(sig_t, sig, NN_t(~valid), 'linear', NaN);

%% Subplot 3: scrollable signal
try
    subplot(3,1,3)
    win_len = 8;
    scrollplot( ...
        {sig_t, sig,         'color', '#0072BD'}, {'X'}, win_len, ...
        {RR_t,  sig(Rpeaks), '+', 'MarkerSize', 10, 'LineWidth', 2, 'color', 'r',}, ...
        {NN_t,  nn_amp,      '+', 'MarkerSize', 11, 'LineWidth', 2, 'color', 'k'});
    scroll = true;
    if strcmp(sigtype, 'ecg')
        title('ECG signal + R peaks (press -> arrow to scroll)');
        ylabel('μV');
    else
        title('PPG signal + Pulse wave peaks (press -> arrow to scroll)');
        ylabel('a.u.');
    end
    xlabel('Time (s)')
    legend('signal', 'before correction', 'after correction');
catch
    warning('Scroll plot failed. Please submit an issue at: https://github.com/amisepa/BrainBeats/issues')
    scroll = false;
end

%% Subplot 1: full signal
if scroll, subplot(3,1,1); else, subplot(2,1,1); end
plot(sig_t, sig,         'color', '#0072BD'); hold on;
plot(RR_t,  sig(Rpeaks), '.', 'MarkerSize', 7, 'color', 'r');
plot(NN_t,  nn_amp,      '.', 'MarkerSize', 7, 'color', [0.6350 0.0780 0.1840]);
axis tight
ylim([-std(sig,'omitnan')*7  std(sig,'omitnan')*7])
if strcmp(sigtype, 'ecg')
    title('Entire ECG time series + R-peaks');
    ylabel('μV');
else
    title('Entire PPG time series + peaks');
    ylabel('a.u.')
end

%% Subplot 2: NN interval time series
if scroll, subplot(3,1,2); else, subplot(2,1,2); end
plot(RR_t, RR, '--', 'color', 'r',        'linewidth', 0.5); hold on;
plot(NN_t, NN, '-',  'color', '#0072BD',  'LineWidth',  1);
title('NN intervals (blue) & RR artifacts before correction (red)');
ylabel('NN intervals (s)'); xlabel('Time (s)');
legend('RR (before)', 'NN (after)');
axis tight; box on

set(findall(gcf,'type','axes'), 'fontSize', 11, 'fontweight', 'bold');