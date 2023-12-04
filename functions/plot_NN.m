%% Plot the filtered/scaled ECG signal with marked R peaks, NN intervals
% with RR artifacts that have been interpolated (if any). On Windows,
% users can scroll through 30-s windows of the ECG signal to inspect R
% peaks.
%
% Copyright (C) - Cedric Cannard, 2023

function plot_NN(sig_t,sig,RR_t,RR,Rpeaks,NN_t,NN,flagged)

figure('color','w');

% Take data points from raw signal that match RR series (only for PPG?)
% if size(sig,2) ~= RR
%     idx = ismember(sig_t,RR_t);
%     sig_t = sig_t(idx);
%     sig = sig(idx);
% end

try
    subplot(3,1,3)
    winSize = 30/sig_t(end);
    scrollplot({sig_t,sig,'color','#0072BD'},{RR_t,sig(Rpeaks),'.','MarkerSize',10,'color','#D95319'}, {'X'},{'Time (s)'},winSize);
    scroll = true;
catch
    warning('Scroll plot failed. Please submit an issue at: https://github.com/amisepa/BrainBeats/issues')
    scroll = false;
end

if scroll, subplot(3,1,1); else, subplot(2,1,1); end
plot(sig_t, sig,'color','#0072BD'); hold on;
plot(RR_t, sig(Rpeaks),'.','MarkerSize',10,'color','#D95319');
axis tight
title('Filtered ECG signal + R peaks'); ylabel('mV');

if scroll, subplot(3,1,2); else, subplot(2,1,2); end
if sum(flagged) == 0
    plot(RR_t,RR,'-','color','#0072BD','linewidth',1);
else
    plot(RR_t,RR,'-','color','#A2142F','linewidth',1);
    hold on; plot(NN_t, NN,'-','color',"#0072BD", 'LineWidth', 1);
    % legend('RR artifacts','NN intervals')
end
title('RR (red) and NN (blue) intervals'); ylabel('RR intervals (s)'); xlabel('Time (s)');
axis tight; box on

% if scroll, subplot(4,1,4); else, subplot(3,1,3); end
% area(sqi(1,:), sqi(2,:)); hold on; plot(xlim,[.9 .9],'r--','linewidth',2)
% title(sprintf('Signal quality index (%g%% of RR series are artifacts)', round(SQI,1)) ); legend('','minimum threshold'); axis tight

set(gcf,'Toolbar','none','Menu','none');  % remove toolbobar and menu
set(gcf,'Name','RR intervals and artifacts','NumberTitle','Off')  % name
set(findall(gcf,'type','axes'),'fontSize',10,'fontweight','bold');
