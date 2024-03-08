%% Plot the filtered/scaled ECG signal with marked R peaks, NN intervals
% with RR artifacts that have been interpolated (if any). On Windows,
% users can scroll through 30-s windows of the ECG signal to inspect R
% peaks.
%
% Copyright (C) - Cedric Cannard, 2023

function plot_NN(sig_t,sig,RR_t,RR,Rpeaks,NN_t,NN,flagged,sigtype)


figure("color","w","toolbar","none","menubar","none",...
    'name','RR intervals and artifacts','numberTitle','off');

% eeglab background color
try icadefs; set(gcf, 'color', BACKCOLOR); catch; end  

try
    subplot(3,1,3)
    winSize = 30/sig_t(end);
    scrollplot({sig_t,sig,'color','#0072BD'},{RR_t,sig(Rpeaks),'.',...
        'MarkerSize',10,'color','#D95319'}, {'X'},{'Time (s)'},winSize);
    scroll = true;
    if strcmp(sigtype,'ecg')
        title('ECG signal + R peaks (press -> arrow to scroll)'); 
    else
        title('PPG signal + R peaks (press -> arrow to scroll)'); 
    end
    ylabel('μV'); xlabel('Time (s)')
catch
    warning('Scroll plot failed. Please submit an issue at: https://github.com/amisepa/BrainBeats/issues')
    scroll = false;
end

if scroll, subplot(3,1,1); else, subplot(2,1,1); end
plot(sig_t, sig,'color','#0072BD'); hold on;
plot(RR_t, sig(Rpeaks),'.','MarkerSize',10,'color','#D95319');
axis tight
if strcmp(sigtype,'ecg')
    title('ECG signal + R peaks'); 
else
    title('PPG signal + R peaks'); 
end
ylabel('μV');

if scroll, subplot(3,1,2); else, subplot(2,1,2); end
if sum(flagged) == 0
    plot(RR_t,RR,'-','color','#0072BD','linewidth',1);
else
    plot(RR_t,RR,'-','color','#A2142F','linewidth',1);
    hold on; plot(NN_t, NN,'-','color',"#0072BD", 'LineWidth', 1);
    % legend('RR artifacts','After correction')
end
title('NN intervals (blue) & RR artifacts before interpolation (red)'); ylabel('NN intervals (s)'); xlabel('Time (s)');
axis tight;
box on

% if scroll, subplot(4,1,4); else, subplot(3,1,3); end
% area(sqi(1,:), sqi(2,:)); hold on; plot(xlim,[.9 .9],'r--','linewidth',2)
% title(sprintf('Signal quality index (%g%% of RR series are artifacts)', round(SQI,1)) ); legend('','minimum threshold'); axis tight

set(findall(gcf,'type','axes'),'fontSize',11,'fontweight','bold');
