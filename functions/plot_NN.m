%% Plot the filtered/scaled ECG signal with marked R peaks, NN intervals
% with RR artifacts that have been interpolated (if any). On Windows,
% users can scroll through 30-s windows of the ECG signal to inspect R
% peaks.
%
% Copyright (C) - Cedric Cannard, 2023

function plot_NN(sig_t, sig, RR_t, RR, Rpeaks, NN_t, NN, Npeaks, sigtype)


figure("color","w","toolbar","none","menubar","none",...
    'name','RR intervals and artifacts','numberTitle','off');

% eeglab background color
try icadefs; set(gcf, 'color', BACKCOLOR); catch; end  

% Smooth RR & NN time series
try
    RR = smooth(RR, 7, 'sgolay');
    NN = smooth(NN, 7, 'sgolay');
catch
    disp("Smoothing failed (for visulazation of RR and NN intervals).")
end


try
    subplot(3,1,3)
    scrollplot({sig_t,sig,'color','#0072BD'},{'X'},15, ...
        {RR_t, sig(Rpeaks),'.','MarkerSize',12,'color',[0.6350 0.0780 0.1840]},...
        {NN_t, sig(Npeaks),'.','MarkerSize',12,'color',[0.4660 0.6740 0.1880]});
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
plot(RR_t, sig(Rpeaks),'.','MarkerSize',7,'color',[0.6350 0.0780 0.1840]);
plot(NN_t, sig(Npeaks),'.','MarkerSize',7,'color',[0.4660 0.6740 0.1880]);
axis tight
if strcmp(sigtype,'ecg')
    title('ECG signal + R peaks'); 
else
    title('PPG signal + R peaks'); 
end
ylabel('μV');

if scroll, subplot(3,1,2); else, subplot(2,1,2); end
% if sum(flagged) == 0
%     plot(RR_t,RR,'-','color','#0072BD','linewidth',0.5);
% else
plot(RR_t,RR,'-','color','#A2142F','linewidth',0.5);
hold on; plot(NN_t, NN,'-','color',"#0072BD", 'LineWidth', 1);
% legend('RR artifacts','After correction')
% end
title('NN intervals (blue) & RR artifacts before interpolation (red)'); ylabel('NN intervals (s)'); xlabel('Time (s)');
axis tight;
box on

% if scroll, subplot(4,1,4); else, subplot(3,1,3); end
% area(sqi(1,:), sqi(2,:)); hold on; plot(xlim,[.9 .9],'r--','linewidth',2)
% title(sprintf('Signal quality index (%g%% of RR series are artifacts)', round(SQI,1)) ); legend('','minimum threshold'); axis tight

set(findall(gcf,'type','axes'),'fontSize',11,'fontweight','bold');
