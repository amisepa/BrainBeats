clear; close all; clc

cd 'C:\Users\Tracy\Documents\MATLAB\PhysioNet-Cardiovascular-Signal-Toolbox'
load('C:\Users\Tracy\Desktop\cardio.mat')

HRVparams = InitializeHRVparams('test');
HRVparams.Fs = CARDIO.srate;

[rr,t,sqi] = Analyze_ABP_PPG_Waveforms(CARDIO.data',{'PPG'},HRVparams,[],'test');
figure; plot(t,rr,'-','color','#A2142F','linewidth',1);

[NN, tNN, tWin, AFWindows,out] = PreparDataForHRVAnlysis(rr,t,[],sqi,HRVparams,'test');
hold on; plot(tNN, NN,'-','color',"#0072BD", 'LineWidth', 1);
legend('RR artifacts','NN intervals'); title('Physionet');
ylabel('RR intervals (s)'); xlabel('Time (s)');
axis tight; box on

