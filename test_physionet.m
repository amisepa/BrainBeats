clear; close all; clc

HRVparams = InitializeHRVparams('test');
HRVparams.windowlength = floor(EEG.xmax)-1;
HRVparams.Fs = EEG.srate;
idx = contains({EEG.chanlocs.labels}, 'EXG5');
% signal = EEG.data(idx,:);
signal = ecg(1,:);

[t, rr, jqrs_ann, SQIvalue , tSQI] = ConvertRawDataToRRIntervals(signal, HRVparams, '0');
[NN, tNN, tWin, AFWindows,out] = PreparDataForHRVAnlysis(rr,t,[],[],HRVparams,subID);

figure; plot(t,rr); hold on; plot(tNN, NN); 

HRVout = [tWin' (tWin+HRVparams.windowlength)'];
HRVtitle = {'t_start' 't_end'};

% HRV time
TimeMetrics = EvalTimeDomainHRVstats(NN,tNN,sqi,HRVparams,tWin);
HRVout = [HRVout cell2mat(struct2cell(TimeMetrics))'];
