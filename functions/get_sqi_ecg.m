% Caculate ECG's signal quality index(SQI) on each time window, by
% comparing two QRS time series. SQI< 0.9 should not be included in analyses.
%
% INPUTS:
%  qrs1: first annotation file
%  qrs2: second annotation file
%
% OUTPUTS:
%   sqi: sqi of each window
%   t_sqi: time index
%
% Reference when using the SQI feature:
%   Vest, Da Poian,Li, Liu, Nemati, Shah & Clifford (2018). An open-source
%   benchmarked toolbox for cardiovascular waveform and interval analysis.
%   Physiological measurement, 39(10), 105004.

% Cedric Cannard, 2022

function [sqi, t_sqi] = get_sqi_ecg(qrs1,raw_ecg,fs)

disp("Calculating signal quality index (SQI) for a ECG time series...")

% % Resample to 125 hz for method 2
% if fs ~= 125
%     [N,D] = rat(125/fs);                                       % Rational Fraction Approximation
%     % Check = [125/fs, N/D]                                      % Approximation Accuracy Check (Delete)
%     signal = resample(raw_ecg, N, D);                                   % Resampled Signal
%     % figure; plot(raw_ecg); hold on; plot(signal); legend('raw', 'resampled');
% end
% gain = 2000; % multiply by 2000 if ECG is in physical units
% qrs2 = get_rwave2(signal*gain, 125, 50);

qrs2 = get_rwave3(raw_ecg,fs);

windowlength = 10;   % default = 10
increment = 1;      % default = 1
threshold = 0.1;  % default = 0.1;
margin = 2;

qrs1 = qrs1(:)./fs;
qrs2 = qrs2(:)./fs;

endtime = max([qrs1(end), qrs2(end)]);
t = (1/fs):(1/fs):endtime;

% Create Windows
nx = floor(t(end));                 % length of sequence
overlap = windowlength-increment;   % number of overlapping elements
Nwinds = fix((nx-overlap)/(windowlength-overlap));      % number of sliding windows
t_sqi = (0:Nwinds-1) * (windowlength-overlap);          % starting index of each windows

% Calculate SQI for each Window
sqi = nan(1,length(t_sqi));
for iSeg = 1:length(t_sqi)

    % Check window for sufficient data
    if ~isnan(t_sqi(iSeg))
        % Isolate data in this window
        idx = qrs1 >= t_sqi(iSeg) & qrs1 < t_sqi(iSeg) + windowlength;

        % Normalize timing of annotation data to the windows
        a1 = qrs1(idx) - t_sqi(iSeg);
        a2 = qrs2 - t_sqi(iSeg);

        tmpSqi = run_sqi(a1,a2,threshold,margin,windowlength,fs);
        if ~isempty(tmpSqi)   % leave NaN if empty
            sqi(iSeg) = tmpSqi;
        end
    end
end

% fprintf('Mean signal quality index (SQI): %g \n', round(mean(sqi,'omitnan'),2))

%% Subfunction

function [F1,Se,PPV,Nb] = run_sqi(refqrs,testqrs,thres,margin,windowlen,fs)
% sqi = run_sqi(refqrs,testqrs,thres,margin,windowlen,fs)
% compare two sets of annotation with one as the reference (refqrs) and one
% as the test (testqrs)
%
% inputs
%     refqrs:       reference qrs annotation (in sec)
%     testqrs:      test qrs annotations (in sec)
%     thres:        threshold (in sec,default 0.05s)
%     margin:       margin time not include in comparison (in sec,default 2s)
%     windowlen:    length of the comparison window (in sec,default 60s)
%     fs:           sampling frequency
%
% output
%     sqi: match proportion according to some criteria you can change
%     depending on what you are looking for (can be Se, PPV or F1 measure).
%     See at the end of the function.
%
% When using this work, then please cite [1] and [2]:
%     [1] Behar Joachim, Oster Julien, Qiao Li, Clifford Gari D. Signal Quality
%     During Arrhythmia and its Application to False Alarm Reduction.
%     IEEE Transactions on Biomedical Engineering. 60(6). 1660-6. 2013.
%
%     [2] Li, Qiao, Roger G. Mark, and Gari D. Clifford. "Robust heart rate estimation
%     from multiple asynchronous noisy sources using signal quality indices and
%     a Kalman filter." Physiological measurement 29.1 (2008): 15.
%
% PCinCC2014, version 1.0, June 2014
% Released under the GNU General Public License
%
% Copyright (C) 2014  Joachim Behar
% Oxford university, Intelligent Patient Monitoring Group - Oxford 2014
% joachim.behar@gmail.com
%
% Updates:
% 07-02-2014
% JB- testes with Octave -> running OK
%
% 02-10-2014
% Bug fix: Dealing with annotations close to the border of the search
% window. Lines 70 - 100
% David Springer
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.

% == managing inputs
if nargin<2; error('bsqi: wrong number of input arguments \n'); end
if nargin<3; thres=0.05; end
if nargin<4; margin=2; end
if nargin<5; windowlen=60; end
if nargin<6; fs=1000; end

if size(refqrs,1)>size(refqrs,2); refqrs=refqrs';end
if size(testqrs,1)>size(testqrs,2); testqrs=testqrs';end

start = margin*fs;
stop = (windowlen-margin)*fs;
refqrs = refqrs.*fs; % convert into samples from time
testqrs = testqrs.*fs; % convert into samples from time

try
    refqrs  = refqrs(refqrs>start & refqrs<stop)'; % reference annotations
    testqrs = testqrs(testqrs>start & testqrs<stop)'; % test annotations

    if ~isempty(refqrs)

        NB_REF  = length(refqrs);
        NB_TEST = length(testqrs);

        % == removing borders problems
        indbord = find(refqrs<thres*fs | refqrs>(windowlen-thres)*fs); % reference QRS at the border
        if ~isempty(indbord)
            [IndMatchBord,DistQRSbord] = dsearchn(testqrs,refqrs(indbord));
            Indeces_below_threshold = DistQRSbord<thres*fs; %Added line to find indices of annotation of interest (refqrs) below threshold (02-10-14)
            IndMatchBord = IndMatchBord(Indeces_below_threshold);
            NB_QRS_BORD = length(indbord); % QRS at the border of the window (0,1,2)
            NB_MATCHING = length(IndMatchBord); % nb of corresponding mathing QRS (0,1,2)
            if isempty(IndMatchBord)
                refqrs(indbord) = [];
            elseif NB_MATCHING<NB_QRS_BORD
                refqrs(indbord(~Indeces_below_threshold)) = []; %Removing other indices not below threshold (02-10-14)
            end
        end

        indbord = find(testqrs<thres*fs | testqrs>(windowlen-thres)*fs); % reference QRS at the border
        if ~isempty(indbord)
            [IndMatchBord,DistQRSbord] = dsearchn(refqrs,testqrs(indbord));
            Indeces_below_threshold = DistQRSbord<thres*fs; %Added line to find indices of annotation of interest (testqrs) below threshold (02-10-14)
            IndMatchBord = IndMatchBord(Indeces_below_threshold);
            NB_QRS_BORD = length(indbord); % QRS at the border of the window (0,1,2)
            NB_MATCHING = length(IndMatchBord); % nb of corresponding mathing QRS (0,1,2)
            if isempty(IndMatchBord)
                testqrs(indbord) = [];
            elseif NB_MATCHING<NB_QRS_BORD
                testqrs(indbord(~Indeces_below_threshold)) = []; %Removing other indices not below threshold (02-10-14)
            end
        end

        % == core function
        [IndMatch,Dist] = dsearchn(refqrs,testqrs);         % closest ref for each point in test qrs
        IndMatchInWindow = IndMatch(Dist<thres*fs);         % keep only the ones within a certain window
        NB_MATCH_UNIQUE = length(unique(IndMatchInWindow)); % how many unique matching
        TP = NB_MATCH_UNIQUE;                               % number of identified ref QRS
        FN = NB_REF-TP;                                     % number of missed ref QRS
        FP = NB_TEST-TP;                                    % how many extra detection?
        Se  = TP/(TP+FN);
        PPV = TP/(FP+TP);
        F1 = 2*Se*PPV/(Se+PPV);                             % accuracy measure

        Nb.TP = TP;
        Nb.FN = FN;
        Nb.FP = FP;
    else
        F1=[];Se=[];PPV=[];Nb=[];
    end
catch
    F1=[];Se=[];PPV=[];Nb=[];
end
