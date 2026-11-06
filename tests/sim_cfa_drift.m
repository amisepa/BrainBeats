function R = sim_cfa_drift(files, clean)
% SIM_CFA_DRIFT - Sensitivity of CFA_MAP's test of a change of the QRS map
% during a recording, by simulation on real recordings.
%
% Each recording is run through BrainBeats (HEP, ECG), then the cardiac
% field is rotated relative to the head for part of the heartbeats: a head
% turn (yaw, about the vertical axis) or nod (pitch, about the left-right
% axis) of 5-30 degrees, for the second half of the heartbeats ('step') or
% for 35-65% of them ('transient'). The field is taken as fixed while the
% electrodes move in it: the rotated heartbeats get (W - I) * erp, W the
% spherical-spline interpolation from the electrode positions to the
% rotated ones and erp the average heartbeat (conductivity changes of the
% neck and torso are ignored). The recording as it is gives the false
% alarms.
%
% Usage:
%   R = sim_cfa_drift                      % the two sample datasets, cleaned
%   R = sim_cfa_drift(files, clean)        % e.g. ds003838 recordings
%
% Inputs:
%   files - cell array of .set files with an 'ECG' channel (default: the
%           tutorial dataset and sub-047 of sample_data)
%   clean - run the EEG cleaning (default true)
%
% Output:
%   R - table: file, scenario, axis, angle (deg), map_r (correlation of the
%       rotated and the original QRS map), p_change, rho (at the largest
%       change), changed
%
% Results (2026-11-04; 12 recordings of 120-370 heartbeats: the two sample
% datasets and sub-062, 033-036, 038-042 of ds003838 'derived', raw and
% cleaned): no false alarm on the 24; 'step' yaw 30 deg detected in 3-4/12
% (map r .51-.82 raw, .66-.81 cleaned), 20 deg in 1-2/12 (.67-.91), 10 deg
% in 0-1/12 (.90-.97); pitch 30 deg in 0-1/12; 'transient' yaw 30 deg in
% 1/12.
%
% Runtime: ~3 min per recording with the cleaning.
%
% Copyright (C) - Cedric Cannard, 2026

bbroot = fileparts(fileparts(mfilename('fullpath')));
if nargin < 1 || isempty(files)
    files = {fullfile(bbroot,'sample_data','dataset.set'), fullfile(bbroot,'sample_data','sub-047_rest.set')};
end
if nargin < 2, clean = true; end
if ~exist('pop_loadset','file'), eeglab; close; end
p = strsplit(path, pathsep);
old = p(contains(p,'BrainBeats','IgnoreCase',true) & ~startsWith(p,bbroot));
if ~isempty(old), rmpath(old{:}); end
addpath(bbroot, fullfile(bbroot,'functions'));

angles = [5 10 20 30];
axName = {'yaw' 'pitch'};
R = table();
for k = 1:numel(files)
    [fp, fn, ext] = fileparts(files{k});
    EEG = pop_loadset('filename',[fn ext],'filepath',fp);
    EEG = pop_select(EEG,'nochannel',intersect({'PPG'},{EEG.chanlocs.labels}));
    H = brainbeats_process(EEG,'analysis','hep','heart_signal','ECG','heart_channels',{'ECG'}, ...
        'clean_eeg',clean,'icamethod',1,'hep_surrogates',0,'hep_cfa',false, ...
        'vis_cleaning',false,'vis_outputs',false,'save',false,'gong',false);
    isE = ~ismember(upper({H.chanlocs.labels}), {'ECG','PPG'});
    X = double(H.data(isE,:,:)); loc = H.chanlocs(isE);
    n = size(X,3); bt = (1:n) * 0.8;      % heartbeat times: only used for the reported minutes
    Xm = X - mean(X(:, H.times < -50, :), 2);
    erp = mean(Xm, 3);
    c0 = cfa_map(X, H.times, bt);
    R = [R; row(fn, 'none', '-', 0, NaN, c0)]; %#ok<AGROW>
    % spherical-spline interpolation matrix to the rotated positions
    M = eeg_emptyset; M.data = eye(numel(loc)); M.nbchan = numel(loc); M.pnts = numel(loc);
    M.trials = 1; M.srate = 1; M.xmin = 0; M.chanlocs = loc; M = eeg_checkset(M);
    for ax = 1:2
        for a = angles
            co = cosd(a); si = sind(a);
            if ax == 1, Rot = [co -si 0; si co 0; 0 0 1]; else, Rot = [co 0 si; 0 1 0; -si 0 co]; end
            L = loc;
            for i = 1:numel(L)
                v = Rot * [L(i).X; L(i).Y; L(i).Z];
                L(i).X = v(1); L(i).Y = v(2); L(i).Z = v(3); L(i).labels = sprintf('rot%d', i);
            end
            L = convertlocs(L, 'cart2all');
            Mi = eeg_interp(M, L, 'spherical');
            [~, j] = ismember({L.labels}, {Mi.chanlocs.labels});
            W = double(Mi.data(j,:));
            for sc = {'step' 'transient'}
                if strcmp(sc{1},'step'), sel = floor(n/2)+1:n; else, sel = round(.35*n)+1:round(.65*n); end
                Xr = X;
                Xr(:,:,sel) = Xr(:,:,sel) + (W - eye(size(W))) * erp;
                c = cfa_map(Xr, H.times, bt);
                R = [R; row(fn, sc{1}, axName{ax}, a, corr(W*c0.map, c0.map), c)]; %#ok<AGROW>
            end
        end
    end
    fprintf('%s: false alarm %d; detected (step, yaw 10/20/30 deg): %s\n', fn, R.changed(find(strcmp(R.file,fn),1)), ...
        mat2str(R.changed(strcmp(R.file,fn) & strcmp(R.scenario,'step') & strcmp(R.axis,'yaw') & R.angle >= 10)'));
end


function T = row(fn, sc, ax, a, mr, c)
T = table({fn}, {sc}, {ax}, a, mr, c.p_change, c.rho, c.changed, ...
    'VariableNames', {'file' 'scenario' 'axis' 'angle' 'map_r' 'p_change' 'rho' 'changed'});
