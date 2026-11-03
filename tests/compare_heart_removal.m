function T = compare_heart_removal(files, outfile, synthetic, methods)
% COMPARE_HEART_REMOVAL - Heart artifact left in the HEP after each removal
% method, and how much of a heartbeat-locked brain response each keeps, for
% several datasets.
%
% For each dataset (EEG + an 'ECG' channel, continuous), runs the HEP
% analysis with the EEG cleaning (Picard ICA) and each 'heart_removal' option
% (methods, default all: 'none' 'ica' 'ica_caracas' 'ica_ecg' 'ecg_regression'
% 'ssp'). Reports the heart-locked EEG
% amplitude (CFA_AMPLITUDE: RMS over channels of the heartbeat-locked
% average, -50 to 100 ms and 200 to 500 ms) before any removal and left
% after each method, and the number of HEP epochs. With synthetic = true,
% each method is run again on the data with a synthetic response added after
% each R-peak on the frontocentral channels (1 uV, Gaussian at 300 ms, SD
% 40 ms, amplitude jitter 30%), and the share of it left in the output is
% reported ('kept', %: projection of the difference of the ROI averages on
% the response, 200-400 ms).
%
% Usage:
%   T = compare_heart_removal({'sub-01.set' 'sub-02.set'})
%   T = compare_heart_removal(files, 'heart_removal.csv', true)
%
% Inputs:
%   files     - cell array of .set files (full paths)
%   outfile   - optional .csv file for the table ('' for none)
%   synthetic - also measure the synthetic response kept (default false;
%               doubles the run time)
%   methods   - cell array of 'heart_removal' options (default all; 'none'
%               is always run first, for kept_vs_none)
%
% Output:
%   T - table, one row per dataset and method: file, method, qrs_before,
%       qrs_after, late_before, late_after (uV), qrs_reduction (%), epochs,
%       kept (%, NaN without synthetic), kept_vs_none (% of the 'none' value)
%
% Notes:
% - 'before' is the cleaned EEG (filter, reference, bad channels) before
%   heart removal; 'after' is the EEG after all cleaning. The 200-500 ms
%   amplitude also holds any neural HEP, so it bounds the residual artifact
%   from above.
% - The random number generator is reset before each run so that the ICA
%   starts the same way with and without the synthetic response; the epochs
%   kept can still differ by a few.
%
% Copyright (C) - Cedric Cannard, 2026

if nargin < 2, outfile = ''; end
if nargin < 3 || isempty(synthetic), synthetic = false; end
if ~exist('pop_loadset','file'), eeglab; close; end
bb = fileparts(fileparts(mfilename('fullpath')));
p = strsplit(path, pathsep);
old = p(contains(p,'BrainBeats','IgnoreCase',true) & ~startsWith(p,bb));
if ~isempty(old), rmpath(old{:}); end
addpath(bb, fullfile(bb,'functions'));

files = cellstr(files);
common = {'heart_signal','ECG','heart_channels',{'ECG'},'clean_eeg',true,'icamethod',1, ...
    'vis_cleaning',false,'vis_outputs',false,'save',false,'gong',false};
roiLabels = {'F1' 'F2' 'F3' 'F4' 'Fz' 'FC1' 'FC2' 'FC3' 'FC4' 'FCz' 'C1' 'C2' 'C3' 'C4' 'Cz'};
if nargin < 4 || isempty(methods), methods = {'none' 'ica' 'ica_caracas' 'ica_ecg' 'ecg_regression' 'ssp'}; end
methods = [{'none'} setdiff(methods, {'none'}, 'stable')];
rows = {};
for k = 1:numel(files)
    [fp, fn, ext] = fileparts(files{k});
    EEG = pop_loadset('filename',[fn ext],'filepath',fp);
    labels = {EEG.chanlocs.labels};
    other = labels(contains(lower(labels), {'ppg' 'pleth' 'eog' 'emg' 'gsr' 'eda'}));
    if ~isempty(other), EEG = pop_select(EEG,'nochannel',other); end

    % Synthetic response after each R-peak (R-peaks from the ECG, as the
    % analysis finds them)
    if synthetic
        fs = EEG.srate;
        iECG = find(strcmpi({EEG.chanlocs.labels},'ECG'));
        [~,~,rp] = get_RR(double(EEG.data(iECG,:)), EEG.times, struct('heart_signal','ecg','fs',fs));
        roiIn = ismember(lower({EEG.chanlocs.labels}), lower(roiLabels));
        t = (round(-0.3*fs):round(0.7*fs))/fs;
        bump = exp(-(t-0.3).^2/(2*0.04^2));
        rng(1); gains = 1 + 0.3*randn(numel(rp),1);
        S = EEG;
        okB = rp + round(t(1)*fs) > 0 & rp + round(t(end)*fs) <= EEG.pnts;
        for b = find(okB(:))'
            idx = rp(b) + round(t*fs);
            S.data(roiIn,idx) = S.data(roiIn,idx) + gains(b) * bump;
        end
        gain = mean(gains(okB));
    end

    for m = methods
        R = run_method(EEG, m{1}, common);
        pre = R.brainbeats.preprocessings;
        if strcmp(m{1},'ecg_regression')   % amplitudes before the regression
            qb = pre.heart_regression.cfa_before; lb = pre.heart_regression.cfa_late_before;
        elseif strcmp(m{1},'ssp')          % before the projection
            qb = pre.heart_ssp.cfa_before; lb = pre.heart_ssp.cfa_late_before;
        else
            qb = pre.cfa_before; lb = pre.cfa_late_before;
        end
        kept = NaN;
        if synthetic
            Rs = run_method(S, m{1}, common);
            kept = kept_share(Rs, R, roiLabels, gain);
        end
        nEp = NaN; if R.trials > 1, nEp = R.trials; end
        rows(end+1,:) = {[fn ext], m{1}, qb, pre.cfa_after, lb, pre.cfa_late_after, nEp, kept}; %#ok<AGROW>
        fprintf('%s %s: QRS %.2f -> %.2f uV, kept %.0f%%\n', fn, m{1}, qb, pre.cfa_after, kept);
    end
end
T = cell2table(rows, 'VariableNames', {'file' 'method' 'qrs_before' 'qrs_after' 'late_before' 'late_after' 'epochs' 'kept'});
T.qrs_reduction = 100 * (1 - T.qrs_after ./ T.qrs_before);
% Response kept relative to no removal: the average reference of the
% cleaning alone scales the injected response by 1 - nROI/nChannels
T.kept_vs_none = nan(height(T),1);
for f = unique(T.file)'
    i = strcmp(T.file, f{1}); k0 = T.kept(i & strcmp(T.method,'none'));
    if ~isempty(k0), T.kept_vs_none(i) = 100 * T.kept(i) / k0; end
end
disp(T)
if ~isempty(outfile)
    writetable(T, outfile);
end
end

function R = run_method(EEG, method, common)
% One BrainBeats run, the ICA starting the same way each time
rng(0, 'twister');
R = brainbeats_process(EEG, 'analysis','hep', common{:}, 'heart_removal',method);
end

function k = kept_share(Rs, R, roiLabels, gain)
% Share (%) of the synthetic response in the difference of the ROI averages
% (with minus without), 200 to 400 ms
roi = @(X) ismember(lower({X.chanlocs.labels}), lower(roiLabels));
d = mean(mean(Rs.data(roi(Rs),:,:),3),1) - mean(mean(R.data(roi(R),:,:),3),1);
tm = R.times;
ref = gain * exp(-(tm-300).^2/(2*40^2));
s = tm >= 200 & tm <= 400;
k = 100 * sum(d(s).*ref(s)) / sum(ref(s).^2);
end
