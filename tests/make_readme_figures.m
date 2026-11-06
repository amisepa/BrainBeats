function make_readme_figures(runs)
% MAKE_README_FIGURES - Regenerate the README figures from the sample data.
%
% Saves the BrainBeats window, then runs HEP (ECG, and PPG with the pulse
% arrival time corrected; channels and ICs) on sample_data/sub-047_rest.set
% (a clearer heartbeat-evoked response than the tutorial dataset) and
% features on sample_data/dataset.set, with the plots on, and saves the
% figures shown in README.md into the figures folder (existing files are
% overwritten). It also serves as the interactive check of the plotting code
% before a release.
%
% Usage:
%   make_readme_figures                          % all figures
%   make_readme_figures({'hep_ecg' 'features'})  % a subset of the runs
%
% Input:
%   runs - any of 'gui', 'hep_ecg', 'hep_ecg_ppg' (ROI HEP from the ECG and
%          the PPG, with and without the pulse arrival time), 'ics_ecg',
%          'ics_ppg' (IC level), 'features' (default: all)
%
% Notes:
% - Works with the desktop or under 'matlab -batch' (figures then drawn
%   invisible, except the main window; MATLAB may not exit after drawing:
%   wrap it in a timeout, and check the figures were saved).
% - Nothing is saved next to the sample data ('save' is off), and the working
%   directory is tempdir during the runs (the 3D headplot writes a spline file).
% - Any other copy of BrainBeats on the path (e.g. eeglab/plugins/BrainBeats*)
%   is removed so this repository's code is what runs.
%
% Runtime: ~15-20 min (ICA, surrogate heartbeat trains, nonlinear EEG
% features).
%
% Copyright (C) - Cedric Cannard, 2026

% Paths: EEGLAB, then this repository instead of any other BrainBeats copy
bbroot = fileparts(fileparts(mfilename('fullpath')));
if ~exist('pop_loadset','file'), eeglab; close; end
p = strsplit(path, pathsep);
old = p(contains(p,'BrainBeats','IgnoreCase',true) & ~startsWith(p,bbroot));
if ~isempty(old), rmpath(old{:}); end
addpath(bbroot, fullfile(bbroot,'functions'));
fprintf('Using %s\n', which('brainbeats_process'));

outdir = fullfile(bbroot,'figures');
dpath = fullfile(bbroot,'sample_data');
olddir = cd(tempdir);
restore = onCleanup(@() cd(olddir));

if nargin < 1 || isempty(runs)
    runs = {'gui' 'hep_ecg' 'hep_ecg_ppg' 'ics_ecg' 'ics_ppg' 'features'};
end
runs = cellstr(runs);
common = {'save',false,'gong',false};

% Under 'matlab -batch', figures are drawn invisible: they would otherwise
% open on the screen, where closing one loses it (exportgraphics does not
% need them visible)
if exist('batchStartupOptionUsed') > 0 && batchStartupOptionUsed %#ok<EXIST>
    vis = get(groot,'DefaultFigureVisible');
    set(groot,'DefaultFigureVisible','off');
    restoreVis = onCleanup(@() set(groot,'DefaultFigureVisible',vis));
end

for r = 1:numel(runs)
    close all force
    fprintf('\n==================== %s ====================\n', runs{r});
    if any(strcmp(runs{r}, {'hep_ecg' 'hep_ecg_ppg' 'ics_ecg' 'ics_ppg'}))
        EEG = pop_loadset('filename','sub-047_rest.set','filepath',dpath);
    else
        EEG = pop_loadset('filename','dataset.set','filepath',dpath);
    end

    switch runs{r}

        case 'gui'
            % The main window, with the sample dataset (drawn on the screen:
            % the snapshot is a screen capture)
            setappdata(0,'brainbeats_main_snapshot',fullfile(outdir,'gui_main.png'));
            vis0 = get(groot,'DefaultFigureVisible');
            set(groot,'DefaultFigureVisible','on');
            try
                brainbeats_main_window(EEG);
            catch ME
                warning('Main window: %s', ME.message)
            end
            set(groot,'DefaultFigureVisible',vis0);
            rmappdata(0,'brainbeats_main_snapshot');
            fprintf('Saved %s\n', fullfile(outdir,'gui_main.png'));

        case 'hep_ecg'
            % Cleaned EEG, default epochs (-300 to 600 ms), frontocentral ROI,
            % 100 surrogate heartbeat trains. vis_cleaning on for the R-peak plot.
            EEG = pop_select(EEG,'nochannel',{'PPG'});
            H = brainbeats_process(EEG,'analysis','hep','heart_signal','ECG', ...
                'heart_channels',{'ECG'},'clean_eeg',true,'icamethod',1, ...
                'hep_surrogates',100,'vis_cleaning',true,'vis_outputs',true,common{:});
            % numbers of the README caption
            S = H.brainbeats.roi.surrogate; t = H.brainbeats.roi.hep_times; p = ['p_' S.stat];
            fprintf('README: ROI HEP differs from the surrogates (%s) at %s ms; HRSP (induced) points %d, HRPC points %d\n', ...
                S.stat, spans(t, S.hep.(p)(1,:) < .05), sum(S.hrsp_induced.(p)(:) < .05), sum(S.hrpc.(p)(:) < .05));
            save_fig(byname('Time series & corresponding RR/NN intervals',1), outdir, 'heartbeats_ecg.png', 'export');
            save_fig(byname('HEP',1), outdir, 'hep_ecg.png', 'export');
            save_fig(byname('HEP time-frequency',1), outdir, 'hep_tf.png', 'export');
            save_fig(byname('QRS map (cardiac field artifact)',1), outdir, 'hep_cfa.png', 'export');

        case 'hep_ecg_ppg'
            % ROI HEP time-locked to the ECG R-peaks, to the PPG pulse onsets
            % shifted back by the pulse arrival time measured from the ECG of
            % the file ('ppg_transit','auto'), and to the PPG pulse onsets
            % themselves ('off'); same cleaning for the three
            opts = {'clean_eeg',true,'icamethod',1,'vis_cleaning',false,'vis_outputs',false,common{:}};
            H = cell(1,3);
            H{1} = brainbeats_process(pop_select(EEG,'nochannel',{'PPG'}),'analysis','hep', ...
                'heart_signal','ECG','heart_channels',{'ECG'},opts{:});
            H{2} = brainbeats_process(EEG,'analysis','hep','heart_signal','PPG','heart_channels',{'PPG'},opts{:});
            H{3} = brainbeats_process(EEG,'analysis','hep','heart_signal','PPG','heart_channels',{'PPG'}, ...
                'ppg_transit','off',opts{:});
            pat = H{2}.brainbeats.preprocessings.ppg_transit.median;
            t = H{1}.brainbeats.roi.hep_times;
            y = cellfun(@(h) h.brainbeats.roi.hep(1,:), H, 'UniformOutput', false);
            rp = cellfun(@(v) corr(v(:), y{1}(:)), y(2:3));
            fig = figure('Color','w','Name','HEP from the ECG and the PPG','Position',[100 100 900 380]);
            hold on
            plot(t, y{1}, 'k', 'LineWidth', 2.5);
            plot(t, y{2}, 'Color', [0 .45 .74], 'LineWidth', 2);
            plot(t, y{3}, '--', 'Color', [.85 .33 .1], 'LineWidth', 2);
            xline(0, 'k:'); yline(0, 'k:'); box on; xlim(t([1 end]))
            xlabel('Latency (ms)'); ylabel('HEP, frontocentral ROI (\muV)')
            legend({'ECG R-peaks', sprintf('PPG pulse onsets - %.0f ms (pulse arrival time from the ECG; r = %.2f)', pat, rp(1)), ...
                sprintf('PPG pulse onsets, not shifted (r = %.2f)', rp(2))}, 'Location','southwest', 'Box','off')
            title('Same EEG, time-locked to the ECG, or to the PPG with and without the pulse arrival time')
            set(gca,'FontSize',11)
            save_fig(fig, outdir, 'hep_ecg_ppg.png', 'export');

        case {'ics_ecg' 'ics_ppg'}
            % IC level: map and HEP against 100 surrogate trains of each IC
            % left after the cleaning (first figure: up to 24 ICs)
            if strcmp(runs{r},'ics_ecg')
                EEG = pop_select(EEG,'nochannel',{'PPG'});
                heart = {'heart_signal','ECG','heart_channels',{'ECG'}};
            else
                heart = {'heart_signal','PPG','heart_channels',{'PPG'}};
            end
            brainbeats_process(EEG,'analysis','hep',heart{:},'clean_eeg',true,'icamethod',1, ...
                'hep_level','ics','hep_surrogates',100,'vis_cleaning',false,'vis_outputs',true,common{:});
            f = sortfigs(findall(groot,'Type','figure','-regexp','Name','^HEP of the independent components'));
            if isempty(f), f = []; else, f = f(1); end
            save_fig(f, outdir, sprintf('hep_%s.png', runs{r}), 'export');

        case 'features'
            % Both feature figures end up named 'Visualization of features':
            % the PSD figure is created first, the scalp topographies second
            EEG = pop_select(EEG,'nochannel',{'PPG'});
            brainbeats_process(EEG,'analysis','features','heart_signal','ECG', ...
                'heart_channels',{'ECG'},'clean_eeg',true,'linenoise',50, ...
                'vis_cleaning',false,'vis_outputs',true,common{:});
            save_fig(byname('Visualization of features',1), outdir, 'features_psd.png', 'export');
            save_fig(byname('Visualization of features',2), outdir, 'features_eeg.png', 'export');

        otherwise
            warning('Unknown run ''%s'' (skipped).', runs{r})
    end
end
fprintf('\nDone. Check the figures in %s, then README.md.\n', outdir);


function fig = byname(name, k)
% k-th figure (in creation order) with this Name; empty if none
figs = sortfigs(findall(groot,'Type','figure','Name',name));
if numel(figs) >= k, fig = figs(k); else, fig = []; end


function s = spans(t, sig)
% Latency intervals where sig is true, as text
d = diff([false sig(:)' false]); on = find(d == 1); off = find(d == -1) - 1;
s = strjoin(arrayfun(@(a,b) sprintf('%g-%g', t(a), t(b)), on, off, 'UniformOutput', false), ', ');
if isempty(s), s = 'none'; end


function figs = sortfigs(figs)
if isempty(figs), figs = []; return, end
[~, order] = sort([figs.Number]);
figs = figs(order);


function save_fig(fig, outdir, file, mode)
% mode 'export': axes content at 150 dpi; 'frame': the figure as on screen,
% for EEGLAB figures made of buttons and scrollable axes
if isempty(fig) || ~isgraphics(fig,'figure')
    names = arrayfun(@(h) h.Name, findall(groot,'Type','figure'), 'UniformOutput', false);
    warning('Figure for %s not found (not saved). Open figures: %s', file, strjoin(names,' | ')); return
end
out = fullfile(outdir, file);
figure(fig); drawnow
if strcmp(mode,'frame')
    F = getframe(fig);
    imwrite(F.cdata, out);
else
    exportgraphics(fig, out, 'Resolution', 150, 'BackgroundColor', 'white');
end
fprintf('Saved %s\n', out);
