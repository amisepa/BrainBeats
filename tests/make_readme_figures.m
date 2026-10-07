function make_readme_figures(runs)
% MAKE_README_FIGURES - Regenerate the README figures from the sample data.
%
% Saves the BrainBeats window, then runs HEP (ECG, and PPG with the pulse arrival time corrected), features, rm_heart and
% coherence on sample_data/dataset.set with the plots on, and saves the
% figures shown in README.md into the figures folder (existing files are
% overwritten). It also serves as the interactive check of the plotting code
% before a release.
%
% Usage:
%   make_readme_figures                          % all figures
%   make_readme_figures({'hep_ecg' 'features'})  % a subset of the runs
%
% Input:
%   runs - any of 'gui', 'hep_ecg', 'hep_ppg', 'features',
%          'rm_heart', 'coherence' (default: all)
%
% Notes:
% - Run it in MATLAB with the desktop: drawing can hang under 'matlab -batch'.
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
    runs = {'gui' 'hep_ecg' 'hep_ppg' 'features' 'rm_heart' 'coherence'};
end
runs = cellstr(runs);
common = {'save',false,'gong',false};

for r = 1:numel(runs)
    close all force
    fprintf('\n==================== %s ====================\n', runs{r});
    EEG = pop_loadset('filename','dataset.set','filepath',dpath);

    switch runs{r}

        case 'gui'
            % The main window, with the sample dataset
            setappdata(0,'brainbeats_main_snapshot',fullfile(outdir,'gui_main.png'));
            try
                brainbeats_main_window(EEG);
            catch ME
                warning('Main window: %s', ME.message)
            end
            rmappdata(0,'brainbeats_main_snapshot');
            fprintf('Saved %s\n', fullfile(outdir,'gui_main.png'));

        case 'hep_ecg'
            % Cleaned EEG, default epochs (-300 to 600 ms), frontocentral ROI,
            % 100 surrogate heartbeat trains. vis_cleaning on for the R-peak plot.
            EEG = pop_select(EEG,'nochannel',{'PPG'});
            brainbeats_process(EEG,'analysis','hep','heart_signal','ECG', ...
                'heart_channels',{'ECG'},'clean_eeg',true,'icamethod',1, ...
                'hep_surrogates',100,'vis_cleaning',true,'vis_outputs',true,common{:});
            save_fig(byname('Time series & corresponding RR/NN intervals',1), outdir, 'heartbeats_ecg.png', 'export');
            save_fig(byname('HEP',1), outdir, 'hep_ecg.png', 'export');
            save_fig(byname('HEP time-frequency',1), outdir, 'hep_tf.png', 'export');

        case 'hep_ppg'
            % PPG beats shifted back by the pulse arrival time, estimated
            % from the EEG cardiac field artifact (ECG removed), same
            % cleaning as the ECG run
            EEG = pop_select(EEG,'nochannel',{'ECG'});
            brainbeats_process(EEG,'analysis','hep','heart_signal','PPG', ...
                'heart_channels',{'PPG'},'clean_eeg',true,'icamethod',1, ...
                'vis_cleaning',false,'vis_outputs',true,common{:});
            save_fig(byname('HEP',1), outdir, 'hep_ppg.png', 'export');

        case 'features'
            % Both feature figures end up named 'Visualization of features':
            % the PSD figure is created first, the scalp topographies second
            EEG = pop_select(EEG,'nochannel',{'PPG'});
            brainbeats_process(EEG,'analysis','features','heart_signal','ECG', ...
                'heart_channels',{'ECG'},'clean_eeg',true,'linenoise',50, ...
                'vis_cleaning',false,'vis_outputs',true,common{:});
            save_fig(byname('Visualization of features',1), outdir, 'features_psd.png', 'export');
            save_fig(byname('Visualization of features',2), outdir, 'features_eeg.png', 'export');

        case 'rm_heart'
            % Same settings as the tutorial (Picard to save time)
            EEG = pop_select(EEG,'nochannel',{'PPG'});
            brainbeats_process(EEG,'analysis','rm_heart','heart_signal','ECG', ...
                'heart_channels',{'ECG'},'clean_eeg',true,'linenoise',50, ...
                'conf_thresh',.65,'icamethod',1,'keep_heart',true, ...
                'vis_cleaning',false,'vis_outputs',true,common{:});
            save_fig(byname('Heart component(s) removed',1), outdir, 'rm_heart_components.png', 'frame');
            save_fig(byname('Heart components removed',1), outdir, 'rm_heart_data.png', 'frame');

        case 'coherence'
            % Unnamed figures, in creation order: spectra of the 4 measures,
            % then coherence, partial, directed and partial directed topos
            EEG = pop_select(EEG,'nochannel',{'PPG'});
            brainbeats_process(EEG,'analysis','coherence','heart_signal','ECG', ...
                'heart_channels',{'ECG'},'clean_eeg',false,'parpool',false, ...
                'vis_cleaning',false,'vis_outputs',true,common{:});
            figs = sortfigs(findall(groot,'Type','figure'));
            if numel(figs) >= 2
                save_fig(figs(1), outdir, 'coherence_spectra.png', 'export');
                save_fig(figs(2), outdir, 'coherence_bands.png', 'export');
            else
                warning('Expected at least 2 coherence figures, found %d.', numel(figs))
            end

        otherwise
            warning('Unknown run ''%s'' (skipped).', runs{r})
    end
end
fprintf('\nDone. Check the figures in %s, then README.md.\n', outdir);


function fig = byname(name, k)
% k-th figure (in creation order) with this Name; empty if none
figs = sortfigs(findall(groot,'Type','figure','Name',name));
if numel(figs) >= k, fig = figs(k); else, fig = []; end


function figs = sortfigs(figs)
[~, order] = sort([figs.Number]);
figs = figs(order);


function save_fig(fig, outdir, file, mode)
% mode 'export': axes content at 150 dpi; 'frame': the figure as on screen,
% for EEGLAB figures made of buttons and scrollable axes
if isempty(fig) || ~isgraphics(fig,'figure')
    warning('Figure for %s not found (not saved).', file); return
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
