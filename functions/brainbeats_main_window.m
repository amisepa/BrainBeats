% BRAINBEATS_MAIN_WINDOW - Main BrainBeats window (EEGLAB menu 'BrainBeats').
%
% Shows the dataset (or loads one, or the sample dataset), and selects the
% analysis, the heart signal and channel(s), for HEP whether the measures are
% computed on the scalp channels, the independent components or both, and the
% plots and saving. 'Next' opens the parameters of the analysis (GETPARAMS_GUI).
%
% Usage:
%   [params, EEG, abort] = brainbeats_main_window(EEG)
%
% Input:
%   EEG    - current EEGLAB dataset (can be empty: load one in the window)
%
% Outputs:
%   params - analysis, heart_signal, heart_channels, hep_level, vis_cleaning,
%            vis_outputs, save
%   EEG    - the dataset (loaded in the window if it was empty)
%   abort  - true if the window was cancelled or closed
%
% For automated tests, setappdata(0,'brainbeats_main_answers',S) makes the
% window return S (a params struct) without opening, and
% setappdata(0,'brainbeats_main_snapshot',file) saves an image of the window
% (with the current dataset) and closes it.
%
% Copyright (C) - Cedric Cannard, 2026

function [params, EEG, abort] = brainbeats_main_window(EEG)

params = struct(); abort = true;
if nargin < 1, EEG = []; end
S = getappdata(0,'brainbeats_main_answers');
if ~isempty(S)
    params = S; abort = false;
    return
end

bbroot = fileparts(which('eegplugin_BrainBeats.m'));
[bg, fg] = eeglab_colors();
analyses = {'HEP, HRSP, HRPC: heartbeat-evoked potentials, spectral perturbations, phase consistency' ...
    'EEG and heart-rate variability (HRV) features' 'Brain-heart coherence (beta)'};
analysisKeys = {'hep' 'features' 'coherence'};
heartTypes = {'ECG' 'PPG'};
domains = {'Scalp channels' 'Independent components' 'Both'};

% Not modal: the channel list and file dialogs opened from here must stay usable
W = 780; H = 590;
scr = get(0,'ScreenSize');
f = figure('Name','BrainBeats','NumberTitle','off','MenuBar','none','ToolBar','none', ...
    'Resize','off','Color',bg,'WindowStyle','normal', ...
    'Position',[(scr(3)-W)/2 max(40,(scr(4)-H)/2) W H], 'CloseRequestFcn',@(~,~) uiresume(gcbf));

% Header: logo, name, one-line description
logo = fullfile(bbroot,'brainbeats_logo2.png');
if isfile(logo)
    try
        ax = axes('Parent',f,'Units','pixels','Position',[20 H-130 115 115]);
        image(ax, load_logo(logo, bg, 230)); axis(ax,'image','off');
    catch
    end
end
txt(f, [150 H-65 600 36], 'BrainBeats', 22, 'bold', fg, bg);
txt(f, [150 H-95 600 24], sprintf('Joint analysis of EEG and heart signals (ECG, PPG)  -  v%s', eegplugin_BrainBeats), 11, 'normal', fg, bg);
txt(f, [150 H-120 600 22], 'Cannard, Wahbeh & Delorme (2024). Journal of Visualized Experiments, 206.', 9, 'normal', [.3 .3 .3], bg);

% Data
y = H-165;
heading(f, y, W, 'Data', fg, bg);
hData = txt(f, [30 y-52 W-340 40], '', 10, 'normal', [0 0 0], bg);
uicontrol(f,'Style','pushbutton','String','Load dataset...','Position',[W-290 y-45 130 30],'Callback',@loadData);
uicontrol(f,'Style','pushbutton','String','Sample dataset','Position',[W-150 y-45 120 30],'Callback',@loadSample);

% Analysis
y = y - 85;
heading(f, y, W, 'Analysis', fg, bg);
txt(f, [30 y-42 170 22], 'Analysis:', 10, 'bold', [0 0 0], bg);
hAna = uicontrol(f,'Style','popupmenu','String',analyses,'Position',[200 y-40 W-230 24],'Callback',@refresh);
txt(f, [30 y-77 170 22], 'Heart signal:', 10, 'bold', [0 0 0], bg);
hSig = uicontrol(f,'Style','popupmenu','String',heartTypes,'Position',[200 y-75 120 24],'Callback',@autoChan);
txt(f, [345 y-77 100 22], 'Channel(s):', 10, 'bold', [0 0 0], bg);
hChan = uicontrol(f,'Style','edit','String','','Position',[450 y-75 240 24],'HorizontalAlignment','left');
uicontrol(f,'Style','pushbutton','String','...','Position',[700 y-75 50 24],'Callback',@pickChan);
hLevTxt = txt(f, [30 y-112 170 22], 'Perform analysis on:', 10, 'bold', [0 0 0], bg);
hIcs = uicontrol(f,'Style','popupmenu','String',domains,'Position',[200 y-110 W-230 24],'Callback',@refresh);
hInfo = txt(f, [30 y-200 W-60 78], '', 9, 'normal', [.25 .25 .25], bg);

% Outputs
y = y - 215;
heading(f, y, W, 'Outputs', fg, bg);
hVisC = uicontrol(f,'Style','checkbox','String','Plot the preprocessing steps','Value',1, ...
    'Position',[30 y-40 240 24],'BackgroundColor',bg,'FontSize',10);
hVisO = uicontrol(f,'Style','checkbox','String','Plot the outputs','Value',1, ...
    'Position',[290 y-40 200 24],'BackgroundColor',bg,'FontSize',10);
hSave = uicontrol(f,'Style','checkbox','String','Save the outputs','Value',1, ...
    'Position',[500 y-40 200 24],'BackgroundColor',bg,'FontSize',10);

% Buttons
uicontrol(f,'Style','pushbutton','String','Help','Position',[20 15 90 32], ...
    'Callback',@(~,~) pophelp('brainbeats_process'));
uicontrol(f,'Style','pushbutton','String','Cancel','Position',[W-310 15 130 32], ...
    'Callback',@(~,~) uiresume(f));
uicontrol(f,'Style','pushbutton','String','Next  >','FontWeight','bold','Position',[W-165 15 145 32], ...
    'Callback',@next);

showData(); autoChan(); refresh();
try finish_figure(f); catch; end

% Screenshot for the documentation (setappdata(0,'brainbeats_main_snapshot',file))
snap = getappdata(0,'brainbeats_main_snapshot');
if ~isempty(snap)
    drawnow; pause(0.5);
    F = getframe(f); imwrite(F.cdata, snap);
    delete(f); abort = true;
    return
end
uiwait(f);
if isgraphics(f), delete(f); end


    function next(~,~)
        if isempty(EEG) || ~isfield(EEG,'data') || isempty(EEG.data)
            errordlg('Load a dataset first.','BrainBeats'); return
        end
        chans = strsplit(strtrim(get(hChan,'String')));
        chans = chans(~cellfun(@isempty, chans));
        if isempty(chans)
            errordlg('Select the ECG/PPG channel(s).','BrainBeats'); return
        end
        miss = chans(~ismember(lower(chans), lower({EEG.chanlocs.labels})));
        if ~isempty(miss)
            errordlg(sprintf('Channel(s) not in the dataset: %s', strjoin(miss,', ')),'BrainBeats'); return
        end
        params.analysis = analysisKeys{get(hAna,'Value')};
        params.heart_signal = lower(heartTypes{get(hSig,'Value')});
        params.heart_channels = chans;
        if strcmp(params.analysis,'hep')
            levels = {'channels' 'ics' 'both'};
            params.hep_level = levels{get(hIcs,'Value')};
        end
        params.vis_cleaning = logical(get(hVisC,'Value'));
        params.vis_outputs = logical(get(hVisO,'Value'));
        params.save = logical(get(hSave,'Value'));
        abort = false;
        uiresume(f);
    end

    function refresh(~,~)
        isHEP = get(hAna,'Value') == 1;
        onoff = {'off' 'on'};
        set([hLevTxt hIcs], 'Enable', onoff{isHEP+1});
        switch get(hAna,'Value')
            case 1
                s = ['Epochs around the heartbeats: heartbeat-evoked potentials (HEP, time domain), heartbeat-' ...
                    'related spectral perturbations (HRSP, time-frequency power) and heartbeat-related phase ' ...
                    'consistency (HRPC, pairwise phase consistency across heartbeats), with the surrogate ' ...
                    'control. Heart components can be removed from the EEG with the cleaning (next window).'];
                if get(hIcs,'Value') ~= 2
                    s = [s ' Channels: all channels, and plots of the frontocentral ROI (can be changed in ' ...
                        'the next window).'];
                end
                if get(hIcs,'Value') ~= 1
                    s = [s ' Components: every independent component, exported with its map and ICLabel ' ...
                        'class, without selecting any (heartbeat-locked components can hold cardiac field ' ...
                        'artifact). Needs ICA (with the EEG cleaning, or run here).'];
                end
            case 2
                s = 'HRV (time, frequency, nonlinear) and EEG (time, frequency, nonlinear) features of the continuous data.';
            case 3
                s = ['Coherence, partial and directed coherence between each EEG channel and the heart ' ...
                    '(NN intervals), from a multivariate autoregressive model. Beta: not validated yet.'];
        end
        set(hInfo,'String',s);
    end

    function autoChan(~,~)
        if isempty(EEG) || ~isfield(EEG,'chanlocs') || isempty(EEG.chanlocs), return; end
        labs = {EEG.chanlocs.labels};
        if get(hSig,'Value') == 1
            m = labs(contains(lower(labs), {'ecg' 'ekg'}));
        else
            m = labs(contains(lower(labs), {'ppg' 'pleth' 'pulse' 'bvp'}));
        end
        set(hChan,'String',strjoin(m,' '));
    end

    function pickChan(~,~)
        if isempty(EEG) || ~isfield(EEG,'chanlocs') || isempty(EEG.chanlocs), return; end
        [~, sel] = pop_chansel({EEG.chanlocs.labels}, 'withindex', 'on');
        if ~isempty(sel), set(hChan,'String',sel); end
    end

    function loadData(~,~)
        try
            E = pop_loadset();
            if ~isempty(E) && isfield(E,'data') && ~isempty(E.data)
                EEG = E; showData(); autoChan();
            end
        catch ME
            errordlg(ME.message,'BrainBeats');
        end
    end

    function loadSample(~,~)
        EEG = pop_loadset('filename','dataset.set','filepath',fullfile(bbroot,'sample_data'));
        showData(); autoChan();
    end

    function showData()
        if isempty(EEG) || ~isfield(EEG,'data') || isempty(EEG.data)
            set(hData,'String',sprintf('No dataset loaded.\nLoad a dataset with EEG and ECG/PPG channels, or the sample dataset.'));
        else
            name = EEG.setname; if isempty(name), name = EEG.filename; end
            set(hData,'String',sprintf('%s\n%d channels, %g Hz, %.1f min%s', name, EEG.nbchan, EEG.srate, ...
                EEG.xmax*EEG.trials/60, ternary(EEG.trials > 1, sprintf(', %d epochs', EEG.trials), '')));
        end
    end
end

function h = txt(f, pos, s, sz, wt, col, bg)
h = uicontrol(f,'Style','text','String',s,'Position',pos,'FontSize',sz,'FontWeight',wt, ...
    'ForegroundColor',col,'BackgroundColor',bg,'HorizontalAlignment','left');
end

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end

function [bg, fg] = eeglab_colors()
% EEGLAB window colors (icadefs is a script: it cannot run in the main
% function, whose workspace is static because of its nested functions)
bg = [.93 .96 1]; fg = [0 0 .4];
try
    icadefs;
    bg = GUIBACKCOLOR; fg = GUITEXTCOLOR;
catch
end
end

function heading(f, y, W, label, fg, bg)
% Section heading: bold label over a thin line
uicontrol(f,'Style','text','String',label,'Position',[15 y W-30 20],'FontSize',10,'FontWeight','bold', ...
    'ForegroundColor',fg,'BackgroundColor',bg,'HorizontalAlignment','left');
uicontrol(f,'Style','text','String','','Position',[15 y-2 W-30 1],'BackgroundColor',fg);
end

function img = load_logo(file, bg, px)
% Logo as an RGB image of about px pixels: cropped to its opaque part, blended
% with the window color, and reduced by block averaging (a large image drawn
% in a small axes, or its transparency, is rendered poorly)
[img, map, alpha] = imread(file);
if ~isempty(map), img = ind2rgb(img, map); end
if isinteger(img), img = double(img) / double(intmax(class(img))); else, img = double(img); end
if size(img,3) == 1, img = repmat(img, 1, 1, 3); end
if isempty(alpha)
    alpha = ones(size(img,1), size(img,2));
else
    alpha = double(alpha) / double(intmax(class(alpha)));
    [r, c] = find(alpha > .05);
    img = img(min(r):max(r), min(c):max(c), :);
    alpha = alpha(min(r):max(r), min(c):max(c));
end
img = img .* alpha + reshape(bg, 1, 1, 3) .* (1 - alpha);
k = max(1, floor(max(size(img,1), size(img,2)) / px));
h = floor(size(img,1)/k); w = floor(size(img,2)/k);
img = reshape(img(1:h*k, 1:w*k, :), k, h, k, w, 3);
img = min(max(squeeze(mean(mean(img, 1), 3)), 0), 1);
end
