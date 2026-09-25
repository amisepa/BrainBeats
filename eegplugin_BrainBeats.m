% EEGPLUGIN_BRAINBEATS - EEGLAB plugin entry point for BrainBeats.
%
% Adds the plugin folders to the MATLAB path and creates the
% Tools > BrainBeats menu, which opens the brainbeats_process GUI.
% Called by EEGLAB at startup.
%
% Usage:
%   vers = eegplugin_BrainBeats(fig, try_strings, catch_strings);
%
% Inputs:
%   fig           - EEGLAB main figure handle
%   try_strings   - EEGLAB try strings for menu callbacks
%   catch_strings - EEGLAB catch strings for menu callbacks
%
% Output:
%   vers          - plugin version (string)
%
% Copyright (C) - Cedric Cannard, 2023
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function vers = eegplugin_BrainBeats(fig,try_strings,catch_strings)

% Plugin version
vers = '1.6';

% Add paths to subfolders
p = fileparts(which('eegplugin_BrainBeats.m'));
addpath(p);
addpath(fullfile(p,'functions'))
addpath(fullfile(p,'sample_data'))

% Find the EEGLAB Tools menu
menu = findobj(fig, 'tag', 'tools');

% Menu callback (runs brainbeats_process on the current dataset and stores
% the command in the EEGLAB history)
process = [try_strings.no_check '[EEG, LASTCOM] = brainbeats_process(EEG);' catch_strings.new_and_hist];

% Create menus
submenu = uimenu(menu, 'Label', 'BrainBeats', 'separator', 'on');
uimenu(submenu, 'Label', '1st level (subject)', 'CallBack', process);
