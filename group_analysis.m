%% GROUP ANALYSIS: HEP
% Download dataset (6 GB) here: https://nemar.org/dataexplorer/detail?dataset_id=ds003944
% Unizp the folder and copy-paste the path below
clear; close all; clc

% Copy path to data folder here
dataDir = 'C:\Users\Tracy\Downloads\ds003944';

% Import participants data
cd(dataDir)
tmp = readtable("participants.tsv","FileType","text",'Delimiter','\t');
ids = tmp.participant_id;
age = round(tmp.age);
grp = tmp.type;
nSub = length(ids);

% Remove participants with A to avoid duplicates
idx = contains(ids,'A');
ids(idx) = []; age(idx) = []; grp(idx) = [];
nSub = length(ids);

% Launch EEGLAB
eeglab; close;

% Install the bva-io EEGLAB plugin 
% Go to File > Manage eeglab extensions > click on bva-io > Install

% Path to electrode locations
locPath = fileparts(which('dipfitdefs.m'));

% Loop through files and process them for HEP
progressbar('Prossessing files for HEP analysis')
for iFile = 1:nSub

    id = ids{iFile};
    filepath = fullfile(dataDir, id, 'eeg');
    filename = sprintf('%s_task-Rest_eeg.vhdr',id);

    % Load .vhdr file
    EEG = pop_loadbv(filepath, filename, [], []);

    % Import electrode labels, type, and locations
    labels = readtable(fullfile(filepath, sprintf('%s_task-Rest_channels.tsv',id)),"FileType","text",'Delimiter','\t');
    % EEG.chanlocs = [];
    % EEG.chanlocs.labels = labels.name;
    for i = 1:size(labels,1)
        EEG.chanlocs(i).labels = labels.name{i};
        EEG.chanlocs(i).type = labels.type{i};
    end
    EEG = pop_chanedit(EEG,'lookup',fullfile(locPath,'standard_BEM','elec','standard_1005.elc'));

    % Remove unnecessary channels
    EEG = pop_select(EEG, 'rmchannel',{'VEOG' 'Misc' 'M2'});

    % Downsample to 125 Hz to accelerate things
    try
        EEG = pop_resample(EEG, 125);
    catch
        warning('This file has a strange sampling rate and fails. Skipping it.');
        continue
    end
    
    % Add filename and filepath if not already present
    EEG.filename = filename;
    EEG.filepath = filepath;
    
    % Fields to edit for STUDY
    EEG.group = grp{iFile};
    EEG.subject = id;

    % Process file with BrainBeats for HEP analysis
    % EEG = brainbeats_process(EEG,'analysis','hep','heart_signal','ECG', ...
    %     'heart_channels',{'ECG'},'clean_eeg',true,'save',true,...
    %     'vis_cleaning',false,'vis_outputs',true); 

    % Process and extract features for feature-based analysis
    EEG = brainbeats_process(EEG,'analysis','features','heart_signal','ECG', ...
        'heart_channels',{'ECG'},'clean_eeg',0,'norm',1, ...
        'hrv_features', {'time' 'frequency' 'nonlinear'}, ...
        'eeg_features', {'time' 'frequency' 'nonlinear'}, ...
        'gpu',0,'parpool',1,'vis_cleaning',0,'vis_outputs',1,'save',0);
    
    % Update progressbar
    progressbar(iFile/nSub)

end
disp('Done processing all files for HEP analysi!'); gong



%% Create STUDY

commands = {};
progressbar('Creating EEGLAB study')
for iSub = 1:nSub

    disp('--------------------------------------------')
    fprintf('Subject %g \n ', iSub)
    disp('--------------------------------------------')
    
    id = ids{iSub};
    filepath = fullfile(dataDir, id, 'eeg');
    filename = sprintf('%s_task-Rest_eeg_HEP.set',id);

    if exist(fullfile(filepath,filename),'file')
        % EEG = pop_loadset('filename',filename,'filepath',filepath);
        % pop_eegplot(EEG,1,1,1);
        
        commands = [ commands(:)' 'index' iSub 'load' fullfile(filepath, filename) ];
        [STUDY, ALLEEG] = std_editset(STUDY,ALLEEG,'name','HEP', ...
            'commands',commands,'updatedat','on','savedat','off','rmclust','off');
        [STUDY, ALLEEG] = std_checkset(STUDY, ALLEEG); 
        EEG = ALLEEG; 
    end


    % Store epoch length (in frames) to homogenize files to smallest later
    nFrames(iSub,:) = size(ALLEEG(iSub).times,2);

    % progressbar
    progressbar(iFile/nSub)

end

% Save
STUDY.filename = 'HEP.study';
STUDY.filepath = dataDir;
pop_savestudy(STUDY,EEG,'filename','study.study','filepath',dataDir);
EEG = ALLEEG; CURRENTSTUDY = 1; CURRENTSET = 1:length(EEG);
eeglab redraw

% Create design with group variable
STUDY = std_makedesign(STUDY,ALLEEG,1,'name','STUDY.design 1','delfiles','off','defaultdesign','off', ...
    'variable1','group','values1',{'Control','Psychosis'},'vartype1','categorical', ...
    'variable2','type','values2',{'R-peak'},'vartype2','categorical', ...
    'subjselect',{ALLEEG.subject});

% Precompute ERP
[STUDY, ALLEEG] = std_precomp(STUDY,ALLEEG,{},'savetrials','on',...
    'recompute','on','erp','on','erpparams',{'rmbase',[-200 -100]});
% [STUDY, ALLEEG] = std_precomp(STUDY,ALLEEG,{},'savetrials','on',...
%     'recompute','on','erp','on');

%% Run stats with LIMO

% Run LIMO 1st level
pop_limo(STUDY,ALLEEG,'method','WLS','measure','daterp', ...
    'timelim',[-100 500],'erase','on','splitreg','off');


%% Run statistics without GLM

low_frame = -100;
nboot = 1000;   % number of permutations for H0
alpha = 0.05;      % NaN to show maps of p-values
dpt = 'idpt';        % 'dpt' for dependent, 'idpt' for independent

% frames (ERP)
frames = ALLEEG(1).times;
[~, low_frame] = min(abs(frames - (low_frame))); % in ms;
high_frame = min(nFrames);   % smallest epoch length
frames = [low_frame high_frame];
times = ALLEEG(1).times(low_frame:high_frame);
nframes = length(low_frame:high_frame);

% number of channels
nChan = ALLEEG(1).nbchan;

% Average ERP across trials and with desired frames
X = nan(nChan,nframes,nSub);
for iSub = 1:nSub
    X(:,:,iSub) = trimmean(ALLEEG(iSub).data(:,low_frame:high_frame,:),20,3);
end

% Separate the two groups
grps = unique({ALLEEG.group});
idx = strcmpi({ALLEEG.group},grps(1));
x1 = X(:,:,idx);    % control group
x2 = X(:,:,~idx);   % psychosis group
x1_name = grps(1);
x2_name = grps(2);

% Run permutation statistics
[tvals,pvals,tvals_H0,pvals_H0] = compute_randomeffect(x1, x2, nboot, 'trimmed mean',dpt);
disp('Saving statistical outputs...')
save(fullfile(dataDir,'group_statistics.mat'),'tvals','pvals','tvals_H0','pvals_H0')

% Get channel neighbors for cluster correction for multiple comparison
chanlocs = ALLEEG(1).chanlocs;
[neighbors, neighbormatrix] = get_channelneighbors(chanlocs);

% Uncorrected
mcctype = 0;
[mask, pcorr, nClust] = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcctype, alpha, neighbormatrix);
plotresults('time', times, tvals, mask, pcorr, alpha, chanlocs, mcctype);
% saveas(gcf, fullfile(dataDir,'group-statistics_uncorrected.fig'));
% print(gcf, fullfile(dataDir,'group-statistics_uncorrected.png'),'-dpng','-r300');   % 300 dpi .png
% print(gcf, fullfile(dataDir,'group-statistics_uncorrected.tiff'),'-dtiff','-r300');  % 300 dpi .tiff

% Max-corrected
mcctype = 1;
[mask, pcorr, nClust] = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcctype, alpha, neighbormatrix);
plotresults('time', times, tvals, mask, pcorr, alpha, chanlocs, mcctype);
% saveas(gcf, fullfile(dataDir,'group-statistics_max-corrected.fig'));
% print(gcf, fullfile(dataDir,'group-statistics_max-corrected.png'),'-dpng','-r300');   % 300 dpi .png
% print(gcf, fullfile(dataDir,'group-statistics_max-corrected.tiff'),'-dtiff','-r300');  % 300 dpi .tiff

% Cluster-corrected
mcctype = 2;
[mask, pcorr, nClust] = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcctype, alpha, neighbormatrix);
plotresults('time', times, tvals, mask, pcorr, alpha, chanlocs, mcctype);
% saveas(gcf, fullfile(dataDir,'group-statistics_cluster-corrected.fig'));
% print(gcf, fullfile(dataDir,'group-statistics_cluster-corrected.png'),'-dpng','-r300');   % 300 dpi .png
% print(gcf, fullfile(dataDir,'group-statistics_cluster-corrected.tiff'),'-dtiff','-r300');  % 300 dpi .tiff

% TFCE-corrected
mcctype = 3;
[mask, pcorr, nClust] = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcctype, alpha, neighbormatrix);
plotresults('time', times, tvals, mask, pcorr, alpha, chanlocs, mcctype);
% saveas(gcf, fullfile(dataDir,'group-statistics_tfce-corrected.fig'));
% print(gcf, fullfile(dataDir,'group-statistics_tfce-corrected.png'),'-dpng','-r300');   % 300 dpi .png
% print(gcf, fullfile(dataDir,'group-statistics_tfce-corrected.tiff'),'-dtiff','-r300');  % 300 dpi .tiff

% % Plot mean + 95% HDI of peak effect
% elecName = 'Cz';
% elecNum = contains({chanlocs.labels},elecName);
% plotHDI(times,squeeze(x1(elecNum,:,:)), squeeze(x2(elecNum,:,:)), ...
%     'trimmed mean',0.05,[],'control','psychosis')


%% Meditation data

% Download the open-source dataset here (4.3 GB): 
% https://nemar.org/dataexplorer/detail?dataset_id=ds001787
% Unzip the folder in the destination of your choice, and add the path to the
% directory below. Find details on the study here: 
% https://cerco.cnrs.fr/wp-content/uploads/2019/04/brandmeyer_t_18_2519.pdf
% data_dir = 'copy-paste you path here';
data_dir = 'C:\Users\Tracy\Documents\MATLAB\data_meditation_tracy';
cd(data_dir)

% EEG = pop_biosig(fullfile('sub-001\ses-01\eeg\sub-001_ses-01_task-meditation_eeg.bdf'));

% Import STUDY using pop_importbids
[STUDY, ALLEEG] = pop_importbids( data_dir, 'eventtype','value', ...
    'bidsevent','on','bidschanloc','on', ...
    'outputdir',fullfile(data_dir,'derivatives') );
[STUDY, ALLEEG] = std_checkset(STUDY, ALLEEG);
EEG = ALLEEG; CURRENTSTUDY = 1; CURRENTSET = 1:length(EEG);

% Add session variable
% STUDY = std_makedesign(STUDY, ALLEEG, 1, 'name','STUDY.design 1','delfiles','off', ...
%     'defaultdesign','off','variable1','group','values1',{'expert','novice'}, ...
%     'vartype1','categorical','variable2','session','values2',{1,2,3},'vartype2', ...
%     'continuous','subjselect',{ALLEEG.subject});

% Loop through each file to process for HEP analysis
for iFile = 1:length(ALLEEG)
    EEG = pop_biosig(fullfile(ALLEEG(iFile).filepath, ALLEEG(iFile).filename));
end


% Extract relevant info from STUDY
ids = {ALLEEG.subject};
nSub = length(ids);  % n
filepaths = {ALLEEG.filepath};
filenames = {ALLEEG.filename};


design = STUDY.design;

% Process files for HEP analysis
for iSub = 1:nSub
    
    % Get subject data
    % EEG = pop_biosig(fullfile(filepath, filename));


    % Probes 
    % Q1: “Please rate the depth of your meditation,” for
    %   which participants evaluated the subjective depth of their
    %   “meditative state” for the moments immediately preceding
    %   the first probe, on a scale from 0 (not meditating at all) to
    %   3 (deep meditative state) by pressing the corresponding key
    %   on the keypad. 
    % Q2: “Please rate the depth of your mind wandering” automatically 
    %   followed. Participants evaluated the subjective depth of their 
    %   “mind wandering” for the period
    %   of time which immediately preceded the first probe, on a
    %   scale from 0 (not mind wandering at all) to 3 (immersed
    %   in their thoughts). 
    % Q3: "Please rate how tired you are,” where participants were asked to
    %   rate the subjective depth of their drowsiness at the time of the
    %   first question, from 0 (not drowsy at all) to 3 (very drowsy).

    % Markers
    % "2": "Response 1 (this may be a response to question 1, 2 or 3)",
    % "4": "Response 2 (this may be a response to question 1, 2 or 3)",
    % "8": "Response 3 (this may be a response to question 1, 2 or 3)",
    % "16": "Indicate involuntary response",
    % "128": "First question onset (most important marker)"

    % Load file
    filepath = fullfile(data_dir, subjects{iSub}, 'ses-01','eeg');
    % filename = sprintf('%s', subjects{iSub}, 'ses-01', );
    % sub-001_ses-01_task-meditation_eeg.bdf


    % Data were then segmented into 10 s-epochs from −10.05 to −0.05 s 
    % prior to onset of Q1 in the experience-sampling probe series.
    
    
    % You will be prompted to install the biosig plugin if not already installed
    % to load the .bdf files


end


