%% Prep BDF data for BrainBEats tutorial
clear; close all; clc
dataDir = 'G:\Shared drives\Grants\Post Award Grants\(736) Bial Full-trance 2017\Research\Data\EEG\BDF_files';
cd 'C:\Users\Tracy\Documents\MATLAB\BrainBeats'
outputDir = 'sample_data';
% cd(outputDir)
eeglab; close;
chanlocpath = fileparts(which('dipfitdefs.m'));

for iSub = 1:13

    disp('---------------------------------------------------------------')
    fprintf('                      File %g/13 \n', iSub)
    disp('---------------------------------------------------------------')

    filename = sprintf('subj%2.2d_1.bdf',iSub);  % only session 1
    EEG = pop_biosig(fullfile(dataDir, filename));
    EEG = pop_select(EEG,'rmchannel',{'EXG1','EXG2','EXG3','EXG4','EXG7','EXG8','GSR1','GSR2','Erg1','Erg2','Resp','Plet','Temp'});
    EEG = pop_chanedit(EEG,'rplurchanloc',1,'lookup',fullfile(chanlocpath,'standard_BEM','elec','standard_1005.elc'));
    EEG = fix_events(EEG, filename);
    % EEG.event = rmfield(EEG.event, 'edftype');
    EEG = pop_resample(EEG,256);

    % EEG = rm_DC(EEG);
    EEG.nbchan = EEG.nbchan+1;
    EEG.data(end+1,:) = zeros(1, EEG.pnts);
    EEG.chanlocs(1,EEG.nbchan).labels = 'initialReference';
    EEG = pop_reref(EEG, []);
    EEG = pop_select(EEG,'nochannel',{'initialReference'});
    dataRank = sum(eig(cov(double(EEG.data(:,:)'))) > 1E-7)

    % MW 1
    idx = find(contains({EEG.event.type}, {'rest1_start' 'rest1_end'}));
    if iSub == 1
        lats = [ EEG.event(idx(1)).latency+EEG.srate*10 EEG.event(idx(2)).latency-EEG.srate*300 ];
    else
        lats = [ EEG.event(idx(1)).latency+EEG.srate*10 EEG.event(idx(2)).latency-EEG.srate*5 ];
    end
    MW1 = pop_select(EEG, 'point', lats);
    fprintf('MW 1 length = %g min \n', MW1.xmax/60);

    % MW 2
    idx = find(contains({EEG.event.type}, {'rest2_start' 'rest2_end'}));
    if iSub == 1
        lats = [ EEG.event(idx(1)).latency+EEG.srate*10 EEG.event(idx(2)).latency-EEG.srate*180 ];
    else
        lats = [ EEG.event(idx(1)).latency+EEG.srate*10 EEG.event(idx(2)).latency-EEG.srate*5 ];
    end
    MW2 = pop_select(EEG, 'point', lats);
    fprintf('MW 2 length = %g min \n', MW2.xmax/60);
    
    % MW 3
    idx = find(contains({EEG.event.type}, {'rest3_start' 'rest3_end'}));
    lats = [ EEG.event(idx(1)).latency+EEG.srate*10 EEG.event(idx(2)).latency-EEG.srate*5 ];
    MW3 = pop_select(EEG, 'point', lats);
    fprintf('MW 2 length = %g min \n', MW3.xmax/60);
    
    % Merge MW trials and save 
    MW = pop_mergeset(MW1, MW2);
    MW = pop_mergeset(MW, MW3);
    MW.subject = sprintf('sub-%2.2d',iSub);
    MW.condition = 'mindwandering';
    MW.saved = 'no';
    MW = eeg_checkset(MW);
    newname = sprintf('%s_mindwandering.set', MW.subject);
    pop_saveset(MW,'filepath',outputDir,'filename',newname);

    % TRANCE 1
    idx = find(contains({EEG.event.type}, {'trance1_start' 'trance1_end'}));
    if iSub == 1
        lats = [ EEG.event(idx(1)).latency+EEG.srate*10 EEG.event(idx(2)).latency-EEG.srate*300 ];
    else
        lats = [ EEG.event(idx(1)).latency+EEG.srate*10 EEG.event(idx(2)).latency-EEG.srate*5 ];
    end
    TR1 = pop_select(EEG, 'point', lats);
    fprintf('TR 1 length = %g min \n', TR1.xmax/60);

    % TRANCE 2
    idx = find(contains({EEG.event.type}, {'trance2_start' 'trance2_end'}));
    if iSub == 1
        lats = [ EEG.event(idx(1)).latency+EEG.srate*10 EEG.event(idx(2)).latency-EEG.srate*180 ];
    else
        lats = [ EEG.event(idx(1)).latency+EEG.srate*10 EEG.event(idx(2)).latency-EEG.srate*5 ];
    end
    TR2 = pop_select(EEG, 'point', lats);
    fprintf('TR 2 length = %g min \n', TR2.xmax/60);
    
    % TRANCE 3
    idx = find(contains({EEG.event.type}, {'trance3_start' 'trance3_end'}));
    lats = [ EEG.event(idx(1)).latency+EEG.srate*10 EEG.event(idx(2)).latency-EEG.srate*5 ];
    TR3 = pop_select(EEG, 'point', lats);
    fprintf('TR 3 length = %g min \n', TR3.xmax/60);
    
    % Merge TRANCE trials and save 
    TR = pop_mergeset(TR1, TR2);
    TR = pop_mergeset(TR, TR3);
    TR.subject = sprintf('sub-%2.2d',iSub);
    TR.condition = 'trance';
    TR.saved = 'no';
    TR = eeg_checkset(TR);
    newname = sprintf('%s_trance.set', TR.subject);
    pop_saveset(TR,'filepath',outputDir,'filename',newname);

    fprintf('  =====>>>> SUBJECT %g: MINDWANDERING length = %g min \n', iSub ,MW.xmax/60);
    fprintf('  =====>>>> SUBJECT %g: TRANCE length = %g min \n', iSub ,TR.xmax/60);

end

disp('Done'); gong
