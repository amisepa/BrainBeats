function EEG = fix_events(EEG, filename)

if strcmp(filename, 'subj01_1.bdf') 
    EEG = pop_selectevent(EEG,'event',7,'renametype','rest1_start','deleteevents','off');   %rest 1 start
    EEG = pop_selectevent(EEG,'event',8,'renametype','rest1_end','deleteevents','off');   %rest 1 end
    EEG = pop_selectevent(EEG,'event',9,'renametype','trance1_start','deleteevents','off');   %trance 1 start
    EEG = pop_selectevent(EEG,'event',12,'renametype','trance1_end','deleteevents','off');   %trance 1 end
    EEG = pop_selectevent(EEG,'event',16,'renametype','trance2_start','deleteevents','off');   %trance 2 start
    EEG = pop_selectevent(EEG,'event',17,'renametype','trance2_end','deleteevents','off');   %trance 2 end
    EEG = pop_selectevent(EEG,'event',18,'renametype','rest2_start','deleteevents','off');   %rest 2 start
    EEG = pop_selectevent(EEG,'event',19,'renametype','rest2_end','deleteevents','off');   %rest 2 end
    EEG = pop_selectevent(EEG,'event',21,'renametype','rest3_start','deleteevents','off');   %rest 3 start
    EEG = pop_selectevent(EEG,'event',22,'renametype','rest3_end','deleteevents','off');   %rest 3 end
    EEG = pop_selectevent(EEG,'event',23,'renametype','trance3_start','deleteevents','off');   %trance 3 start
    EEG = pop_selectevent(EEG,'event',24,'renametype','trance3_end','deleteevents','off');   %trance 3 end

else

    % Fix error for subject 12
    if strcmp(filename, 'subj12_1.bdf')
        EEG = pop_selectevent(EEG,'event',16,'renametype','767','deleteevents','off');
    % elseif  strcmp(filename, 'subj12_2.bdf')
    %     EEG = pop_selectevent(EEG,'event',12,'renametype','4351','deleteevents','off');
    end

    % Convert event names to strings
    for iEvent = 1:length(EEG.event)
        EEG.event(iEvent).type = num2str(EEG.event(iEvent).type);
    end

    % Find start markers
    Rs = find(strcmp({EEG.event.type}, '767'));     %rest start markers
    Ts = find(strcmp({EEG.event.type}, '2303'));    %trance start markers
    if length(Rs) ~= 3 || length(Ts) ~= 3
        warning('Wrong number of start markers: there should be 3 rest_start and 3 trance_start (line 135)');
        warning('Taking first 3 events and ignoring the next ones');
        Rs = Rs(1:3);
        Ts = Ts(1:3);
    end

    % Fix event names
    [EEG.event(1:end).duration] = deal([]);

    % Start of trials
    for iEvent = 1:3
        EEG.event(Rs(iEvent)).type = ['rest' num2str(iEvent) '_start'];      %Rest start events
        EEG.event(Ts(iEvent)).type = ['trance' num2str(iEvent) '_start'];    %Trance start events
    end

    %%%%%%%%%%%%%%%% SESSION 1 %%%%%%%%%%%%%%%%%%%
    %rest end 1
    markers = {EEG.event(Rs(1)+1:Ts(1)-1).type};     %markers between rest1_start and trance1_start
    if length(markers) > 1
        target = find(strcmp(markers, '1279'));
        if length(target) > 1
            target = target(1);
        end
    else
        target = 1;
    end
    duration = ( EEG.event(Rs(1)+target).latency - EEG.event(Rs(1)).latency ) / EEG.srate / 60;
    if duration > 3
        EEG.event(Rs(1)+target).type = 'rest1_end';
        EEG.event(Rs(1)+target).duration = duration;
    else
        cprintf('red', 'cannot find rest 1 end marker'); return
    end

    % rest end 2
    markers = {EEG.event(Rs(2)+1:Rs(3)-1).type};     %markers between rest1_start and trance1_start
    if length(markers) > 1
        target = find(strcmp(markers, '4351'));
        if length(target) > 1
            target = target(1);
        end
    else
        target = 1;
    end
    duration = ( EEG.event(Rs(2)+target).latency - EEG.event(Rs(2)).latency ) / EEG.srate / 60;
    if duration > 3
        EEG.event(Rs(2)+target).type = 'rest2_end';
        EEG.event(Rs(2)+target).duration = duration;
    else
        cprintf('red', 'cannot find rest 2 end marker'); return
    end

    % rest end 3
    markers = {EEG.event(Rs(3)+1:Ts(3)-1).type};     %markers between rest1_start and trance1_start
    if length(markers) > 1
        target = find(strcmp(markers, '1279'));
        if length(target) > 1
            target = target(1);
        end
    else
        target = 1;
    end
    duration = ( EEG.event(Rs(3)+target).latency - EEG.event(Rs(3)).latency ) / EEG.srate / 60;
    if duration > 3
        EEG.event(Rs(3)+target).type = 'rest3_end';
        EEG.event(Rs(3)+target).duration = duration;
    else
        cprintf('red', 'cannot find rest 3 end marker'); return
    end

    % trance end 1
    markers = {EEG.event(Ts(1)+1:Ts(2)-1).type};     %markers between rest1_start and trance1_start
    if length(markers) > 1
        target = find(strcmp(markers, '4351'));
        if length(target) > 1
            target = target(1);
        end
    else
        target = 1;
    end
    duration = ( EEG.event(Ts(1)+target).latency - EEG.event(Ts(1)).latency ) / EEG.srate / 60;
    if duration > 3
        EEG.event(Ts(1)+target).type = 'trance1_end';
        EEG.event(Ts(1)+target).duration = duration;
    else
        cprintf('red', 'cannot find trance 1 end marker'); return
    end

    % trance end 2
    markers = {EEG.event(Ts(2)+1:Rs(2)-1).type};     %markers between rest1_start and trance1_start
    if length(markers) > 1
        target = find(strcmp(markers, '511'));
        if length(target) > 1
            target = target(1);
        end
    else
        target = 1;
    end
    duration = ( EEG.event(Ts(2)+target).latency - EEG.event(Ts(2)).latency ) / EEG.srate / 60;
    if duration > 3
        EEG.event(Ts(2)+target).type = 'trance2_end';
        EEG.event(Ts(2)+target).duration = duration;
    else
        cprintf('red', 'cannot find trance 2 end marker'); return
    end

    % trance end 3
    markers = {EEG.event(Ts(3)+1:end).type};     %markers between rest1_start and trance1_start
    if length(markers) > 1
        target = find(contains(markers, {'4351' '16639' '8447'}));
        if length(target) > 1
            target = target(1);
        end
    else
        target = 1;
    end
    duration = ( EEG.event(Ts(3)+target).latency - EEG.event(Ts(3)).latency ) / EEG.srate / 60;
    if duration > 3
        EEG.event(Ts(3)+target).type = 'trance3_end';
        EEG.event(Ts(3)+target).duration = duration;
    else
        cprintf('red', 'cannot find trance 3 end marker'); return
    end

end

EEG = pop_selectevent(EEG,'type',{'rest1_start','rest1_end','rest2_start',...
    'rest2_end','rest3_start','rest3_end','trance1_start','trance1_end',...
    'trance2_start','trance2_end','trance3_start','trance3_end'},'deleteevents','on');
EEG = eeg_checkset(EEG);

if length(EEG.event) ~= 12
    error('wrong number of trials!.');
end
