% Run heartbeat evoked potentials (HEP)
%
% Cedric Cannard, 2023

function run_HEP(EEG, params, Rpeaks)

% Add heartbeat markers in the signals, taking into account
% existing events
nEv = length(EEG.event);
urevents = num2cell(nEv+1:nEv+length(Rpeaks));
evt = num2cell(Rpeaks);
types = repmat({'R-peak'},1,length(evt));
[EEG.event(1,nEv+1:nEv+length(Rpeaks)).latency] = evt{:};       % assign latencies
[EEG.event(1,nEv+1:nEv+length(Rpeaks)).type] = types{:};        % assign types
[EEG.event(1,nEv+1:nEv+length(Rpeaks)).urevent] = urevents{:};  % assign event index
EEG = eeg_checkset(EEG);

% Calculate minimum distance between R peaks (in s)
minEpoch = min(diff(Rpeaks)/EEG.srate);

% Epoch
EEG = pop_epoch(EEG,{},[-minEpoch minEpoch],'epochinfo','yes');

if params.clean_eeg
    [EEG, params] = clean_eeg(EEG,params);
end

% Save
newname = sprintf('%s_HEP.set', EEG.filename(1:end-4));
pop_saveset(EEG,'filename',newname,'filepath',EEG.filepath);

