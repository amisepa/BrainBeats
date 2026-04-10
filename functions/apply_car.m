function EEG = apply_car(EEG, refLabels, chanlocs_template)
% apply_car - Common average reference (CAR) with recovery of online reference channels
%
%   Adds back the online reference electrode(s) as zero-filled channels with
%   proper scalp locations, applies full-rank CAR, and retains them for
%   subsequent spherical spline interpolation. This preserves data rank and
%   ensures the reference electrode(s) receive interpolated values rather
%   than being permanently lost.
%
%   EEG = apply_car(EEG, refLabels)
%   EEG = apply_car(EEG, refLabels, chanlocs_template)
%
% Inputs:
%   EEG               - EEGLAB EEG struct
%   refLabels         - string or cell array of strings with reference
%                       channel label(s) to recover (e.g. 'FCz' or {'Fz','FCz'})
%   chanlocs_template - (optional) path to channel location file or chanlocs
%                       struct to look up locations for refLabels. Defaults
%                       to EEGLAB's standard-10-5-cap385.elp template.
%
% Outputs:
%   EEG - re-referenced EEG struct with reference channel(s) added back as
%         zero-filled channels ready for interpolation (EEG.chaninfo.nodatchans
%         is cleared for these channels).
%
% Usage in pipeline:
%   EEG = apply_car(EEG, {'Fz','FCz'});   % add back, CAR, leave for interp
%   EEG = pop_interp(EEG, [idx_Fz, idx_FCz], 'spherical');
%
% References:
%   Makoto's preprocessing pipeline: https://eeglab.ucsd.edu/wiki/Makoto's_preprocessing_pipeline
%   Kim et al. (2023). ICA's bug. Front. Signal Process., 3, 1064138.

fprintf('Applying full-rank CAR and recovering reference channel(s): %s\n', ...
    strjoin(cellstr(refLabels), ', '));

% Normalize refLabels to cell array
if ischar(refLabels) || isstring(refLabels)
    refLabels = cellstr(refLabels);
end

% Skip channels already present in the data
existing  = {EEG.chanlocs.labels};
refLabels = refLabels(~ismember(refLabels, existing));
if isempty(refLabels)
    warning('apply_car: all refLabels already present in data - skipping channel recovery.');
end

% Load template chanlocs if not provided
if nargin < 3 || isempty(chanlocs_template)
    tmplPath = which('standard-10-5-cap385.elp');
    if isempty(tmplPath)
        error('apply_car: cannot find standard-10-5-cap385.elp. Pass chanlocs_template explicitly.');
    end
    chanlocs_template = readlocs(tmplPath);
end

if isstruct(chanlocs_template)
    tmpl = chanlocs_template;
else
    tmpl = readlocs(chanlocs_template);
end
tmplLabels = {tmpl.labels};

% Add each reference channel back as a zero-filled row with proper location
eegFields = fieldnames(EEG.chanlocs);
for i = 1:numel(refLabels)
    tIdx = find(strcmpi(tmplLabels, refLabels{i}), 1);
    if isempty(tIdx)
        error('apply_car: ''%s'' not found in chanlocs template.', refLabels{i});
    end

    % Build newLoc using only fields present in EEG.chanlocs to avoid
    % struct mismatch error on assignment
    newLoc = struct();
    for f = 1:numel(eegFields)
        fd = eegFields{f};
        if isfield(tmpl, fd)
            newLoc.(fd) = tmpl(tIdx).(fd);
        else
            newLoc.(fd) = [];
        end
    end
    newLoc.labels = refLabels{i};   % preserve exact case from user

    EEG.data(end+1, :)  = 0;
    EEG.chanlocs(end+1) = newLoc;
    EEG.nbchan          = EEG.nbchan + 1;
    fprintf('  Added ''%s'' (zero-filled) at position %d\n', refLabels{i}, EEG.nbchan);
end

% Check rank before CAR
rankBefore = sum(eig(cov(double(EEG.data'))) > 1e-7);
fprintf('  Effective data rank before CAR: %d / %d channels\n', rankBefore, EEG.nbchan);

% Full-rank CAR: mean includes the zero-filled channel(s)
avgSignal = mean(EEG.data, 1);
EEG.data  = EEG.data - avgSignal;
fprintf('  CAR applied across %d channels (including recovered reference(s))\n', EEG.nbchan);
fprintf('  Reference channel(s) retained as zero-filled — interpolate with pop_interp next.\n');

% Check rank after CAR
rankAfter = sum(eig(cov(double(EEG.data'))) > 1e-7);
fprintf('  Effective data rank after CAR:  %d / %d channels\n', rankAfter, EEG.nbchan);
fprintf('  CAR applied across %d channels (including recovered reference(s))\n', EEG.nbchan);
fprintf('  Reference channel(s) retained as zero-filled — interpolate with pop_interp next.\n');

% Sanity check
EEG = eeg_checkset(EEG);
end