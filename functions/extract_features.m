% EXTRACT_FEATURES - Gather HRV and EEG features from the brainbeats_process
% output into tables (one variable per feature).
% Standalone: not called by the main BrainBeats pipeline.
%
% Usage:
%   [hrv, eeg] = extract_features(Features)
%
% Inputs:
%   Features - structure output by brainbeats_process, with fields HRV
%              and/or EEG
%
% Outputs:
%   hrv - table of HRV features (one row): time domain, mean ULF/VLF/LF/HF/LF-HF
%         power across windows, Poincare SD1/SD2, PRSA, fuzzy entropy, fractal
%         dimension (and MFE peak scale and area under the curve if present).
%         Empty if Features.HRV does not exist.
%   eeg - table of EEG features (one row, channel features as one column per
%         channel in channel order): time domain, band power, IAF, qEEG,
%         fuzzy entropy and fractal dimension, and alpha asymmetry (one
%         variable per electrode pair). Empty if Features.EEG does not exist.
%
% Copyright (C) - Cedric Cannard, 2023

function [hrv, eeg] = extract_features(Features)

hrv = table.empty;
eeg = table.empty;

%% HRV

if isfield(Features,'HRV')


    % Time
    if isfield(Features.HRV, 'time')
        hrv = struct2table(Features.HRV.time);
    end
    
    % Frequency (average across windows)
    if isfield(Features.HRV, 'frequency')

        % ULF
        if isfield(Features.HRV.frequency,'ulf')
            hrv(:,size(hrv,2)+1) = table(mean(Features.HRV.frequency.ulf));
            hrv.Properties.VariableNames(end) = {'ULF-HRV'};
        end

        % VLF
        if isfield(Features.HRV.frequency,'vlf')
            hrv(:,size(hrv,2)+1) = table(mean(Features.HRV.frequency.vlf));
            hrv.Properties.VariableNames(end) = {'VLF-HRV'};
        end

        % LF
        if isfield(Features.HRV.frequency,'lf')
            hrv(:,size(hrv,2)+1) = table(mean(Features.HRV.frequency.lf));
            hrv.Properties.VariableNames(end) = {'LF-HRV'};
        end

        % HF
        if isfield(Features.HRV.frequency,'hf')
            hrv(:,size(hrv,2)+1) = table(mean(Features.HRV.frequency.hf));
            hrv.Properties.VariableNames(end) = {'HF-HRV'};
        end

        % LF/HF
        if isfield(Features.HRV.frequency,'lfhf')
            hrv(:,size(hrv,2)+1) = table(mean(Features.HRV.frequency.lfhf));
            hrv.Properties.VariableNames(end) = {'LF/HF'};
        end
    end

    % Nonlinear
    if isfield(Features.HRV, 'nonlinear')

        % list of fields
        features = fieldnames(Features.HRV.nonlinear);

        % Poincare
        if sum(contains(features,'Poincare')) > 0
            hrv(:,size(hrv,2)+1) = table(Features.HRV.nonlinear.Poincare.SD1);
            hrv.Properties.VariableNames(end) = {'Poincaré: SD1'};
    
            hrv(:,size(hrv,2)+1) = table(Features.HRV.nonlinear.Poincare.SD2);
            hrv.Properties.VariableNames(end) = {'Poincaré: SD2'};
    
            hrv(:,size(hrv,2)+1) = table(Features.HRV.nonlinear.Poincare.SD1SD2);
            hrv.Properties.VariableNames(end) = {'SD1/SD2'};
        end
    
        % PRSA
        if any(strcmp(features,'PRSA_AC')) && any(strcmp(features,'PRSA_DC'))
            hrv(:,size(hrv,2)+1) = table(Features.HRV.nonlinear.PRSA_AC);
            hrv.Properties.VariableNames(end) = {'PRSA_AC'};
    
            hrv(:,size(hrv,2)+1) = table(Features.HRV.nonlinear.PRSA_DC);
            hrv.Properties.VariableNames(end) = {'PRSA_DC'};
        end
    
        % Fuzzy entropy
        if any(strcmp(features,'FE'))
            hrv(:,size(hrv,2)+1) = table(Features.HRV.nonlinear.FE);
            hrv.Properties.VariableNames(end) = {'HRV-FE'};
        end

        % Fractal dimension
        if any(strcmp(features,'FD'))
            hrv(:,size(hrv,2)+1) = table(Features.HRV.nonlinear.FD);
            hrv.Properties.VariableNames(end) = {'HRV-FD'};
        end

        % Multiscale fuzzy entropy (if present)
        if any(strcmp(features,'MFE'))
            % Peak scale number
            [~, peakScale] = max(Features.HRV.nonlinear.MFE);
            hrv(:,size(hrv,2)+1) = table(peakScale);
            hrv.Properties.VariableNames(end) = {'HRV-MFE_peak'};
            
            % Area under the curve (using trapezoidal numerical integration)
            hrv(:,size(hrv,2)+1) = table(trapz(Features.HRV.nonlinear.MFE));
            hrv.Properties.VariableNames(end) = {'HRV-MFE_auc'};
        end
    
    end
end


%% EEG

if isfield(Features,'EEG')

    % Time (one value per channel)
    if isfield(Features.EEG, 'time')
        features = fieldnames(Features.EEG.time);
        for iFeat = 1:length(features)
            eeg = add_var(eeg, Features.EEG.time.(features{iFeat}), sprintf('EEG-%s', features{iFeat}));
        end
    end

    % Frequency
    if isfield(Features.EEG, 'frequency')
        freq = Features.EEG.frequency;

        % Band power (one value per channel)
        bandNames = {'delta' 'theta' 'alpha' 'beta' 'gamma'};
        for iBand = 1:length(bandNames)
            if isfield(freq, bandNames{iBand})
                eeg = add_var(eeg, freq.(bandNames{iBand}), sprintf('EEG-%s', bandNames{iBand}));
            end
        end

        % IAF (per channel and mean)
        if isfield(freq,'IAF')
            eeg = add_var(eeg, freq.IAF, 'EEG-IAF');
        end
        if isfield(freq,'IAF_mean')
            eeg = add_var(eeg, freq.IAF_mean, 'EEG-IAF_mean');
        end

        % Alpha asymmetry (one variable per electrode pair)
        if isfield(freq,'asymmetry') && isfield(freq,'asymmetry_pairs_labels')
            asy = freq.asymmetry;
            pairs = freq.asymmetry_pairs_labels;
            for iPair = 1:length(asy)
                eeg = add_var(eeg, asy(iPair), sprintf('Asy (%s)',pairs{iPair}));
            end
        end
    end

    % qEEG (one value per channel)
    if isfield(Features.EEG, 'qeeg')
        features = fieldnames(Features.EEG.qeeg);
        for iFeat = 1:length(features)
            eeg = add_var(eeg, Features.EEG.qeeg.(features{iFeat}), sprintf('qEEG-%s', features{iFeat}));
        end
    end

    % Nonlinear (one value per channel)
    if isfield(Features.EEG, 'nonlinear')
        features = fieldnames(Features.EEG.nonlinear);
        for iFeat = 1:length(features)
            eeg = add_var(eeg, Features.EEG.nonlinear.(features{iFeat}), sprintf('EEG-%s', features{iFeat}));
        end
    end
end

%% Subfunction

function tbl = add_var(tbl, x, name)
% Append x to the one-row table tbl as variable 'name' (vectors become
% one column per element).
t = table(double(gather(x(:)')), 'VariableNames', {name});
if isempty(tbl)
    tbl = t;
else
    tbl = [tbl t];
end


