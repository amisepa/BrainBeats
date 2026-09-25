function plot_features(Features,params)
% PLOT_FEATURES - Plot HRV and EEG features computed by BrainBeats.
%
% Usage:
%   plot_features(Features)
%   plot_features(Features, params)
%
% Inputs:
%   Features - EEG.brainbeats.features, with HRV and/or EEG substructures
%   params   - (optional) BrainBeats parameters. Missing hrv_*/eeg_* flags are
%              inferred from the fields present in Features. Fields read:
%              hrv_features, hrv_frequency, hrv_spec (HRV PSD method, for
%              the ylabel: 'LombScargle_norm' (default) = normalized,
%              otherwise s^2/Hz), eeg_features, eeg_frequency,
%              eeg_nonlinear, eeg_norm (0 = uV^2/Hz, 1 = dB (default),
%              2 = normalized; ylabel only), and chanlocs (EEG channel
%              locations, required for EEG features)
%
% Plots:
%   - PSD of HRV (areas of the ULF/VLF/LF/HF bands estimated; skipped with a
%     warning if the recording was too short for any) and of EEG (20% trimmed mean across
%     channels; delta, theta, alpha, beta, gamma areas, with the band limits
%     used for the features), side by side if both.
%   - EEG scalp topographies: mean delta/theta/alpha/beta/gamma power, IAF,
%     fuzzy entropy and fractal dimension, and a 3D headplot of alpha asymmetry.
%
% Copyright (C) - Cedric Cannard, 2023

disp('Plotting features...')

% Fill in any missing params from the features themselves (useful when users
% have computed features and want to replot them without redefining all the
% params, e.g. plot_features(Features) or with params.chanlocs only)
if nargin < 2, params = struct(); end
hasHRV = isfield(Features,'HRV') && isstruct(Features.HRV);
hasEEG = isfield(Features,'EEG') && isstruct(Features.EEG);
if ~isfield(params,'hrv_features'),  params.hrv_features  = hasHRV; end
if ~isfield(params,'hrv_frequency'), params.hrv_frequency = hasHRV && isfield(Features.HRV,'frequency'); end
if ~isfield(params,'hrv_nonlinear'), params.hrv_nonlinear = hasHRV && isfield(Features.HRV,'nonlinear'); end
if ~isfield(params,'hrv_norm'),      params.hrv_norm      = false; end   % default in get_hrv_features (units only)
if ~isfield(params,'eeg_features'),  params.eeg_features  = hasEEG; end
if ~isfield(params,'eeg_frequency'), params.eeg_frequency = hasEEG && isfield(Features.EEG,'frequency'); end
if ~isfield(params,'eeg_nonlinear'), params.eeg_nonlinear = hasEEG && isfield(Features.EEG,'nonlinear'); end
if ~isfield(params,'eeg_norm'),      params.eeg_norm      = 1; end       % default in get_eeg_features (units only)

% EEG channel locations
if params.eeg_features && ~isfield(params,'chanlocs')
    errordlg('Sorry, you need to load your EEG channel locations into params.chanlocs to plot EEG features (see tutorial).')
    return
end

% Pull features data
if params.hrv_features
    HRV = Features.HRV;
end
if params.eeg_features
    EEG = Features.EEG;
end

% abort if empty
if ~params.hrv_features && ~params.eeg_features
    fprintf('No features to plot. \n');
    return
end

%% Power spectral density (PSD)

% HRV spectrum available? (empty, and band powers NaN, when the recording is
% too short for any HRV band)
plotHRVpsd = params.hrv_features && params.hrv_frequency && isfield(HRV,'frequency') ...
    && isfield(HRV.frequency,'pwr') && ~isempty(HRV.frequency.pwr) && any(~isnan(HRV.frequency.pwr(:)));
if params.hrv_features && params.hrv_frequency && ~plotHRVpsd
    warning('No HRV power spectrum to plot (the recording is likely too short for HRV frequency features).')
end

% PSD - HRV
if plotHRVpsd

    figure('color','w');
    try icadefs; set(gcf, 'color', BACKCOLOR); catch; end     % eeglab background color

    % Subplot mode if EEG PSD was also extracted
    if params.eeg_features && params.eeg_frequency
        nexttile([2 3])
    end

    hold on
    pwr = HRV.frequency.pwr; % spectrum of the last window of each band
    freqs = HRV.frequency.pwr_freqs;
    bandNames = {'ULF' 'VLF' 'LF' 'HF'};
    bandColors = {"#A2142F" "#D95319" "#EDB120" "#0072BD"};
    bands = HRV.frequency.bands;
    baseval = min(pwr);

    % Area of each band that was estimated (field present and not NaN)
    plotted = false(1,length(bandNames));
    for iBand = 1:length(bandNames)
        fld = lower(bandNames{iBand});
        if isfield(HRV.frequency, fld) && ~all(isnan(HRV.frequency.(fld)))
            x = freqs >= bands(iBand,1) & freqs <= bands(iBand,2);
            if any(x)
                area(freqs(x),pwr(x),'BaseValue',baseval,'FaceColor',bandColors{iBand},'FaceAlpha',.7);
                plotted(iBand) = true;
            end
        end
    end
    warning off
    legend(bandNames(plotted))
    warning on
    axis tight; box on
    xlabel('Frequency (Hz)');

    % Units depend on the PSD method (params.hrv_spec, default 'LombScargle_norm')
    if isfield(params,'hrv_spec') && ~strcmp(params.hrv_spec,'LombScargle_norm')
        ylabel('Power (s^2/Hz)');
    else
        ylabel('Power (normalized)');
    end
    title(sprintf('Power spectral density - HRV'))
    set(gcf,'Toolbar','none','Menu','none');  % remove toolbar and menu
    set(gcf,'Name','Visualization of features','NumberTitle','Off')  % name
    set(findall(gcf,'type','axes'),'fontSize',12,'fontweight','bold');

end

% PSD - EEG
if params.eeg_features && params.eeg_frequency

    if plotHRVpsd
        nexttile([2 3])
    else
        figure('color','w')
        try icadefs; set(gcf, 'color', BACKCOLOR); catch; end     % eeglab background color
    end

    % PSD units for ylabel
    if isfield(params,'eeg_norm')
        if params.eeg_norm == 0
            units = 'Power (uV^2/Hz)';
        elseif params.eeg_norm == 1
            units = 'Power (db)';
        elseif params.eeg_norm == 2
            units = 'Power (normalized)';
        end
    else
        units = 'Power (db)';  % assume default
    end        

    hold on
    pwr = trimmean(EEG.frequency.pwr,20,1); % 20% trimmed mean across channels
    freqs = EEG.frequency.freqs;
    if isfield(EEG.frequency,'bands')
        bands = EEG.frequency.bands;    % band limits used by get_eeg_features
    else
        bands = [0 3; 3 7; 7 13; 13 30; 30 max(freqs)];
    end
    baseval = min(pwr);

    % delta
    x = freqs >= bands(1,1) & freqs <= bands(1,2);
    y = pwr(x);
    area(freqs(x),y,'BaseValue',baseval,'FaceColor',"#A2142F",'FaceAlpha',.7)

    % theta
    x = freqs >= bands(2,1) & freqs <= bands(2,2);
    y = pwr(x);
    area(freqs(x),y,'BaseValue',baseval,'FaceColor',"#D95319",'FaceAlpha',.7)

    % alpha
    x = freqs >= bands(3,1) & freqs <= bands(3,2);
    y = pwr(x);
    area(freqs(x),y,'BaseValue',baseval,'FaceColor',"#EDB120",'FaceAlpha',.7)

    % beta
    x = freqs >= bands(4,1) & freqs <= bands(4,2);
    y = pwr(x);
    area(freqs(x),y,'BaseValue',baseval,'FaceColor',"#0072BD",'FaceAlpha',.7)

    % gamma
    x = freqs >= bands(5,1) & freqs <= bands(5,2);
    y = pwr(x);
    area(freqs(x),y,'BaseValue',baseval,'FaceColor',"#4DBEEE",'FaceAlpha',.7)

    title('Power spectral density - EEG'); 
    warning off
    legend({'delta' 'theta' 'alpha' 'beta' 'gamma'})
    warning on
    xlabel('Frequency (Hz)');
    ylabel(units); 
    axis tight;

    set(gcf,'Toolbar','none','Menu','none');  % remove toolbar and menu
    set(gcf,'Name','Visualization of features','NumberTitle','Off')  % name
    set(findall(gcf,'type','axes'),'fontSize',12,'fontweight','bold'); % font
    pause(0.1)  % to allow plot before next plot
end

%% Scalp topos 2D

mode = 1;  % 1 for 2D, 2 for 3D (also used by the nonlinear topos below)

if params.eeg_features && params.eeg_frequency


    % Create figure
    figure('color','w','Units','Normalized','OuterPosition', [0 0 1 1],...
        'Toolbar','none','Menu','none','Name','EEG features','NumberTitle','Off')

    try icadefs; set(gcf, 'color', BACKCOLOR); catch; end     % eeglab background color

    % Mean power per channel: delta, theta, alpha, beta, gamma
    warning('off','all')
    nexttile
    plot_topo(gather(mean(EEG.frequency.delta,2)),params.chanlocs,mode,'psd');
    cb = colorbar;
    ylabel(cb,units,'Rotation',270,'fontSize',12,'fontweight','bold')

    title('Delta power');
    nexttile
    plot_topo(gather(mean(EEG.frequency.theta,2)),params.chanlocs,mode,'psd');
    cb = colorbar; 
    ylabel(cb,units,'Rotation',270,'fontSize',12,'fontweight','bold')
    title('Theta power');

    nexttile
    plot_topo(gather(mean(EEG.frequency.alpha,2)),params.chanlocs,mode,'psd');
    cb = colorbar; 
    ylabel(cb,units,'Rotation',270,'fontSize',12,'fontweight','bold')
    title('Alpha power');

    nexttile
    plot_topo(gather(mean(EEG.frequency.beta,2)),params.chanlocs,mode,'psd');
    cb = colorbar; 
    ylabel(cb,units,'Rotation',270,'fontSize',12,'fontweight','bold')
    title('Beta power');

    nexttile
    plot_topo(gather(mean(EEG.frequency.gamma,2)),params.chanlocs,mode,'psd');
    cb = colorbar; 
    ylabel(cb,units,'Rotation',270,'fontSize',12,'fontweight','bold')
    title('Gamma power');
    warning on

    % IAF
    try
        nexttile
        plot_topo(gather(EEG.frequency.IAF),params.chanlocs,mode,'psd');
        cb = colorbar; ylabel(cb,'Frequency (Hz)','Rotation',270,'fontSize',12,'fontweight','bold')
        title('Individual alpha frequency (IAF)');
    catch
        warning('Could not plot the individual alpha frequency (IAF). IAF estimation may have failed (can happen if no clear alpha peak distribution is present)')
    end

    set(findall(gcf,'type','axes'),'fontSize',12,'fontweight','bold');

end

% Nonlinear features (added to the current figure)
if params.eeg_features && params.eeg_nonlinear
    nexttile
    plot_topo(gather(EEG.nonlinear.FE),params.chanlocs,mode,'entropy');
    cb = colorbar; ylabel(cb,'Fuzzy entropy','Rotation',270,'fontSize',12,'fontweight','bold')
    title('Fuzzy entropy');

    nexttile
    plot_topo(gather(EEG.nonlinear.FD),params.chanlocs,mode,'entropy');
    cb = colorbar; ylabel(cb,'Fractal dimension','Rotation',270,'fontSize',12,'fontweight','bold')
    title('Fractal dimension');
end

if ~isempty(get(groot,'CurrentFigure'))    % (no figure if there was nothing to plot)
    set(findall(gcf,'type','axes'),'fontSize',12,'fontweight','bold');
    set(gcf,'Toolbar','none','Menu','none');  % remove toolbar and menu
    set(gcf,'Name','Visualization of features','NumberTitle','Off')  % name
end

% Alpha asymmetry 3D headplot (last, because of colorbar issues)
if params.eeg_features && params.eeg_frequency
    try
        warning('off','all')
        nexttile
        view = [-85 20];  % 'left'
        asy = EEG.frequency.asymmetry;
        pairNums = EEG.frequency.asymmetry_pairs_num;
        headplotparams = { 'meshfile','mheadnew.mat','transform',...
            [0.664455 -3.39403 -14.2521 -0.00241453 0.015519 -1.55584 11 10.1455 12],...
            'colormap',parula,'maplimits','absmax','cbar',1,...
            'electrodes','off','material','metal','verbose','off'};
        brainbeats_headplot('setup',params.chanlocs(pairNums(:,1)),...
            'tmp.spl',headplotparams{:}); % Generate temporary spline file
        brainbeats_headplot(asy,'tmp.spl','view',view,headplotparams{:});  % 3D headplot of asymmetry
        title('Alpha asymmetry')
        delete 'tmp.spl'
        set(findall(gcf,'type','axes'),'fontSize',12,'fontweight','bold');
        warning on

    catch
        warning('Failed to plot the 3D headplot of alpha asymmetry. This may happen if your EEG data are low-density (i.e., few EEG channels only)')
    end
end
