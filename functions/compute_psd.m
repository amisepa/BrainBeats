% COMPUTE_PSD - Power spectral density (or power spectrum) of each channel
% with Welch's method (pwelch).
%
% Usage:
%   [pwr, pwr_db, f] = compute_psd(eegData, winSize, taperM, overlap, nfft, Fs, fRange, type, useGPU)
%   [pwr, pwr_db, f] = compute_psd(EEG.data, EEG.srate*2, 'hamming', 50, [], EEG.srate, [1 40], 'psd', false)
%
% Inputs:
%   eegData - data (channels x samples)
%   winSize - window length in samples (default = Fs*2, i.e. 2 s)
%   taperM  - taper window: 'hamming' (default), 'hann', 'blackman', 'rectwin'
%   overlap - window overlap in % (default = 50)
%   nfft    - number of FFT points (default = next power of 2 of winSize)
%   Fs      - sample rate in Hz (required)
%   fRange  - frequency range to keep in Hz (default = [1/(Fs/2) Fs/2])
%   type    - 'psd' (default, uV^2/Hz) or 'power' (PSD scaled by the
%             equivalent noise bandwidth of the window: power at each frequency, uV^2)
%   useGPU  - compute on GPU (true) or not (false, default)
%
% Outputs:
%   pwr     - PSD or power (channels x frequencies)
%   pwr_db  - same in decibels, 10*log10(pwr)
%   f       - frequencies (Hz, column vector) within fRange
%
% Copyright (C) - Cedric Cannard, 2021

function [pwr, pwr_db, f] = compute_psd(eegData,winSize,taperM,overlap,nfft,Fs,fRange,type,useGPU)

% Error if no sampling rate provided
if ~exist('Fs', 'var') || isempty(Fs)
    errordlg('You need to provide the sampling rate Fs to use this function.'); return;
end

% Window size
if ~exist('winSize', 'var') || isempty(winSize) 
    disp('Window size not provided: 2-s windows');
    winSize = Fs*2;
end

% Taper
if ~exist('taperM', 'var') || isempty(taperM) 
    taperM = 'hamming'; %hamming (default); hann; blackman; rectwin
end
fh = str2func(taperM);

% Overlap
if ~exist('overlap', 'var') || isempty(overlap)
    overlap = 50;
end
overlap = winSize/(100/overlap); %get overlap in samples

% Frequency range default
if ~exist('fRange', 'var') || isempty(fRange)
    nyquist = Fs/2;
    fRange = [1/nyquist nyquist]; 
end

% Power type default
if ~exist('type', 'var') || isempty(type)
    type = 'psd';
end

% GPU default
if ~exist('useGPU', 'var') || isempty(useGPU)
    useGPU = false;
end

% nfft (next power of 2 of the window length in samples)
if ~exist('nfft', 'var') || isempty(nfft)
    samplesPerWindow = (winSize/Fs)*Fs;
    nfft = 2^nextpow2(samplesPerWindow);
end

% Power spectral density (PSD)
for iChan = 1:size(eegData,1)
    if useGPU
        signal = gpuArray(eegData(iChan,:));
    else
        signal = eegData(iChan,:);
    end
    [pwr(iChan,:), f] = pwelch(signal,fh(winSize),overlap,nfft,Fs,type);
end

% Truncate PSD to frequency range of interest
freq = f>= fRange(1) & f<=fRange(2);
f = f(freq);
pwr = pwr(:,freq);

% Convert to decibels (dB)
pwr_db = 10*log10(pwr);

