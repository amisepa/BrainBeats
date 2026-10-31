% CFA_AMPLITUDE - Heart-locked EEG amplitude: RMS over channels of the
% heartbeat-locked average, relative to its mean from -200 to -100 ms, in the
% QRS window (-50 to 100 ms) and in the late window (200 to 500 ms: the T wave,
% where the amplitude also holds any neural HEP, so it is an upper bound of
% the residual artifact).
%
% Usage:
%   a = cfa_amplitude(D, beats, fs)
%
% Inputs:
%   D     - continuous data (channels x samples)
%   beats - R-peak sample indices
%   fs    - sampling rate (Hz)
%
% Output:
%   a     - [QRS late] amplitudes in the units of D (NaN with fewer than 10
%           heartbeats)
%
% Copyright (C) - Cedric Cannard, 2026

function a = cfa_amplitude(D, beats, fs)
w = round(-0.2*fs):round(0.5*fs);
beats = beats(beats + w(1) > 0 & beats + w(end) <= size(D,2));
if numel(beats) < 10
    a = [NaN NaN]; return
end
avg = zeros(size(D,1), numel(w));
for i = 1:numel(beats)
    avg = avg + double(D(:, beats(i) + w));
end
avg = avg / numel(beats);
t = w / fs;
avg = avg - mean(avg(:, t < -0.1), 2);
a = [sqrt(mean(mean(avg(:, t >= -0.05 & t <= 0.1).^2))) sqrt(mean(mean(avg(:, t >= 0.2 & t <= 0.5).^2)))];
