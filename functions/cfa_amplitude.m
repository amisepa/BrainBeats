% CFA_AMPLITUDE - Heart-locked EEG amplitude: RMS over channels of the
% heartbeat-locked average from -50 to 100 ms around the R-peak (the QRS of
% the cardiac field artifact), relative to its mean from -200 to -100 ms.
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
%   a     - amplitude in the units of D (NaN with fewer than 10 heartbeats)
%
% Copyright (C) - Cedric Cannard, 2026

function a = cfa_amplitude(D, beats, fs)
w = round(-0.2*fs):round(0.4*fs);
beats = beats(beats + w(1) > 0 & beats + w(end) <= size(D,2));
if numel(beats) < 10
    a = NaN; return
end
avg = zeros(size(D,1), numel(w));
for i = 1:numel(beats)
    avg = avg + double(D(:, beats(i) + w));
end
avg = avg / numel(beats);
t = w / fs;
avg = avg - mean(avg(:, t < -0.1), 2);
a = sqrt(mean(mean(avg(:, t >= -0.05 & t <= 0.1).^2)));
