function sim = sim_brainheart(cfg)
% SIM_BRAINHEART - Simulated EEG, ECG and respiration with a known
% brain-heart coupling (ground truth for the coupling analysis).
%
% Slow processes on a 2-Hz grid drive the signals:
%   - respiration r: ~0.25 Hz, slowly varying rate and depth;
%   - alpha log-amplitude a: AR(1), plus respiration (c_ra) and the heart
%     period one second earlier (heart -> brain, c_ha);
%   - heart rate modulation m: alpha log-amplitude two seconds earlier
%     (brain -> heart, c_ah), respiration (respiratory sinus arrhythmia,
%     c_rh), a 0.1-Hz (Mayer) oscillation and noise;
%   - theta and beta log-amplitudes: independent AR(1) (no coupling).
% Heartbeats come from an integral pulse frequency modulation (IPFM)
% model of m (Mateo & Laguna, 2003), step by step, so that heart -> brain
% uses the heart periods already produced. EEG channels mix band-limited
% oscillations (amplitudes from the slow processes) with 1/f background;
% the ECG is a sum of Gaussian waves at the beats, its QRS amplitude
% modulated by respiration; the cardiac field artifact (CFA) adds the ECG
% to the EEG with a spatial gain pattern.
%
% Usage:
%   sim = sim_brainheart(cfg)
%
% cfg fields (defaults): dur (300 s), fs (250 Hz), seed (0), c_ah (0),
%   c_ha (0), c_rh (0.6), c_ra (0.3), cfa_uv (10: R-peak amplitude of the
%   CFA in the EEG, uV, largest channel), noise_uv (10: 1/f background,
%   uV RMS), qrs_resp (0.15: QRS amplitude modulation by respiration)
%
% Output sim: data (nChan x N, uV), fs, chanlocs (labels, theta, radius),
%   ecg (1 x N, mV), resp (1 x N), beats (sample indices), NN, NN_t (s),
%   slow (struct of the 2-Hz latent series), cfg
%
% Copyright (C) - Cedric Cannard, BrainBeats 2026

if nargin < 1, cfg = struct; end
d = struct('dur',300,'fs',250,'seed',0,'c_ah',0,'c_ha',0,'c_rh',0.6,'c_ra',0.3, ...
    'cfa_uv',10,'noise_uv',10,'qrs_resp',0.15);
fn = fieldnames(d);
for k = 1:numel(fn), if ~isfield(cfg,fn{k}), cfg.(fn{k}) = d.(fn{k}); end, end
rng(cfg.seed);
fs = cfg.fs; N = round(cfg.dur*fs); t = (0:N-1)/fs;
fsl = 2; K = floor(cfg.dur*fsl); step = fs/fsl;

% Respiration (fine rate): rate 0.25 Hz +/- slow drift, depth modulated
fr = 0.25 + 0.02*smooth_noise(N, fs, 0.02);
resp = sin(2*pi*cumsum(fr)/fs) .* (1 + 0.2*smooth_noise(N, fs, 0.05));
rsl = resp(1:step:end); rsl = rsl(1:K);

% Mayer wave (0.1 Hz) AR(2) at 2 Hz
w = 2*pi*0.1/fsl; lf = filter(1, [1 -2*0.95*cos(w) 0.95^2], randn(1,K)); lf = lf/std(lf);

% Slow loop with the IPFM heart
T0 = 0.85;                                  % mean heart period (s)
a = zeros(1,K); m = zeros(1,K); h = zeros(1,K);
th = filter(1,[1 -0.8],randn(1,K)); th = th/std(th);
be = filter(1,[1 -0.8],randn(1,K)); be = be/std(be);
beats = []; phase = rand; lastBeat = -T0; NNlast = T0;
for k = 1:K
    km = @(L) max(k-L, 1);
    a(k) = 0.8*a(km(1)) + cfg.c_ra*rsl(km(1)) + cfg.c_ha*h(km(2)) + 0.6*randn;
    m(k) = cfg.c_ah*a(km(4)) + cfg.c_rh*rsl(km(1)) + 0.5*lf(k) + 0.4*randn;
    hr = (1 + 0.06*m(k))/T0;                % beats per second during this step
    for j = 1:step                          % IPFM integration at the fine rate
        phase = phase + hr/fs;
        if phase >= 1
            phase = phase - 1;
            tb = ((k-1)*step + j - 1)/fs;
            beats(end+1) = (k-1)*step + j; %#ok<AGROW>
            NNlast = tb - lastBeat; lastBeat = tb;
        end
    end
    h(k) = (NNlast - T0)/0.05;              % standardized heart period, for heart -> brain
end
a = a/std(a);

% ECG: Gaussian waves at the beats; QRS amplitude follows respiration
w8 = [-0.2 0.1 0.025; -0.03 -0.15 0.01; 0 1 0.012; 0.03 -0.25 0.01; 0.25 0.3 0.05];  % [lat amp width]
ecg = zeros(1,N);
L = round(0.45*fs); tt = (-L:L)/fs;
for b = beats
    idx = b + (-L:L); ok = idx >= 1 & idx <= N;
    s = 1 + cfg.qrs_resp*resp(b);
    wave = zeros(size(tt));
    for i = 1:size(w8,1)
        amp = w8(i,2); if i >= 2 && i <= 4, amp = amp*s; end
        wave = wave + amp*exp(-(tt - w8(i,1)).^2/(2*w8(i,3)^2));
    end
    ecg(idx(ok)) = ecg(idx(ok)) + wave(ok);
end

% Channels: locations, band source weights and CFA gains
labels = {'Fz' 'F3' 'F4' 'Cz' 'C3' 'C4' 'T7' 'T8' 'Pz' 'P3' 'P4' 'O1' 'O2'};
theta  = [0 -39 39 0 -90 90 -90 90 180 -141 141 -162 162];
radius = [.25 .33 .33 0 .25 .25 .5 .5 .25 .33 .33 .51 .51];
wAlpha = [.3 .3 .3 .5 .5 .5 .4 .4 1 1 1 1 1];
wTheta = [1 .8 .8 .7 .5 .5 .3 .3 .4 .3 .3 .2 .2];
wBeta  = [.5 .5 .5 1 1 1 .6 .6 .5 .5 .5 .4 .4];
gCFA   = [.3 .1 .5 .2 -.2 .6 -.6 1 .1 -.3 .4 -.4 .5];        % left-right gradient, like the QRS field
nCh = numel(labels);
chanlocs = struct('labels', labels, 'theta', num2cell(theta), 'radius', num2cell(radius));

% Band oscillations: band-limited noise times the slow amplitude
up = @(x) interp1((0:K-1)/fsl, x, t, 'pchip', 'extrap');
osc = @(band) band_noise(N, fs, band);
Aal = exp(0.4*up(a)); Ath = exp(0.4*up(th)); Abe = exp(0.4*up(be));
sAl = 12*Aal.*osc([8 13]); sTh = 8*Ath.*osc([4 8]); sBe = 5*Abe.*osc([13 30]);
data = wAlpha'*sAl + wTheta'*sTh + wBeta'*sBe;
for c = 1:nCh
    data(c,:) = data(c,:) + cfg.noise_uv*pink(N);
end
data = data + (cfg.cfa_uv*gCFA')*ecg;      % cardiac field artifact

beats = beats(:)';
NN = diff(beats)/fs; NN_t = (beats(2:end)-1)/fs;
sim = struct('data', data, 'fs', fs, 'chanlocs', chanlocs, 'ecg', ecg, 'resp', resp, ...
    'beats', beats, 'NN', NN, 'NN_t', NN_t, ...
    'slow', struct('t', (0:K-1)/fsl, 'alpha', a, 'm', m, 'h', h, 'resp', rsl, 'theta', th, 'beta', be), 'cfg', cfg);


function x = smooth_noise(N, fs, fc)
% Unit-variance noise low-passed at fc (Hz)
n = 2^nextpow2(N); f = (0:n-1)*fs/n; f(f > fs/2) = fs - f(f > fs/2);
X = fft(randn(1,n)) .* (f <= fc);
x = real(ifft(X)); x = x(1:N); x = x/std(x);


function x = band_noise(N, fs, band)
n = 2^nextpow2(N); f = (0:n-1)*fs/n; f(f > fs/2) = fs - f(f > fs/2);
X = fft(randn(1,n)) .* (f >= band(1) & f <= band(2));
x = real(ifft(X)); x = x(1:N); x = x/std(x);


function x = pink(N)
n = 2^nextpow2(N); f = (0:n-1)/n; f(f > 0.5) = 1 - f(f > 0.5); f(1) = f(2);
X = fft(randn(1,n)) ./ sqrt(f);
x = real(ifft(X)); x = x(1:N); x = x/std(x);
