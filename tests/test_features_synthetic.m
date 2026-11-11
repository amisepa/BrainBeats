function test_features_synthetic
% TEST_FEATURES_SYNTHETIC - HRV and EEG features on synthetic data with known
% values, and the same through BRAINBEATS_PROCESS ('features').
%
% 1. GET_HRV_FEATURES on NN series (70 bpm) with LF (0.10 Hz, 30 ms) and HF
%    (0.25 Hz, 40 ms) sinusoids: band powers (A^2/2) with each PSD method,
%    LF/HF, normalized power, SDNN and RMSSD (analytic), Poincare identities,
%    PRSA signs, short recordings (90 s: HF only; 50 s: NaN), and pNN50 at
%    exactly 50 ms (> 50 ms, Task Force 1996).
% 2. GET_EEG_FEATURES on 19 channels of white noise plus band-center sines
%    (2.5, 6, 10, 20, 35 Hz), alpha x1.5 on the right: absolute band powers
%    (A^2/2 + noise), relative powers, units of eeg_norm 0/1/2, peak and
%    individual alpha frequency, alpha asymmetry (ln 2.25), its normalized
%    version, individualized bounds, variance; median and 90% edge
%    frequencies on white noise.
% 3. A synthetic EEG + ECG dataset (ECG built from the beat times of part 1)
%    through BRAINBEATS_PROCESS: interval timing and HRV of the detected
%    beats vs the true NN, beat count, heart rate, EEG features of the 19 EEG
%    channels.
%
% Parts 2 and 3 need EEGLAB (readlocs, dipfit's BESA template); ~1 min.
% Values on 2026-11-08: HRV band powers within 2% (300 s), SDNN and RMSSD
% within 0.5%; EEG band powers within 3.2% (120 s), asymmetry 0.803-0.819
% (expected 0.808); median and 90% edge frequencies within 0.03 Hz (they were
% half a bin, 0.24 Hz, low before 2026-11-08); pipeline HRV within 0.4% of
% the true-NN HRV here and within 1% at 500 and 1000 Hz. Until 2026-11-09 the
% R-peak times were on the 4-ms sample grid (intervals off by 1.6 ms RMS,
% pNN50 21.7 vs 19.7%); GET_RR now refines them to a fraction of a sample.
%
% Copyright (C) - Cedric Cannard, 2026

bbroot = fileparts(fileparts(mfilename('fullpath')));
p = strsplit(path, pathsep);
old = p(contains(p,'BrainBeats','IgnoreCase',true) & ~startsWith(p,bbroot));
if ~isempty(old), rmpath(old{:}); end
addpath(bbroot, fullfile(bbroot,'functions'));

%% 1. HRV features on synthetic NN series
m = 60/70; aL = 0.030; aH = 0.040; fL = 0.10; fH = 0.25; sn = 0.002;
nyq = 1/(2*m);                                  % Nyquist of the beat sampling
eL = (aL^2/2 + sn^2*0.11/nyq)*1e6;              % expected LF and HF (ms^2), white noise included
eH = (aH^2/2 + sn^2*0.25/nyq)*1e6;
eSDNN = sqrt(aL^2/2 + aH^2/2 + sn^2)*1e3;
eRMSSD = sqrt(2*aL^2*sin(pi*fL*m)^2 + 2*aH^2*sin(pi*fH*m)^2 + 2*sn^2)*1e3;
rng(1);
[NN, t] = nn_series(300, 0, m, aL, aH, fL, fH, sn);
for spec = {'LombScargle' 'welch' 'fft' 'LombScargle_norm'}
    H = hrv(NN, t, struct('hrv_spec', spec{1}));
    F = H.frequency;
    assert(abs(F.lfhf - eL/eH) < 0.05, '%s: LF/HF %.3f, expected %.3f', spec{1}, F.lfhf, eL/eH)
    if ~strcmp(spec{1}, 'LombScargle_norm')
        assert(abs(F.lf/eL - 1) < 0.05 && abs(F.hf/eH - 1) < 0.05, '%s: LF %.1f HF %.1f ms^2, expected %.1f and %.1f', spec{1}, F.lf, F.hf, eL, eH)
    end
    fprintf('HRV %-16s LF %6.1f (%.3f) HF %6.1f (%.3f) LF/HF %.2f\n', spec{1}, F.lf, F.lf/eL, F.hf, F.hf/eH, F.lfhf);
end
T = H.time; P = H.nonlinear.Poincare;
assert(abs(T.SDNN/eSDNN - 1) < 0.02 && abs(T.RMSSD/eRMSSD - 1) < 0.02, 'SDNN %.1f RMSSD %.1f ms, expected %.1f and %.1f', T.SDNN, T.RMSSD, eSDNN, eRMSSD)
assert(abs(P.SD1^2/(T.RMSSD^2/2) - 1) < 0.01 && abs(P.SD2^2/(2*T.SDNN^2 - P.SD1^2) - 1) < 0.01, 'Poincare SD1 %.2f SD2 %.2f do not match SDNN and RMSSD', P.SD1, P.SD2)
ac = H.nonlinear.PRSA_AC; dc = H.nonlinear.PRSA_DC;
assert(ac < 0 && dc > 0 && abs(-ac/dc - 1) < 0.2, 'PRSA AC %.2f DC %.2f ms', ac, dc)
H = hrv(NN, t, struct('hrv_norm', true));
F = H.frequency;
assert(abs(F.lf + F.hf - 1) < 1e-3 && abs(F.lf - eL/(eL+eH)) < 0.02, 'normalized LF %.4f HF %.4f, expected LF %.4f', F.lf, F.hf, eL/(eL+eH))
fprintf('HRV time: SDNN %.1f (expected %.1f), RMSSD %.1f (%.1f) ms; normalized LF %.3f (%.3f); PRSA AC %.1f DC %.1f ms\n', ...
    T.SDNN, eSDNN, T.RMSSD, eRMSSD, F.lf, eL/(eL+eH), ac, dc);

% Short recordings: 90 s gives HF from the whole recording, no LF, and no
% normalized power; 50 s gives NaN band powers
[NN, t] = nn_series(90, 0, m, aL, aH, fL, fH, sn);
H = hrv(NN, t, struct()); F = H.frequency; hf90 = F.hf;
assert(~isfield(F,'lf') && abs(hf90/eH - 1) < 0.1, '90 s: LF present or HF %.1f ms^2 (expected %.1f)', hf90, eH)
lastwarn('');
H = hrv(NN, t, struct('hrv_norm', true), 'frequency'); F = H.frequency;   % (nonlinear features warn too)
assert(isnan(F.hf) && contains(lastwarn, 'without LF'), '90 s, hrv_norm: HF %g (expected NaN with a warning)', F.hf)
[NN, t] = nn_series(50, 0, m, aL, aH, fL, fH, sn);
lastwarn('');
H = hrv(NN, t, struct(), 'frequency'); F = H.frequency;
assert(all(isnan([F.ulf F.vlf F.lf F.hf])) && contains(lastwarn, 'too short'), '50 s: band powers not NaN')
fprintf('HRV short recordings OK: 90 s HF %.1f ms^2 (expected %.1f), normalized NaN; 50 s NaN\n', hf90, eH);

% pNN50: differences of exactly 50 ms are not counted, whatever their
% floating-point value (0.85 - 0.80 < 0.05 < 0.90 - 0.85 in double)
for pair = [0.80 0.85; 0.85 0.90; 0.90 0.95; 1.00 1.05]'
    NN = repmat(pair, 60, 1);
    H = hrv(NN, cumsum(NN), struct(), 'time');
    assert(H.time.pNN50 == 0, 'pNN50 %g%% for successive differences of exactly 50 ms (NN %g/%g s)', H.time.pNN50, pair)
end
NN = repmat([0.85; 0.901], 60, 1);
H = hrv(NN, cumsum(NN), struct(), 'time');
assert(H.time.pNN50 == 100, 'pNN50 %g%% for differences of 51 ms', H.time.pNN50)
disp('HRV pNN50 OK: 0% at exactly 50 ms, 100% at 51 ms')

%% 2. EEG features on synthetic EEG
if ~exist('readlocs','file'), eeglab; close; end
fs = 250; rng(3);
[X, chanlocs, A, isRight, iNoise] = eeg_series(fs, 120);
dens = 4/(fs/2);                                 % white noise PSD (uV^2/Hz), SD 2 uV
expAsy = log((2.25*A(3)^2/2 + dens*5)/(A(3)^2/2 + dens*5));
for norm = [0 1 2]
    [~, E] = evalc('get_eeg_features(X, struct(''fs'', fs, ''chanlocs'', chanlocs, ''eeg_time'', true, ''eeg_frequency'', true, ''eeg_nonlinear'', false, ''eeg_norm'', norm))');
    fq = E.frequency; q = E.qeeg; f = fq.freqs; df = mean(diff(f)); B = fq.bands;
    absv = [q.delta_abs q.theta_abs q.alpha_abs q.beta_abs q.gamma_abs];
    rel = [q.delta_rel q.theta_rel q.alpha_rel q.beta_rel q.gamma_rel];
    band = [fq.delta fq.theta fq.alpha fq.beta fq.gamma];
    masks = arrayfun(@(i) f >= B(i,1) & (f < B(i,2) | (i == 5 & f <= B(i,2))), 1:5, 'UniformOutput', false);
    sig = setdiff(1:numel(chanlocs), iNoise);
    a2 = repmat(A.^2/2, numel(chanlocs), 1); a2(isRight,3) = 2.25*A(3)^2/2; a2(iNoise,:) = 0;
    expAbs = a2 + dens*df*cellfun(@sum, masks);
    err = abs(absv(sig,:)./expAbs(sig,:) - 1);
    assert(max(err(:)) < 0.04, 'eeg_norm %d: absolute band power off by up to %.1f%%', norm, 100*max(err(:)))
    assert(max(abs(sum(rel,2) - 1)) < 1e-9, 'relative band powers do not sum to 1')
    switch norm
        case 0, ref = cell2mat(cellfun(@(k) mean(fq.pwr(:,k),2), masks, 'UniformOutput', false));
        case 1, ref = cell2mat(cellfun(@(k) mean(fq.pwr(:,k),2), masks, 'UniformOutput', false));   % pwr in dB
        case 2, ref = rel;
    end
    assert(max(abs(band(:) - ref(:))) < 1e-9, 'eeg_norm %d: band values are not in the documented units', norm)
end
assert(all(abs(q.paf(sig) - 10) <= df), 'peak alpha frequency %s, expected 10 Hz', mat2str(q.paf(sig)'))
assert(abs(fq.IAF_mean - 10) < 0.2, 'IAF %.2f Hz, expected 10', fq.IAF_mean)
assert(numel(fq.asymmetry) == 8 && max(abs(fq.asymmetry - expAsy)) < 0.03, 'asymmetry %s, expected %.3f on 8 pairs', mat2str(fq.asymmetry'), expAsy)
v = var(X, 0, 2); ev = sum(a2,2) + 4;
assert(max(abs(E.time.var./ev - 1)) < 0.03 && max(abs(E.time.var - v)) < 1e-9, 'variance off')
[~, E2] = evalc('get_eeg_features(X, struct(''fs'', fs, ''chanlocs'', chanlocs, ''eeg_time'', false, ''eeg_frequency'', true, ''eeg_nonlinear'', false, ''asy_norm'', true))');
pr = fq.asymmetry_pairs_num; tot = q.total_abs;
assert(max(abs(E2.frequency.asymmetry - (fq.asymmetry - log(tot(pr(:,2))./tot(pr(:,1)))))) < 2e-3, 'normalized asymmetry is not ln(relative alpha right / left)')
[~, E3] = evalc('get_eeg_features(X, struct(''fs'', fs, ''chanlocs'', chanlocs, ''eeg_time'', false, ''eeg_frequency'', true, ''eeg_nonlinear'', false, ''eeg_freqbounds'', ''individualized''))');
Bi = E3.frequency.bands;
assert(Bi(3,1) > 7 && Bi(3,1) < 10 && Bi(3,2) > 10 && Bi(3,2) < 14 && Bi(2,2) == Bi(3,1) && Bi(4,1) == Bi(3,2), 'individualized bands %s', mat2str(Bi))
fprintf('EEG bands OK: absolute power within %.1f%%, paf %.2f Hz, IAF %.2f Hz, asymmetry %.3f-%.3f (expected %.3f), individualized alpha %.2f-%.2f Hz\n', ...
    100*max(err(:)), median(q.paf(sig)), fq.IAF_mean, min(fq.asymmetry), max(fq.asymmetry), expAsy, Bi(3,:));

% Median and 90% edge frequencies on white noise: flat spectrum from the
% lower edge of the first bin to the upper edge of the last
W = 2*randn(10, fs*300);   %#ok<NASGU> used in evalc
[~, E] = evalc('get_eeg_features(W, struct(''fs'', fs, ''chanlocs'', chanlocs(1:10), ''eeg_time'', false, ''eeg_frequency'', true, ''eeg_nonlinear'', false))');
f = E.frequency.freqs; df = mean(diff(f)); lo = f(1) - df/2; hi = f(end) + df/2;
mf = mean(E.qeeg.median); sef = mean(E.qeeg.sef90);
assert(abs(mf - (lo+hi)/2) < 0.1 && abs(sef - (lo + 0.9*(hi-lo))) < 0.1, 'white noise: median %.2f Hz (expected %.2f), SEF90 %.2f (%.2f)', mf, (lo+hi)/2, sef, lo + 0.9*(hi-lo))
fprintf('EEG spectral edges OK on white noise: median %.2f Hz (expected %.2f), SEF90 %.2f (%.2f)\n', mf, (lo+hi)/2, sef, lo + 0.9*(hi-lo));

%% 3. Synthetic EEG + ECG through brainbeats_process
rng(2);
[NN, t, r] = nn_series(300, 1, m, aL, aH, fL, fH, sn);
N = round((r(end) + 1)*fs);
ecg = ecg_from_beats(r, fs, N);
rng(4);
X = eeg_series(fs, N/fs);
EEG = eeg_emptyset;
EEG.data = single([X; ecg]); EEG.srate = fs; EEG.nbchan = size(EEG.data,1); EEG.pnts = N;
EEG.trials = 1; EEG.xmin = 0; EEG.chanlocs = chanlocs; EEG.chanlocs(end+1).labels = 'ECG';
EEG.setname = 'synthetic';
[~, EEG] = evalc('eeg_checkset(EEG)');   %#ok<ASGLU> used in evalc
[~, EEG] = evalc(['brainbeats_process(EEG,''analysis'',''features'',''heart_signal'',''ECG'',''heart_channels'',{''ECG''},' ...
    '''clean_eeg'',false,''eeg_features'',{''time'' ''frequency''},''vis_cleaning'',false,''vis_outputs'',false,''save'',false)']);
Hp = EEG.brainbeats.features.HRV; Ep = EEG.brainbeats.features.EEG;
Ht = hrv(NN, t, struct());
nn = EEG.brainbeats.preprocessings.NN;
assert(numel(nn) == numel(NN), '%d NN intervals detected, %d true', numel(nn), numel(NN))
e = 1000*sqrt(mean((nn(:) - NN(:)).^2));       % interval error, ms RMS (sample grid: ~1.6 at 250 Hz)
assert(e < 0.5 && abs(Hp.time.pNN50 - Ht.time.pNN50) < 1, 'intervals off by %.2f ms RMS, pNN50 %.1f vs %.1f%%', e, Hp.time.pNN50, Ht.time.pNN50)
d = [Hp.time.SDNN/Ht.time.SDNN Hp.time.RMSSD/Ht.time.RMSSD Hp.frequency.lf/Ht.frequency.lf Hp.frequency.hf/Ht.frequency.hf ...
    Hp.nonlinear.PRSA_AC/Ht.nonlinear.PRSA_AC Hp.nonlinear.PRSA_DC/Ht.nonlinear.PRSA_DC] - 1;
assert(max(abs(d)) < 0.02, 'pipeline vs true NN: SDNN, RMSSD, LF, HF, PRSA AC, DC differ by %s%%', mat2str(round(100*d,1)))
assert(abs(Hp.time.heart_rate - 60/mean(NN)) < 0.3, 'heart rate %.1f bpm, true %.1f', Hp.time.heart_rate, 60/mean(NN))
assert(size(Ep.qeeg.alpha_abs,1) == numel(chanlocs), 'EEG features on %d channels, expected %d (ECG excluded)', size(Ep.qeeg.alpha_abs,1), numel(chanlocs))
aR = Ep.qeeg.alpha_abs(isRight); aLft = Ep.qeeg.alpha_abs(~isRight & ~ismember(1:numel(chanlocs), iNoise)');
assert(abs(median(aR)/median(aLft) - 2.25) < 0.1, 'pipeline alpha power right/left %.2f, expected 2.25', median(aR)/median(aLft))
fprintf('Pipeline OK: %d beats, intervals within %.2f ms RMS; SDNN, RMSSD, LF, HF, PRSA within %.1f%% of the true NN; pNN50 %.1f (true %.1f); HR %.1f bpm; EEG features on %d channels\n', ...
    numel(nn)+1, e, 100*max(abs(d)), Hp.time.pNN50, Ht.time.pNN50, Hp.time.heart_rate, size(Ep.qeeg.alpha_abs,1));
end

function H = hrv(NN, t, opt, domain)   %#ok<INUSL> NN, t used in evalc
% get_hrv_features with the given options, every domain unless one is named
if nargin < 4
    p = struct('hrv_time', true, 'hrv_frequency', true, 'hrv_nonlinear', true);
else
    p = struct('hrv_time', strcmp(domain,'time'), 'hrv_frequency', strcmp(domain,'frequency'), 'hrv_nonlinear', false);
end
for k = fieldnames(opt)', p.(k{1}) = opt.(k{1}); end
[~, H] = evalc('get_hrv_features(NN, t, p)');
end

function [NN, t, r] = nn_series(dur, t0, m, aL, aH, fL, fH, sn)
% NN intervals (s) modulated at fL and fH, each set by the time of the beat
% that starts it; t: time of the beat that ends it; r: all beat times
r = t0; NN = []; t = [];
while r(end) < t0 + dur
    x = m + aL*sin(2*pi*fL*r(end)) + aH*sin(2*pi*fH*r(end) + 1) + sn*randn;
    r(end+1,1) = r(end) + x; NN(end+1,1) = x; t(end+1,1) = r(end); %#ok<AGROW>
end
end

function [X, chanlocs, A, isRight, iNoise] = eeg_series(fs, dur)
% 19 channels (10-20, BESA template): white noise (SD 2 uV) plus sines at
% 2.5, 6, 10, 20 and 35 Hz, alpha x1.5 on the right; Pz noise only
labs = {'Fp1' 'Fp2' 'F7' 'F3' 'Fz' 'F4' 'F8' 'T7' 'C3' 'Cz' 'C4' 'T8' 'P7' 'P3' 'Pz' 'P4' 'P8' 'O1' 'O2'};
[~, L] = evalc('readlocs(fullfile(fileparts(which(''eeglab'')),''plugins'',''dipfit'',''standard_BESA'',''standard-10-5-cap385.elp''))');
[ok, ix] = ismember(lower(labs), lower({L.labels}));
assert(all(ok), 'channel labels missing from the BESA template')
chanlocs = L(ix);
f0 = [2.5 6 10 20 35]; A = [6 4 8 3 1.5];
isRight = [chanlocs.Y]' < -1e-3;
iNoise = find(strcmp(labs, 'Pz'));
N = round(fs*dur); tt = (0:N-1)/fs;
X = 2*randn(numel(labs), N);
for c = setdiff(1:numel(labs), iNoise)
    a = A; if isRight(c), a(3) = 1.5*A(3); end
    for b = 1:5
        X(c,:) = X(c,:) + a(b)*sin(2*pi*f0(b)*tt + 2*pi*rand);
    end
end
end

function ecg = ecg_from_beats(r, fs, N)
% ECG (uV): P, Q, R, S, T waves as Gaussians at each beat time, respiratory
% baseline wander, white noise and a DC offset
w = [0.15 -0.10 1.0 -0.25 0.30]*1000; mu = [-0.18 -0.025 0 0.025 0.25]; sd = [0.025 0.008 0.010 0.010 0.040];
tt = (0:N-1)/fs; ecg = zeros(1, N);
for k = 1:numel(r)
    i0 = max(1, round((r(k)-0.4)*fs)); i1 = min(N, round((r(k)+0.6)*fs));
    seg = tt(i0:i1) - r(k);
    for j = 1:5
        ecg(i0:i1) = ecg(i0:i1) + w(j)*exp(-(seg - mu(j)).^2/(2*sd(j)^2));
    end
end
ecg = ecg + 80*sin(2*pi*0.25*tt) + 20*randn(1, N) + 500;
end
