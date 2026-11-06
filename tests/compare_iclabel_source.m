function compare_iclabel_source(files)
% COMPARE_ICLABEL_SOURCE - For HEP, ICA is fitted on a broadband continuous
% copy of the data (1-Hz high-pass, no low-pass, line noise notched, lax
% ASR) and ICLabel classifies the components on that copy; the weights then
% clean the 0.5-30 Hz data. Would ICLabel on the 0.5-30 Hz data (same
% weights) flag other components, and does it matter? (Removing components
% before or after moving the weights is the same linear projection: only the
% data ICLabel sees differ.)
%
% Runs the HEP cleaning (Picard, default options) on each recording with a
% copy of BrainBeats (in tempdir) whose clean_eeg saves the epochs before the
% component removal; then classifies the same components on those epochs and
% compares the components flagged with the HEP thresholds (muscle .99, eye
% .9, line noise .99, channel noise .99): the ROI HEP (-300 to 600 ms), the
% frontal-pole RMS (eye) and the 20-30 Hz power (muscle) after each removal.
%
% Usage:
%   compare_iclabel_source({'sample_data/sub-047_rest.set'})
%
% Results (2026-11-05; sub-047, dataset.set, ds003838 sub-062, 040, 033,
% 070; Picard): the same components removed on 5 of 6 (1-2 each). On
% sub-033 the copy also removed IC 27 (eye .94 on the copy, .74 on the data):
% frontal-pole RMS 10.69 vs 11.35 uV, ROI HEP within 0.12 uV. On the copy,
% more ICs have muscle as their largest class (2-11 vs 0-5), but none
% reached the .99 muscle threshold: no extra removal. Classifying on the
% copy (as BrainBeats does) is kept.
%
% Copyright (C) - Cedric Cannard, 2026

if ~exist('pop_loadset','file'), eeglab; close; end
bb = fileparts(fileparts(mfilename('fullpath')));
p = strsplit(path, pathsep);
old = p(contains(p,'BrainBeats','IgnoreCase',true));
if ~isempty(old), rmpath(old{:}); end
% copy of BrainBeats whose clean_eeg saves the epochs before the removal
cp = fullfile(tempdir, 'bb_iclabel_copy');
if exist(cp,'dir'), rmdir(cp,'s'); end
mkdir(fullfile(cp,'functions'));
copyfile(fullfile(bb,'*.m'), cp); copyfile(fullfile(bb,'functions','*'), fullfile(cp,'functions'));
src = fileread(fullfile(cp,'functions','clean_eeg.m'));
key = 'badComp = find(EEG.reject.gcompreject);';
assert(numel(strfind(src, key)) == 1, 'clean_eeg: removal point not found')
src = strrep(src, key, [key newline 'if ~isempty(getenv(''BB_DUMP_ICA'')), EEGpre = EEG; save(getenv(''BB_DUMP_ICA''), ''EEGpre'', ''-v7.3''); clear EEGpre; end']);
fid = fopen(fullfile(cp,'functions','clean_eeg.m'), 'w'); fwrite(fid, src); fclose(fid);
addpath(cp, fullfile(cp,'functions'));
thr = [NaN .99 .9 NaN .99 .99 NaN];                        % brain muscle eye heart line channel other
cls = {'Brain' 'Muscle' 'Eye' 'Heart' 'Line' 'Channel' 'Other'};
roiLab = {'F1' 'F2' 'F3' 'F4' 'Fz' 'FC1' 'FC2' 'FC3' 'FC4' 'FCz' 'C1' 'C2' 'C3' 'C4' 'Cz'};
files = cellstr(files);
for k = 1:numel(files)
    [fp, fn, ext] = fileparts(files{k});
    EEG = pop_loadset('filename',[fn ext],'filepath',fp);
    EEG = pop_select(EEG,'nochannel',intersect({EEG.chanlocs.labels}, {'PPG'}));
    dump = fullfile(tempdir, sprintf('icadump_%s.mat', fn));
    setenv('BB_DUMP_ICA', dump);
    brainbeats_process(EEG,'analysis','hep','heart_signal','ECG','heart_channels',{'ECG'},'clean_eeg',true, ...
        'icamethod',1,'hep_surrogates',0,'hep_cfa','off','vis_cleaning',0,'vis_outputs',0,'save',0,'gong',0);
    setenv('BB_DUMP_ICA', '');
    S = load(dump); E = S.EEGpre;
    Pc = E.etc.ic_classification.ICLabel.classifications;            % on the copy
    D = pop_iclabel(E, 'default'); Pd = D.etc.ic_classification.ICLabel.classifications;   % on the data
    A = find(any(Pc >= thr, 2)); B = find(any(Pd >= thr, 2));
    W = E.icaweights * E.icasphere; ch = E.icachansind;
    X = double(E.data(ch,:,:)); [nC, nT, nE] = size(X); X = reshape(X, nC, []);
    act = W * X;
    clean = @(c) reshape(X - E.icawinv(:,c) * act(c,:), nC, nT, nE);
    XA = clean(A); XB = clean(B);
    lab = {E.chanlocs(ch).labels};
    win = E.times >= -300 & E.times < 600;
    roi = ismember(lab, roiLab); fpole = ismember(lab, {'Fp1' 'Fp2' 'Fpz' 'AF7' 'AF8'});
    hep = @(Y) mean(mean(Y(roi,win,:),1),3);
    dH = hep(XA) - hep(XB);
    eyeRms = @(Y) sqrt(mean(Y(fpole,win,:).^2, 'all'));
    f = (0:sum(win)-1) * E.srate / sum(win); band = f >= 20 & f < 30;
    hf = @(Y) mean(abs(fft(Y(:,win,:), [], 2)).^2 .* band, 'all') / mean(band);
    fprintf(['ICLSRC %s: %d ICs; flagged on the copy %d (%s), on the data %d (%s); ROI HEP differs by max %.3f uV ' ...
        '(RMS %.1f%% of the HEP); frontal-pole RMS %.2f vs %.2f uV; 20-30 Hz power ratio data/copy %.2f\n'], ...
        fn, size(W,1), numel(A), mat2str(A(:)'), numel(B), mat2str(B(:)'), max(abs(dH)), ...
        100*sqrt(mean(dH.^2))/sqrt(mean(hep(XA).^2)), eyeRms(XA), eyeRms(XB), hf(XB)/hf(XA));
    for c = setxor(A(:)', B(:)')
        [pc, ic] = max(Pc(c,:)); [pd, id] = max(Pd(c,:));
        fprintf('ICLSRC   IC %d: copy %s %.2f (muscle %.2f, eye %.2f, brain %.2f) | data %s %.2f (muscle %.2f, eye %.2f, brain %.2f)\n', ...
            c, cls{ic}, pc, Pc(c,2), Pc(c,3), Pc(c,1), cls{id}, pd, Pd(c,2), Pd(c,3), Pd(c,1));
    end
    nb = [sum(Pc(:,1) == max(Pc,[],2)) sum(Pd(:,1) == max(Pd,[],2))];
    nm = [sum(Pc(:,2) == max(Pc,[],2)) sum(Pd(:,2) == max(Pd,[],2))];
    fprintf('ICLSRC   largest class brain: copy %d, data %d; muscle: copy %d, data %d\n', nb, nm);
    delete(dump);
end
rmpath(cp, fullfile(cp,'functions'));
end
