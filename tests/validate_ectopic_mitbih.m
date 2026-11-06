function [S, T] = validate_ectopic_mitbih(dbPath, records, variants)
% VALIDATE_ECTOPIC_MITBIH - Ectopic-beat step of CLEAN_RR on the MIT-BIH
% Arrhythmia Database (Moody & Mark, 2001), whose beats were labeled by two
% cardiologists: how many premature beats each rule replaces, and how many
% normal beats it replaces by mistake.
%
% Each record goes through CLEAN_RR once per rule ('none': steps 2-5 only;
% 'ratio': the v1.6 rule; 'adaptive': the default) with two inputs:
%   annotated - the annotated beat times (all beat types): the rules alone
%   detected  - the R-peaks GET_RR finds on the MLII lead, with their
%               amplitudes and the filtered ECG, as BRAINBEATS_PROCESS runs
%               them; matched one-to-one to the annotated beats within 150 ms
% A beat is replaced when it is not among the real beats CLEAN_RR returns
% (removed, or moved to the middle of its two intervals). Premature beats:
% V (ventricular), A, a, J, S (supraventricular). False positives: normal
% beats (N, L, R) in sinus rhythm ('(N') between two normal beats, replaced
% by the rule but not by 'none' (by the ectopic step). The first and last
% beat of each record are left out, and the paced records (102, 104, 107,
% 217), as in AAMI EC57.
%
% Usage:
%   [S, T] = validate_ectopic_mitbih(dbPath)
%   [S, T] = validate_ectopic_mitbih(dbPath, records)   % e.g. {'100','119'}
%   [S, T] = validate_ectopic_mitbih(dbPath, [], {'name', {clean_rr args}; ...})
%
% Inputs:
%   dbPath  - folder of the database (.hea, .dat, .atr and RECORDS files),
%             https://physionet.org/content/mitdb/1.0.0/
%   records - record names (default: the 44 records that are not paced)
%   variants - more CLEAN_RR settings to compare, one row each: a name and
%             its name-value pairs (default none)
%
% Outputs:
%   S - struct array, one element per record: name, t (s), label, rhythm,
%       normal (normal beat in sinus rhythm between normal beats), detected,
%       repl_ann (beats x rules, replaced from the annotated beats), kept_det
%       (beats x rules, detected and kept from the detected beats), rmssd
%       (ms: NN of the annotations, all intervals, then each rule on the
%       annotated beats), burden (share of abnormal beats within +/-45
%       beats), alternation (RR_ALTERNATION after cleaning; rows: annotated,
%       detected beats; one column per rule)
%   T - summary table: input, rule, class, beats, count (replaced, or
%       detected and kept in the NN series for 'detected')
%
% Results (2026-11-06; 44 records, 100,645 beats; 166 s). Annotation counts
% match the published ones (record 100: 2273 beats, 33 A, 1 V).
% Annotated beats, replaced by none / ratio / adaptive:
%   normal, sinus (66,838)      ectopic step 19 (0.28 per 1000; 17 in 219) /
%                               44 (0.66; 11 records, 17 in 209, 9 in 108);
%                               min shift 0: 433 (6.5 per 1000)
%   V (6,901)                   3.0 / 46.1 / 16.6%
%   A (2,544)                   3.2 / 8.5 / 15.4%  (1381 in 232)
%   isolated A, a, J, S (639)   5.5 / 21.6 / 52.1%
%   premature beats by the share of abnormal beats within +/-45 beats:
%     0-5% (505)                19.2 / 48.5 / 75.8%
%     5-15% (1,216)             8.8 / 52.7 / 55.3%
%     15-30% (2,688)            2.4 / 53.5 / 9.6%
%     30-100% (5,271)           1.1 / 21.7 / 5.3%
%   The adaptive thresholds (5.2 x the quartile deviation of |dRR| over 91
%   beats) grow with frequent ectopic beats: the rule assumes sparse ones.
%   RMSSD |error| vs the NN series of the annotations, 13 records >= 95% in
%   sinus rhythm with >= 1% premature beats: median uncleaned 53.4 ms, none
%   26.8, ratio 14.7, adaptive 6.0; record 119 (bigeminy): true 34.5,
%   adaptive 407.9, ratio 276.4.
% Detected beats (GET_RR): normal sinus 99.7% detected, V 91.3%, a 52%;
%   337 false detections (153 in 207, 86 in 108); ectopic-step false
%   positives 0.30 (ratio) and 0.74 (adaptive) per 1000 sinus beats.
% RR_ALTERNATION after 'adaptive' on the detected beats: above 10% (the
%   warning of BRAINBEATS_PROCESS) on 14 records: all 7 with atrial
%   fibrillation or flutter (201 202 203 210 219 221 222; 11-30%) and 7 of
%   the 10 others with over 15% abnormal beats (106 119 200 208 223 228 233;
%   15-44%); not on 207 (ventricular flutter, false detections; 4.5%), 213
%   (fusion beats; RMSSD error 6 ms) and 232 (atrial premature beats with
%   long pauses; 0.5%). At most 4.3% on the other 27 records (108, noisy).
%
% Reference:
%   Moody, G. B., & Mark, R. G. (2001). The impact of the MIT-BIH Arrhythmia
%   Database. IEEE Engineering in Medicine and Biology Magazine, 20(3), 45-50.
%
% Copyright (C) - Cedric Cannard, 2026

if nargin < 2 || isempty(records)
    records = strsplit(strtrim(fileread(fullfile(dbPath, 'RECORDS'))));
    records = setdiff(records, {'102' '104' '107' '217'}, 'stable');
end
bb = fileparts(fileparts(mfilename('fullpath')));
p = strsplit(path, pathsep);
old = p(contains(p,'BrainBeats','IgnoreCase',true) & ~startsWith(p,bb));
if ~isempty(old), rmpath(old{:}); end
addpath(fullfile(bb,'functions'));

rules = {'none' 'ratio' 'adaptive'};
ruleArgs = cellfun(@(r) {'ectopic_rule', r}, rules, 'UniformOutput', false);
if nargin > 2 && ~isempty(variants)
    rules = [rules variants(:,1)'];  ruleArgs = [ruleArgs variants(:,2)'];
end
nR = numel(rules);
tol = 0.15;                                  % s, beat matching (AAMI EC57)
isNorm = @(l) ismember(l, 'NLR');
rmssd = @(x) 1000*sqrt(mean(diff(x).^2));

S = struct([]);
for r = 1:numel(records)
    rec = records{r};
    [lab, smp, rhyAt, rhyLab] = read_annotations(fullfile(dbPath, [rec '.atr']));
    [ecg, fs] = read_mlii(fullfile(dbPath, rec));
    tb = smp / fs;  nb = numel(tb);
    ir = sum(rhyAt(:)' <= smp(:), 2);
    rhy = repmat({''}, nb, 1);  rhy(ir > 0) = rhyLab(ir(ir > 0));
    inner = [false; true(nb-2, 1); false];
    nrm = isNorm(lab);
    normal = inner & nrm & [false; nrm(1:end-1)] & [nrm(2:end); false] & strcmp(rhy, '(N');

    % annotated beats
    repl = false(nb, nR);  rm = zeros(1, nR);  alt = zeros(2, nR);
    for k = 1:nR
        [nn, nt, ii] = run_clean(tb(2:end), diff(tb), [], smp(2:end) + 1, ruleArgs{k});
        repl(:,k) = ~ismember(tb, nt(~ii));
        rm(k) = rmssd(nn);
        alt(1,k) = rr_alternation(nn);
    end
    iv = diff(tb);  nnOK = nrm(1:end-1) & nrm(2:end);
    dNN = diff(iv);  dNN = dNN(nnOK(1:end-1) & nnOK(2:end));
    rmTrue = 1000*sqrt(mean(dNN.^2));

    % detected beats
    tms = (0:numel(ecg)-1) / fs * 1000;
    [~, RR, RR_t, pk, sigF] = evalc('get_RR(ecg, tms, struct(''heart_signal'',''ecg'',''fs'',fs))');
    td = [RR_t(1) - RR(1); RR_t(:)];
    [dA, iA] = nearest_in(tb, td);
    [~, iD]  = nearest_in(td, tb);
    det = dA <= tol & iD(iA) == (1:nb)';
    kept = false(nb, nR);
    for k = 1:nR
        [nn, nt, ii] = run_clean(RR_t(:), RR(:), sigF(pk), pk(:), [{'ecg_signal', sigF, 'fs', fs} ruleArgs{k}]);
        kd = ismember(td, nt(~ii));
        kept(det,k) = kd(iA(det));
        alt(2,k) = rr_alternation(nn);
    end

    S(r).name = rec;  S(r).t = tb;  S(r).label = lab;  S(r).rhythm = rhy;
    S(r).normal = normal;  S(r).inner = inner;  S(r).detected = det;
    S(r).repl_ann = repl;  S(r).kept_det = kept;  S(r).rmssd = [rmTrue rmssd(iv) rm];
    S(r).burden = movmean(double(~nrm), 91);    % share of abnormal beats, 91-beat window
    S(r).alternation = alt;                     % RR_ALTERNATION after cleaning; rows: annotated, detected
    S(r).n_false_det = numel(td) - sum(det);
    fprintf('%s: %d beats (%d premature), %d detected, %d false detections; RMSSD NN %.1f, all %.1f, %s ms\n', ...
        rec, sum(inner), sum(inner & ismember(lab, 'VAaJS')), sum(det & inner), S(r).n_false_det, ...
        rmTrue, rmssd(iv), named(rules, rm));
end

%% Summary
lab = vertcat(S.label);  rhy = vertcat(S.rhythm);  inner = vertcat(S.inner);
normal = vertcat(S.normal);  det = vertcat(S.detected);
repl = vertcat(S.repl_ann);  kept = vertcat(S.kept_det);
iso = arrayfun(@(s) [false; isNorm(s.label(1:end-2)) & isNorm(s.label(3:end)); false], S, 'UniformOutput', false);
iso = vertcat(iso{:});                       % both neighbors normal
classes = {'V', 'V', 'isolated V', 'V'; 'A', 'A', 'isolated S (A a J S)', 'AaJS'; 'a', 'a', '', ''; ...
    'J', 'J', '', ''; 'S', 'S', '', ''; 'F', 'F', '', ''; 'E', 'E', '', ''; 'j', 'j', '', ''; 'e', 'e', '', ''};
T = table();
fprintf('\n%d records, %d beats (first and last of each left out)\n', numel(S), sum(inner));
fprintf('\nAnnotated beats: replaced by each rule (%s)\n', strjoin(rules, ' / '));
for c = 1:size(classes, 1)
    sel = inner & lab == classes{c,2};
    if ~any(sel), continue; end
    T = [T; summary_rows('annotated', rules, classes{c,1}, sel, repl)]; %#ok<AGROW>
    fprintf('  %-22s %6d  %s\n', classes{c,1}, sum(sel), pct(repl(sel,:)));
    if ~isempty(classes{c,3})
        sel = inner & ismember(lab, classes{c,4}) & iso;
        T = [T; summary_rows('annotated', rules, classes{c,3}, sel, repl)]; %#ok<AGROW>
        fprintf('  %-22s %6d  %s\n', classes{c,3}, sum(sel), pct(repl(sel,:)));
    end
end
sel = normal;
T = [T; summary_rows('annotated', rules, 'normal, sinus', sel, repl)];
fprintf('  %-22s %6d  %s\n', 'normal, sinus', sum(sel), pct(repl(sel,:)));
for k = 2:nR
    fp = normal & repl(:,k) & ~repl(:,1);
    rr = arrayfun(@(s) sum(s.normal & s.repl_ann(:,k) & ~s.repl_ann(:,1)), S);
    fprintf('    %s, ectopic step: %d of %d (%.2f per 1000) in %d of %d records (most: %s)\n', rules{k}, sum(fp), sum(normal), ...
        1000*sum(fp)/sum(normal), sum(rr > 0), numel(S), most(S, rr));
end
sel = inner & isNorm(lab) & ismember(rhy, {'(AFIB' '(AFL'});
T = [T; summary_rows('annotated', rules, 'normal, AF/flutter', sel, repl)];
fprintf('  %-22s %6d  %s\n', 'normal, AF/flutter', sum(sel), pct(repl(sel,:)));
burden = vertcat(S.burden);  prem = inner & ismember(lab, 'VAaJS');
edges = [0 .05 .15 .30 1.01];
fprintf('  premature beats by the share of abnormal beats within +/-45 beats:\n');
for e = 1:numel(edges)-1
    sel = prem & burden >= edges(e) & burden < edges(e+1);
    cls = sprintf('premature, %.0f-%.0f%%', 100*edges(e), min(100, 100*edges(e+1)));
    T = [T; summary_rows('annotated', rules, cls, sel, repl)]; %#ok<AGROW>
    fprintf('  %-22s %6d  %s\n', cls, sum(sel), pct(repl(sel,:)));
end
fprintf('  per record (>= 10 premature beats): premature beats, abnormal %%, replaced %%\n');
for r = 1:numel(S)
    s = S(r);  pr = s.inner & ismember(s.label, 'VAaJS');
    if sum(pr) < 10, continue; end
    fprintf('    %s %5d %5.1f%%  %s\n', s.name, sum(pr), 100*mean(~isNorm(s.label)), pct(s.repl_ann(pr,:)));
end

fprintf('\nDetected beats (GET_RR, MLII): detected; then detected and kept in the NN series (%s)\n', strjoin(rules, ' / '));
for c = 1:size(classes, 1)
    sel = inner & lab == classes{c,2};
    if ~any(sel), continue; end
    T = [T; summary_rows('detected', rules, classes{c,1}, sel, kept)]; %#ok<AGROW>
    fprintf('  %-22s %6d  %5.1f%%  %s\n', classes{c,1}, sum(sel), 100*mean(det(sel)), pct(kept(sel,:)));
    if ~isempty(classes{c,3})
        sel = inner & ismember(lab, classes{c,4}) & iso;
        T = [T; summary_rows('detected', rules, classes{c,3}, sel, kept)]; %#ok<AGROW>
        fprintf('  %-22s %6d  %5.1f%%  %s\n', classes{c,3}, sum(sel), 100*mean(det(sel)), pct(kept(sel,:)));
    end
end
sel = normal;
T = [T; summary_rows('detected', rules, 'normal, sinus', sel, kept)];
fprintf('  %-22s %6d  %5.1f%%  %s\n', 'normal, sinus', sum(sel), 100*mean(det(sel)), pct(kept(sel,:)));
for k = 2:nR
    fp = normal & det & ~kept(:,k) & kept(:,1);
    rr = arrayfun(@(s) sum(s.normal & s.detected & ~s.kept_det(:,k) & s.kept_det(:,1)), S);
    fprintf('    %s, ectopic step: %d of %d detected (%.2f per 1000) in %d records (most: %s)\n', rules{k}, sum(fp), ...
        sum(normal & det), 1000*sum(fp)/sum(normal & det), sum(rr > 0), most(S, rr));
end
fprintf('  false detections: %d\n', sum([S.n_false_det]));

% RMSSD: records in sinus rhythm (>= 95% of beats) with >= 1% premature beats
fprintf('\nRMSSD (ms), annotated beats; records >= 95%% in sinus rhythm with >= 1%% premature beats\n');
fprintf('  %-6s %6s %8s %8s %8s %8s %8s\n', 'record', 'prem%', 'NN', 'all', rules{:});
E = [];
for r = 1:numel(S)
    s = S(r);  pm = mean(ismember(s.label(s.inner), 'VAaJS'));
    if mean(strcmp(s.rhythm(s.inner), '(N')) < .95 || pm < .01, continue; end
    fprintf('  %-6s %6.1f %8.1f %8.1f %s\n', s.name, 100*pm, s.rmssd(1), s.rmssd(2), sprintf('%8.1f ', s.rmssd(3:end)));
    E = [E; s.rmssd(2:end) - s.rmssd(1)]; %#ok<AGROW>
end
if ~isempty(E)
    fprintf('  median |error| vs NN: all %.1f, %s ms (%d records)\n', median(abs(E(:,1))), ...
        named(rules, median(abs(E(:,2:end)), 1)), size(E,1));
end
% short-long alternation left after cleaning (BRAINBEATS_PROCESS warns above 10%)
k = find(strcmp(rules, 'adaptive'), 1);
fprintf(['\nRR_ALTERNATION after ''adaptive'' (%% of beats; annotated, detected), abnormal beats %%, ' ...
    'rhythms other than sinus, RMSSD error (annotated) ms\n']);
warned = false(numel(S), 1);
for r = 1:numel(S)
    s = S(r);  warned(r) = s.alternation(2,k) > 0.10;
    oth = setdiff(unique(s.rhythm(s.inner)), {'(N' ''});
    fprintf('  %s %5.1f %5.1f   %5.1f   %-24s %+7.1f%s\n', s.name, 100*s.alternation(:,k), 100*mean(~isNorm(s.label)), ...
        strjoin(oth', ' '), s.rmssd(2+k) - s.rmssd(1), repmat('  warned', 1, warned(r)));
end
fprintf('  warned (> 10%%, detected beats): %d of %d records (%s)\n', sum(warned), numel(S), strjoin({S(warned).name}, ' '));
fprintf('validate_ectopic_mitbih done\n');
end

%% Helpers
function [nn, nt, ii] = run_clean(rr_t, rr, amp, pk, args)
[~, nn, nt, ~, ~, ii] = evalc('clean_rr(rr_t, rr, amp, pk, args{:})');
end

function [d, i] = nearest_in(x, y)
% for each x, the index of the nearest y (y sorted) and the distance
i = interp1(y, 1:numel(y), x, 'nearest', 'extrap');
d = abs(x - y(i));
end

function s = named(names, v)
c = [names(:)'; num2cell(v(:)')];
s = sprintf('%s %.1f ', c{:});
end

function s = pct(m)
s = sprintf('%6.1f%% ', 100*mean(m, 1));
end

function R = summary_rows(input, rules, cls, sel, m)
R = table(repmat({input}, numel(rules), 1), rules(:), repmat({cls}, numel(rules), 1), ...
    repmat(sum(sel), numel(rules), 1), sum(m(sel,:), 1)', ...
    'VariableNames', {'input','rule','class','beats','count'});
end

function s = most(S, n)
[n, o] = sort(n, 'descend');  o = o(n > 0);  n = n(n > 0);
s = strjoin(arrayfun(@(k) sprintf('%s %d', S(o(k)).name, n(k)), 1:min(5, numel(o)), 'UniformOutput', false), ', ');
end

function [lab, smp, rhyAt, rhyLab] = read_annotations(f)
% MIT annotation format: 16-bit little-endian words, the code in the top 6
% bits and the time from the previous annotation in the bottom 10; SKIP
% (59) is followed by a 32-bit interval (high word first), AUX (63) by its
% text, which for rhythm annotations ('+', 28) is the rhythm, e.g. '(N'.
sym = blanks(41);
sym([1:14 16 18:41]) = 'NLRaVFJASEj/Q~|sT*D"=pB^t+u?![]en@xf()r';
isBeat = false(1, 41);  isBeat([1:13 25 30 34 35 38 41]) = true;
fid = fopen(f, 'r', 'ieee-le');
w = fread(fid, inf, 'uint16=>double');
fclose(fid);
n = numel(w);  code = zeros(n, 1);  at = zeros(n, 1);  aux = cell(n, 1);
t = 0;  i = 1;  m = 0;
while i <= n
    a = floor(w(i) / 1024);  I = mod(w(i), 1024);
    if a == 0 && I == 0, break; end
    switch a
        case 59                                   % SKIP
            t = t + w(i+1)*65536 + w(i+2);  i = i + 3;
        case 63                                   % AUX, attached to the last annotation
            nw = ceil(I / 2);
            b = typecast(uint16(w(i+1:i+nw)), 'uint8');   % a row for a single word
            if m > 0, s = char(b(1:I)); aux{m} = deblank(s(:)'); end
            i = i + 1 + nw;
        case {60, 61, 62}                         % NUM, SUB, CHN
            i = i + 1;
        otherwise
            t = t + I;  m = m + 1;  code(m) = a;  at(m) = t;  i = i + 1;
    end
end
code = code(1:m);  at = at(1:m);  aux = aux(1:m);
b = code >= 1 & code <= 41;  b(b) = isBeat(code(b));
lab = sym(code(b))';  smp = at(b);
r = code == 28 & ~cellfun(@isempty, aux);
rhyAt = at(r);  rhyLab = strrep(aux(r), char(0), '');
end

function [x, fs] = read_mlii(base)
% Format 212 (two 12-bit samples in three bytes), the MLII lead in mV
h = strsplit(strtrim(fileread([base '.hea'])), newline);
h0 = sscanf(h{1}, '%*s %d %f %d');
nsig = h0(1);  fs = h0(2);
lines = strtrim(h(2:nsig+1));
c = find(endsWith(lines, 'MLII'), 1);
if isempty(c), c = 1; end
f = strsplit(lines{c});                      % file format gain bits zero ...
gain = str2double(regexp(f{3}, '^[\d.]+', 'match', 'once'));  zero = str2double(f{5});
fid = fopen([base '.dat'], 'r');
B = fread(fid, [3 inf], 'uint8=>double');
fclose(fid);
s = [B(1,:) + 256*mod(B(2,:), 16); B(3,:) + 256*floor(B(2,:) / 16)];
s(s > 2047) = s(s > 2047) - 4096;
x = (s(c,:) - zero) / gain;
end
