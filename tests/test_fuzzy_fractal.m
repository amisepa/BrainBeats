function test_fuzzy_fractal
% TEST_FUZZY_FRACTAL - BB_FUZZY_FRACTAL (fuzzy entropy and fractal dimension
% from the ASCENT plugin): white noise is more irregular than a noisy sine,
% the values match ASCENT's own functions, the change-of-version warning is
% given, and without ASCENT the values are NaN with a warning.
% Needs EEGLAB and ASCENT 2.1 or later on the path.
%
% Old (BrainBeats 1.6) vs ASCENT values on 66 recordings: see the v1.7 entry
% of the README version history.
%
% Copyright (C) - Cedric Cannard, 2026

bb = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(bb,'functions'));
assert(exist('compute_FuzzEn','file') > 0, 'ASCENT is not on the path')
rng(1);
x = [randn(1,2000); sin(2*pi*(1:2000)/1000) + 0.01*randn(1,2000)]; %#ok<NASGU> used in evalc: noise, slow sine (2 cycles)
lastwarn('');
[~, fe, fd, vers] = evalc('bb_fuzzy_fractal(x, 2, .15, 2, 1, false)');
[~, id] = lastwarn;
assert(strcmp(id, 'BrainBeats:fuzzyEntropy'), 'no change-of-version warning (last: %s)', id)
assert(isequal(size(fe), [2 1]) && isequal(size(fd), [2 1]), 'outputs are not column vectors')
assert(fe(1) > fe(2) && fd(1) > fd(2), 'noise not more irregular than the sine: FE %s, FD %s', mat2str(fe), mat2str(fd))
[~, ref] = evalc('compute_FuzzEn(x, ''Parallel'', false)');
assert(max(abs(fe - round(ref,3))) == 0, 'FE differs from compute_FuzzEn')
[~, ref] = evalc('compute_FracDim(x, ''Parallel'', false)');
assert(max(abs(fd - round(ref,3))) == 0, 'FD differs from compute_FracDim')
fprintf('bb_fuzzy_fractal OK (ASCENT %s): FE noise %.3f, sine %.3f; FD noise %.3f, sine %.3f\n', vers, fe, fd);

% Without ASCENT (every copy on the path): NaN and a warning
removed = {};
for f = {'compute_FuzzEn' 'eegplugin_ascent'}
    while exist(f{1}, 'file')
        removed{end+1} = fileparts(which(f{1})); %#ok<AGROW>
        rmpath(removed{end});
    end
end
c = onCleanup(@() addpath(removed{:}));
[~, fe, fd] = evalc('bb_fuzzy_fractal(x, 2, .15, 2, 1, false)');
[~, id] = lastwarn;
assert(all(isnan([fe; fd])) && strcmp(id, 'BrainBeats:ascent'), 'without ASCENT: %s, warning %s', mat2str([fe; fd]), id)
fprintf('bb_fuzzy_fractal OK without ASCENT: NaN and warning\n');
end
