function ok = run_all_tests(tier)
% RUN_ALL_TESTS  Run the matlab/crb test suite.
%   run_all_tests            quick tier (Phase 0 + Phase 1 identities, < 5 min)
%   run_all_tests('full')    adds the Monte-Carlo closed-form tests (~4 min)
%   run_all_tests('slow')    full + the legacy end-to-end MCRB regression (~15 min)
% Prints PASS/FAIL per check; errors out (nonzero exit in -batch) on failure.
    if nargin < 1, tier = 'quick'; end
    tests = {'test_units', 'test_kernels', 'test_transforms', 'test_io', ...
             'test_channels', 'test_rng', 'test_crb_units', ...
             'test_fim_fd', 'test_fim_regression_old', 'test_fim_unitarity', 'test_fim_rcp', ...
             'test_fim_mt40_liu', 'test_rx', 'test_mcrb'};
    if any(strcmpi(tier, {'full', 'slow'}))
        tests = [tests, {'test_fim_expected', 'test_mcrb_expected'}];
    end
    if strcmpi(tier, 'slow')
        tests = [tests, {'test_regression_old'}];
    end
    all_r = [];
    t0 = tic;
    for i = 1:numel(tests)
        if exist(tests{i}, 'file') ~= 2
            fprintf('== %s: (not yet written, skipped)\n', tests{i});  continue;
        end
        fprintf('== %s\n', tests{i});
        try
            r = feval(tests{i});
        catch ME
            fprintf('  ERROR: %s\n', ME.message);
            r = tc([], [tests{i} ' (threw)'], Inf, 0);
        end
        all_r = [all_r, r]; %#ok<AGROW>
    end
    npass = sum([all_r.pass]);  ntot = numel(all_r);
    fprintf('\n%d / %d checks passed in %.1f s\n', npass, ntot, toc(t0));
    ok = npass == ntot;
    if ~ok
        fails = all_r(~[all_r.pass]);
        fprintf('FAILED:\n');  fprintf('  %s\n', fails.name);
        error('run_all_tests:fail', '%d checks failed', ntot - npass);
    end
end
