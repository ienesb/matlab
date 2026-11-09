function ok = run_manuscript_checks()
% Scope-restricted validation: every dependency stays in matlab/crb.
% The complete legacy suite is deliberately not invoked because several of
% its tests add matlab/old or matlab/impl, outside this manuscript task.
    names = {'test_units','test_channels','test_rng','test_crb_units', ...
        'test_fim_fd','test_fim_unitarity','test_fim_rcp','test_fim_mt40_liu', ...
        'test_mcrb_scoped','test_mcrb_rx_stats'};
    all_r = []; t0=tic;
    for k=1:numel(names)
        fprintf('== %s\n',names{k});
        r=feval(names{k}); all_r=[all_r r]; %#ok<AGROW>
    end
    ok=all([all_r.pass]);
    fprintf('\nMANUSCRIPT CHECKS: %d/%d passed in %.2f s\n',nnz([all_r.pass]),numel(all_r),toc(t0));
    assert(ok,'run_manuscript_checks:fail','Manuscript numerical checks failed');
end
