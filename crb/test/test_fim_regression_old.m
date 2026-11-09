function r = test_fim_regression_old()
% New FIM builders vs the verified old/eva/utils builders on the legacy preset
% (power-sorted pairing, N=16, M=8, fc=2 GHz, Tcp=0) -- Phase 1 regression.
    r = [];
    oracle_paths('add', 'old');  c = onCleanup(@() oracle_paths('remove', 'old'));
    N = 16; M = 8; df = 15e3;  Ts = 1/(N*df);  sigma2 = 1.3;
    g0 = make_grid(N, M, df, 0, 'cp');
    gl = make_grid(N, M, df, 0, 'rcp');
    ch = chan_eva(g0, 'legacy_preset', true);
    cr = chan_realize(ch, g0, []);
    rs = crn_stream(5, 0, 1);
    X = gen_symbols(N, M, '4qam', rs);  x = gen_symbols(N, M, '4qam', rs);

    for ici = [false true]
        if ici, fo = @getFIM_OFDM_multi_ICI; ft = @getFIM_CP_OTFS_multi_ICI;
        else,   fo = @getFIM_OFDM_multi_noICI; ft = @getFIM_CP_OTFS_multi_noICI; end
        Jo = fo(X, df, g0.TS, sigma2, cr.alpha, cr.tau, cr.nu);
        r = tc(r, sprintf('cp-ofdm  vs old (ici=%d, P=3)', ici), relerr(fim(X, 'cp-ofdm', g0, cr, ici, sigma2), Jo), 1e-12);
        Jt = ft(x, df, g0.TS, sigma2, cr.alpha, cr.tau, cr.nu);
        r = tc(r, sprintf('cp-otfs  vs old (ici=%d, P=3)', ici), relerr(fim(x, 'cp-otfs', g0, cr, ici, sigma2), Jt), 1e-12);
        % single target through the P=1 path
        c1.P = 1; c1.alpha = cr.alpha(2); c1.tau = cr.tau(2); c1.nu = cr.nu(2); c1.eta_bar = eta_pack(c1.alpha, c1.tau, c1.nu);
        r = tc(r, sprintf('cp-ofdm  vs old (ici=%d, P=1)', ici), ...
            relerr(fim(X, 'cp-ofdm', g0, c1, ici, sigma2), fo(X, df, g0.TS, sigma2, c1.alpha, c1.tau, c1.nu)), 1e-12);
        % RCP-OFDM
        Jr = getFIM_RCP_OFDM_multi(X, df, g0.TS, sigma2, cr.alpha, cr.tau, cr.nu, ici);
        r = tc(r, sprintf('rcp-ofdm vs old (ici=%d)', ici), relerr(fim(X, 'rcp-ofdm', gl, cr, ici, sigma2), Jr), 1e-12);
    end
    % CP with Tcp > 0 vs old Ts = T + Tcp
    g3 = make_grid(N, M, df, 3*Ts, 'cp');
    Jo = getFIM_OFDM_multi_ICI(X, df, g3.TS, sigma2, cr.alpha, cr.tau, cr.nu);
    r = tc(r, 'cp-ofdm vs old with Tcp > 0 (old Ts = T + Tcp)', relerr(fim(X, 'cp-ofdm', g3, cr, true, sigma2), Jo), 1e-12);

    % sampled RCP-OTFS: new use_ici = true (exact sampled model) == old use_ici = false
    cs = cr;  cs.tau = [0; 1; 2]*Ts;  cs.eta_bar = eta_pack(cs.alpha, cs.tau, cs.nu);
    Jold = getFIM_RCP_OTFS_multi(x, df, g0.TS, sigma2, cs.alpha, cs.tau, cs.nu, false);
    r = tc(r, 'rcp-otfs (exact sampled) vs old getFIM_RCP_OTFS_multi(use_ici=false)', ...
        relerr(fim(x, 'rcp-otfs', gl, cs, true, sigma2), Jold), 1e-12);
    % the old use_ici=true hybrid is NOT reproduced; report the gap for the record
    Jh = getFIM_RCP_OTFS_multi(x, df, g0.TS, sigma2, cs.alpha, cs.tau, cs.nu, true);
    fprintf('    (info) old hybrid rcp-otfs(use_ici=true) vs new exact sampled: rel %.2e; vs new cont: %.2e\n', ...
        relerr(Jh, fim(x, 'rcp-otfs', gl, cs, true, sigma2)), relerr(Jh, fim(x, 'rcp-otfs-cont', gl, cs, true, sigma2)));
end
