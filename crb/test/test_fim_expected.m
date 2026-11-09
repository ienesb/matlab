function r = test_fim_expected()
% Monte-Carlo realization FIM vs closed-form expected FIM (port of
% old/eva/verify_crb.m; 'full' tier): OFDM and CP-OTFS, {no-ICI, ICI} x {P=1, 3}
% on the legacy preset, plus the RCP-OTFS sampled single-target closed form
% and its O(1/M) Doppler-correction ordering (rcp_otfs_fim.tex).
    r = [];
    oracle_paths('add', 'old');  c = onCleanup(@() oracle_paths('remove', 'old'));
    N = 16; M = 8; df = 15e3;  sigma2 = 1.3;  nT = 20000;
    g0 = make_grid(N, M, df, 0, 'cp');
    ch = chan_eva(g0, 'legacy_preset', true);  cr = chan_realize(ch, g0, []);
    for ici = [false true]
        for P = [1 3]
            cp.P = P; cp.alpha = cr.alpha(1:P); cp.tau = cr.tau(1:P); cp.nu = cr.nu(1:P); cp.eta_bar = eta_pack(cp.alpha, cp.tau, cp.nu);
            Jo = 0; Jt = 0;
            rs = crn_stream(41, ici + 2*(P==3), 1);
            for t = 1:nT
                X = gen_symbols(N, M, '4qam', rs);  x = gen_symbols(N, M, '4qam', rs);
                Jo = Jo + fim(X, 'cp-ofdm', g0, cp, ici, sigma2);
                Jt = Jt + fim(x, 'cp-otfs', g0, cp, ici, sigma2);
            end
            Jo = Jo/nT;  Jt = Jt/nT;
            Je = fim_expected('cp-ofdm', g0, cp, ici, sigma2);
            if ici, Jold = getFIM_OFDM_multi_ICI_expected(N, M, df, g0.TS, sigma2, cp.alpha, cp.tau, cp.nu);
            else,   Jold = getFIM_OFDM_multi_noICI_expected(N, M, df, g0.TS, sigma2, cp.alpha, cp.tau, cp.nu); end
            r = tc(r, sprintf('fim_expected == old expected (ici=%d, P=%d)', ici, P), relerr(Je, Jold), 1e-12);
            r = tc(r, sprintf('OFDM MC avg == expected (ici=%d, P=%d)', ici, P), relerr(Jo, Je), 2e-2);
            r = tc(r, sprintf('OTFS MC avg == expected (ici=%d, P=%d)', ici, P), relerr(Jt, Je), 2e-2);
            if ~ici
                % no-ICI expected FIM does not depend on tau or nu
                cq = cp;  cq.tau = cq.tau + 0.37*g0.Ts;  cq.nu = cq.nu + 0.2/(M*g0.TS);
                Jq = fim_expected('cp-ofdm', g0, cq, false, sigma2);
                if P == 1
                    r = tc(r, 'no-ICI expected FIM invariant to (tau, nu), P=1', relerr(Jq, Je), 1e-12);
                end
            end
        end
    end
    % ICI expected FIM depends on nu only (not tau), single target
    cp.P = 1; cp.alpha = cr.alpha(2); cp.tau = cr.tau(2); cp.nu = cr.nu(2); cp.eta_bar = eta_pack(cp.alpha, cp.tau, cp.nu);
    cq = cp;  cq.tau = cq.tau + 0.4*g0.Ts;
    r = tc(r, 'ICI expected FIM invariant to tau, P=1', relerr(fim_expected('cp-ofdm', g0, cq, true, sigma2), fim_expected('cp-ofdm', g0, cp, true, sigma2)), 1e-12);
    % nu -> 0: the kernel itself collapses (S_a -> P_a).  The derivative kernels
    % do NOT vanish at nu = 0 (I'(q) = T/(j2pi q) for q ~= 0), which is exactly
    % why 'no-ICI' (I' := 0) and 'nu = 0' are different regimes.
    mo = ici_moments(N, 0, g0);
    r = tc(r, 'ici_moments at nu=0: S_a = P_a', ...
        max([abs(mo.S0 - N), abs(mo.S1 - power_sums(N,1)), abs(mo.S2 - power_sums(N,2))]) / power_sums(N,2), 1e-12);
    r = tc(r, 'ici_moments at nu=0: U > 0 (derivative kernels survive)', double(mo.U <= 0), 0);

    % RCP-OTFS sampled single-target closed form vs MC (integer tap), both regimes
    gr = make_grid(N, M, df, 0, 'rcp');
    cs.P = 1;  cs.alpha = cr.alpha(2);  cs.tau = 2*g0.Ts;  cs.nu = cr.nu(2);  cs.eta_bar = eta_pack(cs.alpha, cs.tau, cs.nu);
    for ici = [true false]
        Jm = 0;  rs = crn_stream(43, ici, 1);
        for t = 1:nT
            x = gen_symbols(N, M, '4qam', rs);
            Jm = Jm + fim(x, 'rcp-otfs', gr, cs, ici, sigma2);
        end
        Jm = Jm/nT;
        Je = fim_expected('rcp-otfs', gr, cs, ici, sigma2);
        r = tc(r, sprintf('RCP-OTFS sampled MC avg == closed form (ici=%d)', ici), relerr(Jm, Je), 3e-2);
    end
    Jc = fim_expected('cp-otfs', gr, cs, false, sigma2);   % CP form with TS = T (Tcp = 0)
    Jr = fim_expected('rcp-otfs', gr, cs, true, sigma2);
    ratio = Jr(4,4)/Jc(4,4);  pred = 1 + 2*g0.Ts*g0.T*power_sums(N,1)*power_sums(M,1)/(g0.T^2*N*power_sums(M,2)) ...
                                     + g0.Ts^2*M*power_sums(N,2)/(g0.T^2*N*power_sums(M,2));
    r = tc(r, 'RCP/CP J_nunu ratio = 1 + 2TsT P1Q1/(T^2 N Q2) + Ts^2 M P2/(T^2 N Q2)', abs(ratio - pred), 1e-12);
    fprintf('    (info) J_nunu RCP/CP = %.5f at M = %d (leading 1 + 3/(2M) = %.5f)\n', ratio, M, 1 + 3/(2*M));
end
