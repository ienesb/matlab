function r = test_fim_rcp()
% RCP limits (Phase 1): M = 1 reduces RCP-OFDM to CP-OFDM(Tcp = 0); the RCP-vs-CP
% departure is O(tau/T); tau = 0 recovers the CP mean but not the CP FIM.
    r = [];
    N = 16; df = 15e3;  Ts = 1/(N*df);  sigma2 = 1;
    cr.P = 2;  cr.alpha = [0.9*exp(1j*0.3); 0.6*exp(-1j*1.2)];  cr.tau = [1.3; 2.6]*Ts;
    g1 = make_grid(N, 1, df, 0, 'cp');  g1r = make_grid(N, 1, df, 0, 'rcp');
    cr.nu = [0.3; -0.45]/(8*g1.TS);  cr.eta_bar = eta_pack(cr.alpha, cr.tau, cr.nu);
    X1 = gen_symbols(N, 1, '4qam', crn_stream(31, 0, 1));
    for ici = [false true]
        r = tc(r, sprintf('rcp-ofdm(M=1) == cp-ofdm(Tcp=0) (ici=%d)', ici), ...
            relerr(fim(X1, 'rcp-ofdm', g1r, cr, ici, sigma2), fim(X1, 'cp-ofdm', g1, cr, ici, sigma2)), 1e-12);
    end
    M = 8;
    g0 = make_grid(N, M, df, 0, 'cp');  gr = make_grid(N, M, df, 0, 'rcp');
    X = gen_symbols(N, M, '4qam', crn_stream(31, 0, 2));
    Jc = fim(X, 'cp-ofdm', g0, cr, true, sigma2);  Jr = fim(X, 'rcp-ofdm', gr, cr, true, sigma2);
    fprintf('    (info) RCP-OFDM vs CP-OFDM(Tcp=0) FIM departure: %.3e (tau/T up to %.3f)\n', relerr(Jr, Jc), max(cr.tau)/g0.T);
    r = tc(r, 'RCP vs CP departure bounded (< 0.5)', relerr(Jr, Jc), 0.5);
    cz = cr;  cz.tau(:) = 0;  cz.eta_bar = eta_pack(cz.alpha, cz.tau, cz.nu);
    r = tc(r, 'tau = 0: RCP mean == CP mean', relerr(io_rcp_ofdm(X, gr, cz, true), io_cp_ofdm(X, g0, cz, true)), 1e-13);
    d = relerr(fim(X, 'rcp-ofdm', gr, cz, true, sigma2), fim(X, 'cp-ofdm', g0, cz, true, sigma2));
    fprintf('    (info) tau = 0: RCP vs CP FIM differ by %.3e (boundary term, expected nonzero)\n', d);
end
