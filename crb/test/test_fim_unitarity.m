function r = test_fim_unitarity()
% OTFS FIMs by unitarity vs independent direct-DD oracles (Phase 1).
    r = [];
    addpath(fullfile(fileparts(mfilename('fullpath')), 'oracles'));
    N = 16; M = 8; df = 15e3;  Ts = 1/(N*df);  sigma2 = 0.8;
    g  = make_grid(N, M, df, 3*Ts, 'cp');
    gr = make_grid(N, M, df, 3*Ts, 'rcp');
    x = gen_symbols(N, M, '4qam', crn_stream(21, 0, 1));
    cr.P = 2;  cr.alpha = [0.9*exp(1j*0.3); 0.6*exp(-1j*1.2)];
    cr.tau = [1.3; 2.6]*Ts;  cr.nu = [0.3; -0.45]/(M*g.TS);  cr.eta_bar = eta_pack(cr.alpha, cr.tau, cr.nu);
    for ici = [false true]
        Ju = fim(x, 'cp-otfs', g, cr, ici, sigma2);
        r = tc(r, sprintf('cp-otfs == cp-ofdm(ISFFT x) (ici=%d)', ici), relerr(Ju, fim(isfft_dd(x), 'cp-ofdm', g, cr, ici, sigma2)), 1e-13);
        Jd = fim_cp_otfs_direct(x, g, cr, ici, sigma2);
        r = tc(r, sprintf('cp-otfs == direct DD oracle (A/B or K kernels) (ici=%d)', ici), relerr(Ju, Jd), 1e-9);
    end
    % continuous RCP-OTFS branch
    Jc = fim(x, 'rcp-otfs-cont', gr, cr, true, sigma2);
    r = tc(r, 'rcp-otfs-cont == rcp-ofdm(ISFFT x)', relerr(Jc, fim(isfft_dd(x), 'rcp-ofdm', gr, cr, true, sigma2)), 1e-13);
    cs = cr;  cs.tau = [1; 3]*Ts;  cs.eta_bar = eta_pack(cs.alpha, cs.tau, cs.nu);
    Js = fim(x, 'rcp-otfs', gr, cs, true, sigma2);
    Jc = fim(x, 'rcp-otfs-cont', gr, cs, true, sigma2);
    fprintf('    (info) sampled vs continuous RCP-OTFS FIM, integer taps: rel gap %.3e (different models)\n', relerr(Js, Jc));
    r = tc(r, 'sampled and continuous RCP-OTFS FIMs both finite', double(~all(isfinite([Js(:); Jc(:)]))), 0);
end
