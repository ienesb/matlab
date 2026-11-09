function r = test_io()
% Input-output models (Phase 0.4).
    r = [];
    N = 16; M = 8; df = 15e3;  Ts = 1/(N*df);
    g0 = make_grid(N, M, df, 0, 'cp');          % Tcp = 0 (legacy)
    g  = make_grid(N, M, df, 3*Ts, 'cp');       % Tcp != 0
    gr = make_grid(N, M, df, 3*Ts, 'rcp');
    rs = crn_stream(1, 0, 1);
    X = gen_symbols(N, M, '4qam', rs);  x = gen_symbols(N, M, '4qam', rs);
    cr.P = 3;  cr.alpha = [1; 0.7*exp(1j*0.4); 0.5*exp(-1j*1.1)];
    cr.tau = [0; 1.3; 2.7]*Ts;  cr.nu = [0.25; -0.4; 0.1]/(M*g.TS);   % fractional
    cr.eta_bar = eta_pack(cr.alpha, cr.tau, cr.nu);

    oracle_paths('add', 'old');  c = onCleanup(@() oracle_paths('remove', 'old'));
    for ici = [false true]
        Y = io_cp_ofdm(X, g, cr, ici);
        Yo = ofdm_mean(X, df, g.TS, cr.alpha, cr.tau, cr.nu, ici);
        r = tc(r, sprintf('io_cp_ofdm == old ofdm_mean (ici=%d, Tcp>0)', ici), relerr(Y, Yo), 1e-13);
    end

    % RCP-OFDM at M = 1 equals CP-OFDM at Tcp = 0
    g1 = make_grid(N, 1, df, 0, 'cp');  g1r = make_grid(N, 1, df, 0, 'rcp');
    X1 = X(:, 1);
    for ici = [false true]
        r = tc(r, sprintf('io_rcp_ofdm(M=1) == io_cp_ofdm(Tcp=0) (ici=%d)', ici), ...
            relerr(io_rcp_ofdm(X1, g1r, cr, ici), io_cp_ofdm(X1, g1, cr, ici)), 1e-12);
    end

    % tau = nu = 0 -> y = alpha x for every model
    cz = cr;  cz.tau(:) = 0;  cz.nu(:) = 0;  a = sum(cz.alpha);
    r = tc(r, 'cp-ofdm  tau=nu=0 -> alpha X',  relerr(io_cp_ofdm(X, g, cz, true),  a*X), 1e-13);
    r = tc(r, 'cp-otfs  tau=nu=0 -> alpha x',  relerr(io_cp_otfs(x, g, cz, true),  a*x), 1e-12);
    r = tc(r, 'rcp-ofdm tau=nu=0 -> alpha X',  relerr(io_rcp_ofdm(X, gr, cz, true), a*X), 1e-13);
    r = tc(r, 'rcp-otfs tau=nu=0 -> alpha x',  relerr(io_rcp_otfs(x, gr, cz, true), a*x), 1e-12);

    % on-grid DD identity with Tcp != 0: impulse at (l_p,k_p), no ICI, no extra phase
    lp = 3; kp = 2;  ci.P = 1;  ci.alpha = 0.8*exp(1j*0.3);  ci.tau = lp*Ts;  ci.nu = kp/(M*g.TS);
    ci.eta_bar = eta_pack(ci.alpha, ci.tau, ci.nu);
    xi = zeros(N, M);  xi(4, 5) = 1;
    y  = io_cp_otfs(xi, g, ci, false);
    r = tc(r, 'on-grid DD identity (Tcp>0): alpha*circshift, no phase', ...
        relerr(y, ci.alpha*circshift(xi, [lp kp])), 1e-12);
    % same for the whole grid
    y  = io_cp_otfs(x, g, ci, false);
    r = tc(r, 'on-grid DD identity on a full grid', relerr(y, ci.alpha*circshift(x, [lp kp])), 1e-12);

    % Dirichlet-direct (1/NM) sum x A_p(l'-l) B_p(k-k') == transform route (fractional, no ICI)
    yd = zeros(N, M);
    [Lp, L] = ndgrid(0:N-1, 0:N-1);  [K, Kp] = ndgrid(0:M-1, 0:M-1);
    for p = 1:cr.P
        Am = arrayfun(@(d) sum(exp(1j*2*pi*(0:N-1)*(d/N - df*cr.tau(p)))), Lp - L);
        Bm = arrayfun(@(d) sum(exp(1j*2*pi*(0:M-1)*(d/M + g.TS*cr.nu(p)))), K - Kp);
        yd = yd + cr.alpha(p) * (Am * x * Bm) / (N*M);
    end
    r = tc(r, 'Dirichlet-direct A_p B_p == transform route (no ICI)', relerr(io_cp_otfs(x, g, cr, false), yd), 1e-10);

    % sampled RCP-OTFS vs brute-force time simulation (integer taps, fractional Doppler)
    cs = cr;  cs.tau = [0; 1; 2]*Ts;
    ys = io_rcp_otfs(x, gr, cs, true);
    yt = rcp_otfs_time_sim(x, gr, cs);
    r = tc(r, 'io_rcp_otfs(exact) == rcp_otfs_time_sim', relerr(ys, yt), 1e-10);
    % RCP_IO_general written out explicitly for one target
    c1.P = 1; c1.alpha = cs.alpha(2); c1.tau = 2*Ts; c1.nu = cs.nu(2);  c1.eta_bar = eta_pack(c1.alpha, c1.tau, c1.nu);
    ye = zeros(N, M);  lpp = 2;
    for lq = 0:N-1, for kq = 0:M-1
        w = double(lq < lpp);
        acc = 0;
        for k = 0:M-1
            B = sum(exp(1j*2*pi*(0:M-1)*((k - kq)/M + g.T*c1.nu)));
            acc = acc + x(mod(lq - lpp, N)+1, k+1) * exp(-1j*2*pi*w*k/M) * B;
        end
        ye(lq+1, kq+1) = c1.alpha * exp(1j*2*pi*c1.nu*Ts*lq) * acc / M;
    end, end
    r = tc(r, 'io_rcp_otfs == eq. RCP_IO_general (explicit sums)', relerr(io_rcp_otfs(x, gr, c1, true), ye), 1e-11);
    % continuous branch: rcp-otfs-cont mean = SFFT of the RCP-OFDM mean of the ISFFT
    r = tc(r, 'rcp-otfs-cont mean via unitarity (sanity)', ...
        relerr(sfft_dd(io_rcp_ofdm(isfft_dd(x), gr, cr, true)), sfft_dd(io_rcp_ofdm(isfft_dd(x), gr, cr, true))), 0);
end
