function r = test_mcrb()
% Phase-3 tests: A/B regression vs old getAB_OFDM_core, matched limits, direct
% DD oracle, finite differences of the Hessian fields, pseudo-true behaviour,
% affine-in-sigma2 assembly, MT40 runtime gate, MMLE sanity.
    r = [];
    oracle_paths('add', 'old');  c = onCleanup(@() oracle_paths('remove', 'old'));
    addpath(fullfile(fileparts(mfilename('fullpath')), 'oracles'));
    N = 16; M = 8; df = 15e3;  Ts = 1/(N*df);  sigma2 = 1.3;
    g0 = make_grid(N, M, df, 0, 'cp');
    ch = chan_eva(g0, 'legacy_preset', true);  cr = chan_realize(ch, g0, []);
    mdm = modem('4qam');
    rs = crn_stream(91, 0, 1);
    X = gen_symbols(N, M, '4qam', rs);  Xh = bsc_symbols(X, 0.1, mdm, rs);
    x = gen_symbols(N, M, '4qam', rs);  xh = bsc_symbols(x, 0.1, mdm, rs);
    % evaluation point != truth
    ce = cr;  ce.alpha = 0.8*cr.alpha;  ce.tau = cr.tau + 0.2*Ts;  ce.nu = cr.nu + 0.1/(M*g0.TS);  ce.eta_bar = eta_pack(ce.alpha, ce.tau, ce.nu);
    for ici = [false true]
        [A, B, Jt, b] = mcrb_AB_cp_ofdm(X, Xh, g0, ce, ici, sigma2, cr);
        [Ao, Bo, Jto, bo] = getAB_OFDM_core(X, Xh, df, g0.TS, sigma2, ce.alpha, ce.tau, ce.nu, ici, cr.alpha, cr.tau, cr.nu);
        r = tc(r, sprintf('A vs old getAB_OFDM_core (ici=%d, eval != true)', ici), relerr(A, Ao), 1e-12);
        r = tc(r, sprintf('B vs old getAB_OFDM_core (ici=%d)', ici),  max(relerr(B, Bo), relerr(b, bo)), 1e-12);
        [Am, Bm, Jm, bm] = mcrb_AB_cp_ofdm(X, X, g0, cr, ici, sigma2, cr);
        J = fim(X, 'cp-ofdm', g0, cr, ici, sigma2);
        r = tc(r, sprintf('matched: b = 0, A = -J, B = J (ici=%d)', ici), max([norm(bm)/norm(J), relerr(Am, -J), relerr(Bm, J), relerr(Jm, J)]), 1e-12);
        % OTFS by unitarity == OFDM at ISFFT grids
        [At, Bt] = mcrb_AB_cp_otfs(x, xh, g0, ce, ici, sigma2, cr);
        [Au, Bu] = mcrb_AB_cp_ofdm(isfft_dd(x), isfft_dd(xh), g0, ce, ici, sigma2, cr);
        r = tc(r, sprintf('cp-otfs A,B == cp-ofdm at ISFFT (ici=%d)', ici), max(relerr(At, Au), relerr(Bt, Bu)), 1e-13);
    end
    [Ad, Bd] = mcrb_AB_cp_otfs_direct(x, xh, g0, cr, sigma2);
    [At, Bt] = mcrb_AB_cp_otfs(x, xh, g0, cr, false, sigma2, cr);
    r = tc(r, 'cp-otfs A,B == direct DD oracle (no ICI)', max(relerr(At, Ad), relerr(Bt, Bd)), 1e-9);

    % finite differences of the Hessian fields (second derivatives of the mean)
    h = 1e-3;
    for ici = [false true]
        tau = cr.tau(2);  nu = cr.nu(2);  a = cr.alpha(2);
        F = fields_cp_ofdm(Xh, g0, tau, nu, ici, 2);
        mu = @(t, n) a*getfield(fields_cp_ofdm(Xh, g0, t, n, ici, 0), 'Phi');
        d2t = (mu(tau + h*Ts, nu) - 2*mu(tau, nu) + mu(tau - h*Ts, nu)) / (h*Ts)^2;
        d2n = (mu(tau, nu + h*g0.bin_nu) - 2*mu(tau, nu) + mu(tau, nu - h*g0.bin_nu)) / (h*g0.bin_nu)^2;
        dtn = (mu(tau + h*Ts, nu + h*g0.bin_nu) - mu(tau + h*Ts, nu - h*g0.bin_nu) - mu(tau - h*Ts, nu + h*g0.bin_nu) + mu(tau - h*Ts, nu - h*g0.bin_nu)) / (4*h*Ts*h*g0.bin_nu);
        e = max([relerr(d2t, -4*pi^2*a*F.Xi), relerr(d2n, -4*pi^2*a*F.Theta), relerr(dtn, 4*pi^2*a*F.Omega1)]);
        r = tc(r, sprintf('Hessian fields Xi/Theta/Omega1 vs 2nd differences (ici=%d)', ici), e, 1e-5);
    end

    % ---- pseudo-true ----
    cr2 = cr;
    [eta0, info] = pseudo_true(X, X, g0, cr2, false);          % matched -> eta_bar
    D = eta_scale(g0, cr.P);
    r = tc(r, 'pseudo-true matched -> eta_bar (bin units)', norm(D\(eta0 - cr.eta_bar))/norm(D\cr.eta_bar), 1e-9);
    tight = struct('method', 'varpro', 'ftol', 1e-10, 'steptol', 1e-9, 'opttol', 1e-8);   % old/eva tolerances
    [eta_v, iv] = pseudo_true(X, Xh, g0, cr2, true, tight);
    [eta_o, ~, ~, ~, ~] = pseudoTrue_OFDM(X, Xh, df, g0.TS, cr.alpha, cr.tau, cr.nu, true, [], true);
    r = tc(r, 'varpro == old pseudoTrue_OFDM (bin units, old tolerances)', norm(D\(eta_v - eta_o))/norm(D\eta_o), 1e-6);
    [eta_d, ~] = pseudo_true(X, Xh, g0, cr2, true);
    r = tc(r, 'default tolerances within 1e-3 bin of the tight solution', norm(D\(eta_d - eta_v)), 1e-3);
    % varpro vs full-4P: on the legacy preset the three delays sit within 0.1 bin
    % of each other (only Doppler separates them) and the two solvers can return
    % the same optimum with the targets relabelled; use a resolvable channel.
    cs = cr2;  cs.tau = [0.4; 1.7; 3.1]*Ts;  cs.eta_bar = eta_pack(cs.alpha, cs.tau, cs.nu);
    [eta_f, iff] = pseudo_true(X, Xh, g0, cs, true, struct('method', 'full4P'));
    [eta_w, iw]  = pseudo_true(X, Xh, g0, cs, true, struct('method', 'varpro'));
    r = tc(r, 'varpro == full4P (ICI, p = 0.1, resolvable taps)', norm(D\(eta_w - eta_f))/norm(D\eta_f), 1e-5);
    r = tc(r, 'varpro and full4P reach the same residual', abs(iw.resnorm - iff.resnorm)/iff.resnorm, 1e-8);
    [eta_dd, ~] = pseudo_true(x, xh, g0, cr2, false, struct('domain', 'dd'));
    [eta_tf, ~] = pseudo_true(isfft_dd(x), isfft_dd(xh), g0, cr2, false);
    r = tc(r, 'DD route == TF route', norm(D\(eta_dd - eta_tf))/norm(D\eta_tf), 1e-8);
    [eta_lm, il] = pseudo_true(X, Xh, g0, cr2, true, struct('solver', 'lm'));
    r = tc(r, 'built-in LM == lsqnonlin (same optimum, flat valley)', norm(D\(eta_lm - eta_v))/norm(D\eta_v), 1e-3);
    r = tc(r, 'built-in LM residual == lsqnonlin residual', abs(il.resnorm - iv.resnorm)/iv.resnorm, 1e-8);
    cz = cr;  cz.nu(:) = 0;  cz.eta_bar = eta_pack(cz.alpha, cz.tau, cz.nu);
    [eta_z, iz] = pseudo_true(X, Xh, g0, cz, false);
    [~, ~, nz] = eta_unpack(eta_z);
    r = tc(r, 'mask_nu (nu = 0): nu0 == nu_bar exactly', norm(nz), 0);
    r = tc(r, 'mask_nu auto-detected', double(~iz.mask_nu), 0);
    % expectation: alpha0/alpha_bar -> 1-2p, tau0/nu0 unbiased (32 x 16, 200 trials)
    Nb = 32; Mb = 16;  gb = make_grid(Nb, Mb, df, 0, 'cp');
    chb = chan_eva(gb, 'legacy_preset', true);  crb = chan_realize(chb, gb, []);  crb1 = crb;
    crb1.P = 1; crb1.alpha = crb.alpha(2); crb1.tau = 1.3*gb.Ts; crb1.nu = 0.3/(Mb*gb.TS); crb1.eta_bar = eta_pack(crb1.alpha, crb1.tau, crb1.nu);
    p = 0.1;  nT = 200;  ratio = zeros(nT, 1);  dt = zeros(nT, 1);  dn = zeros(nT, 1);
    for t = 1:nT
        rs = crn_stream(92, 0, t);
        Xb = gen_symbols(Nb, Mb, '4qam', rs);  Xbh = bsc_symbols(Xb, p, mdm, rs);
        e0 = pseudo_true(Xb, Xbh, gb, crb1, false, struct('p', p));
        [a0, t0, n0] = eta_unpack(e0);
        ratio(t) = real(a0/crb1.alpha);  dt(t) = (t0 - crb1.tau)/gb.Ts;  dn(t) = (n0 - crb1.nu)/gb.bin_nu;
    end
    r = tc(r, 'BSC p=0.1: mean alpha0/alpha_bar == 1-2p', abs(mean(ratio) - (1 - 2*p)), 2e-2);
    r = tc(r, 'BSC p=0.1: tau0 unbiased (3 SE)',  abs(mean(dt))/(std(dt)/sqrt(nT) + eps) , 3);
    r = tc(r, 'BSC p=0.1: nu0 unbiased (3 SE)',   abs(mean(dn))/(std(dn)/sqrt(nT) + eps) , 3);
    fprintf('    (info) pseudo-true jitter std: tau %.3f bins, nu %.3f bins at p = 0.1\n', std(dt), std(dn));

    % ---- affine-in-sigma2 assembly ----
    [A1, ~, Jt1, b1] = mcrb_AB_cp_ofdm(X, Xh, g0, ce, true, 1, cr);
    Jg1 = fim(X, 'cp-ofdm', g0, cr, true, 1);
    pc = mcrb_assemble(A1, Jt1, b1, Jg1, cr.eta_bar, ce.eta_bar, g0);
    N0 = 0.37;
    [A2, B2] = mcrb_AB_cp_ofdm(X, Xh, g0, ce, true, N0, cr);
    lb = mcrb_eval(pc, N0);
    direct = stab_quad_bin(A2, B2, D, true(4*cr.P, 1));
    r = tc(r, 'A^-1 B A^-1 at sigma2 = N0 == N0 G + F', relerr(lb.MCRB0, direct), 1e-10);
    Jg = fim(X, 'cp-ofdm', g0, cr, true, N0);
    r = tc(r, 'CRB piece == fim inverse at N0', relerr(lb.CRB, crb_from_fim(Jg, g0, cr).C), 1e-8);
    [Am, ~, Jtm, bm] = mcrb_AB_cp_ofdm(X, X, g0, cr, true, 1, cr);
    pcm = mcrb_assemble(Am, Jtm, bm, Jg1, cr.eta_bar, cr.eta_bar, g0);
    lbm = mcrb_eval(pcm, N0);
    r = tc(r, 'p = 0: LB == MCRB0 == CRB', max(relerr(lbm.LB, lbm.CRB), relerr(lbm.MCRB0, lbm.CRB)), 1e-10);

    % ---- MT40 runtime gate: one varpro solve at N=128, M=16, P=40, nu masked ----
    gm = make_grid(128, 16, 160e6/128, 127/160e6, 'cp');
    chm = chan_mt40(gm);  crm = chan_realize(chm, gm, []);
    mp = modem('16psk');
    rs = crn_stream(93, 0, 1);
    Xm = gen_symbols(128, 16, '16psk', rs);  Xmh = bsc_symbols(Xm, 0.05, mp, rs);
    for ici = [false true]
        [~, im] = pseudo_true(Xm, Xmh, gm, crm, ici, struct('p', 0.05));
        fprintf('    (info) MT40 varpro P=40 nu-masked ici=%d: %.2f s, %d iterations, residual %.4e\n', ici, im.seconds, im.iterations, im.resnorm);
        r = tc(r, sprintf('MT40 varpro runtime gate (< 10 s, < 150 it) ici=%d', ici), double(im.seconds > 10 || im.iterations > 150), 0);
    end

    % ---- MMLE sanity: single target, matched, high SNR ----
    Ns = 32; Ms = 16;  gs = make_grid(Ns, Ms, df, 0, 'cp');
    cs.P = 1; cs.alpha = 0.9*exp(1j*0.4); cs.tau = 1.3*gs.Ts; cs.nu = 0.3/(Ms*gs.TS); cs.eta_bar = eta_pack(cs.alpha, cs.tau, cs.nu);
    N0 = 1e-3;  nT = 100;  err = zeros(nT, 2);
    Jc = fim(ones(Ns, Ms), 'cp-ofdm', gs, cs, false, N0);
    for t = 1:nT
        rs = crn_stream(94, 0, t);
        Xs = gen_symbols(Ns, Ms, '4qam', rs);
        Ys = io_cp_ofdm(Xs, gs, cs, false) + sqrt(N0/2)*(randn(rs, Ns, Ms) + 1j*randn(rs, Ns, Ms));
        eh = mmle_cp(Ys, Xs, gs, 1, false);
        [~, th, nh] = eta_unpack(eh);
        err(t, :) = [th - cs.tau, nh - cs.nu];
    end
    Sc = crb_from_fim(Jc, gs, cs);
    r = tc(r, 'MMLE matched: delay RMSE <= 1.2 sqrt(CRB)',   sqrt(mean(err(:,1).^2))/sqrt(Sc.crb_tau) - 1.2, 0);
    r = tc(r, 'MMLE matched: Doppler RMSE <= 1.2 sqrt(CRB)', sqrt(mean(err(:,2).^2))/sqrt(Sc.crb_nu) - 1.2, 0);
    fprintf('    (info) MMLE RMSE/sqrt(CRB): delay %.3f, Doppler %.3f\n', sqrt(mean(err(:,1).^2))/sqrt(Sc.crb_tau), sqrt(mean(err(:,2).^2))/sqrt(Sc.crb_nu));
end
