function r = test_crb_units()
% Conditioned CRB extraction (Phase 0.6): bin-unit route vs Jacobi route,
% Schur/EFI vs block of the inverse, structurally null nu block, kappa invariance,
% conditioning of the full 9-path EVA FIM.
    r = [];
    N = 512; M = 128; df = 15e3;  Ts = 1/(N*df);
    g = make_grid(N, M, df, 20*Ts, 'cp');
    ch = chan_eva(g, 'reduced', true);
    cr = chan_realize(ch, g, crn_stream(9000, 0, 1));
    X = gen_symbols(N, M, '4qam', crn_stream(9001, 0, 1));
    J = fim(X, 'cp-ofdm', g, cr, false, 1);
    Sb = crb_from_fim(J, g, cr, 'method', 'bin');
    Sj = crb_from_fim(J, g, cr, 'method', 'jacobi');
    r = tc(r, 'bin route == jacobi route (EVA-reduced 512x128)', relerr(Sb.C, Sj.C), 1e-8);
    Se = crb_from_fim(J, g, cr, 'method', 'bin', 'efi', true);
    idx = eta_index(cr.P);
    r = tc(r, 'EFI/Schur == kinematic block of the full inverse', relerr(Se.C(idx.kin, idx.kin), Sb.C(idx.kin, idx.kin)), 1e-10);
    r = tc(r, 'SI cond ~ 1e20 vs bin cond < 1e8', double(~(cond(J) > 1e12 && Sb.kappa < 1e8)), 0);
    fprintf('    cond(J_SI) = %.2e, cond(J_bin) = %.2e\n', cond(J), Sb.kappa);

    % nu mask: nu = 0 no-ICI FIM -> masked inverse == 3P inverse
    cz = cr;  cz.nu(:) = 0;  cz.eta_bar = eta_pack(cz.alpha, cz.tau, cz.nu);
    Jz = fim(X, 'cp-ofdm', g, cz, false, 1);
    Sz = crb_from_fim(Jz, g, cz);
    keep = true(4*cr.P, 1);  keep(idx.nu) = false;
    D = eta_scale(g, cr.P);  Jb = D*Jz*D;
    C3 = D(keep, keep) * inv(Jb(keep, keep)) * D(keep, keep);
    r = tc(r, 'nu mask == 3P inverse (nu = 0)',     relerr(Sz.C(keep, keep), C3), 1e-10);
    r = tc(r, 'nu mask -> crb_nu NaN',              double(~all(isnan(Sz.crb_nu))), 0);
    r = tc(r, 'nu mask kept by keep_nu',            double(any(crb_from_fim(Jz, g, cz, 'keep_nu', true).active == false)), 0);

    % kappa invariant to sigma2 and |alpha| scaling
    J2 = fim(X, 'cp-ofdm', g, cr, false, 3.7);
    cs = cr;  cs.alpha = 2.5*cs.alpha;  cs.eta_bar = eta_pack(cs.alpha, cs.tau, cs.nu);
    J3 = fim(X, 'cp-ofdm', g, cs, false, 1);
    S2 = crb_from_fim(J2, g, cr);  S3 = crb_from_fim(J3, g, cs);
    r = tc(r, 'CRB scales with sigma2',             relerr(S2.C, 3.7*Sb.C), 1e-10);
    r = tc(r, 'kappa invariant to sigma2',          abs(S2.kappa - Sb.kappa)/Sb.kappa, 1e-8);
    idxk = eta_index(cr.P).kin;
    r = tc(r, 'kinematic CRB scales with 1/|alpha|^2', relerr(S3.C(idxk, idxk), Sb.C(idxk, idxk)/2.5^2), 1e-10);

    % stab_quad_bin consistency: A = -J, C = J  ->  A^-1 C A^-1 = J^-1
    Y = stab_quad_bin(-J, J, D, true(4*cr.P, 1));
    r = tc(r, 'stab_quad_bin(-J, J) == CRB', relerr(Y, Sb.C), 1e-8);

    % full 9-path EVA with ICI: conditioning of the bin-unit FIM over a few draws
    ch9 = chan_eva(g);  kmax = 0;  bad = 0;
    for t = 1:3
        c9 = chan_realize(ch9, g, crn_stream(9000, 1, t));
        J9 = fim(X, 'cp-ofdm', g, c9, true, 1);
        S9 = crb_from_fim(J9, g, c9);
        kmax = max(kmax, S9.kappa);  bad = bad + ~S9.valid;
    end
    fprintf('    full EVA (ICI) max cond(J_bin) over 3 draws = %.2e\n', kmax);
    r = tc(r, 'full EVA 9-path FIM invertible (3 draws)', bad, 0);
end
