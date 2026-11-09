function r = test_mcrb_expected()
% Port of old/eva/verify_mcrb.m ('full' tier): Monte-Carlo A and B vs the
% closed-form expectations, OFDM and CP-OTFS, {no-ICI, ICI} x {P = 1, 3}, on the
% legacy preset with p = 0.1.  E{A} agrees across waveforms; E{B} agrees within
% each waveform (Cov(b) differs across them by construction -- reported).
    r = [];
    oracle_paths('add', 'old');  c = onCleanup(@() oracle_paths('remove', 'old'));
    N = 16; M = 8; df = 15e3;  sigma2 = 1.3;  p = 0.1;  nT = 40000;
    g0 = make_grid(N, M, df, 0, 'cp');
    ch = chan_eva(g0, 'legacy_preset', true);  cr = chan_realize(ch, g0, []);
    mdm = modem('4qam');
    for ici = [false true]
        for P = [1 3]
            cp.P = P; cp.alpha = cr.alpha(1:P); cp.tau = cr.tau(1:P); cp.nu = cr.nu(1:P); cp.eta_bar = eta_pack(cp.alpha, cp.tau, cp.nu);
            Ao = 0; Bo = 0; At = 0; Bt = 0;
            rs = crn_stream(95, ici + 2*(P == 3), 1);
            for t = 1:nT
                X = gen_symbols(N, M, '4qam', rs);  Xh = bsc_symbols(X, p, mdm, rs);
                x = gen_symbols(N, M, '4qam', rs);  xh = bsc_symbols(x, p, mdm, rs);
                [A, B] = mcrb_AB_cp_ofdm(X, Xh, g0, cp, ici, sigma2);   Ao = Ao + A;  Bo = Bo + B;
                [A, B] = mcrb_AB_cp_otfs(x, xh, g0, cp, ici, sigma2);   At = At + A;  Bt = Bt + B;
            end
            Ao = Ao/nT; Bo = Bo/nT; At = At/nT; Bt = Bt/nT;
            [Aeo, Beo] = mcrb_AB_expected('cp-ofdm', g0, cp, ici, sigma2, p);
            [Aet, Bet] = mcrb_AB_expected('cp-otfs', g0, cp, ici, sigma2, p);
            if ici, [Aold, Bold] = getAB_OFDM_multi_ICI_expected(N, M, df, g0.TS, sigma2, cp.alpha, cp.tau, cp.nu, p);
            else,   [Aold, Bold] = getAB_OFDM_multi_noICI_expected(N, M, df, g0.TS, sigma2, cp.alpha, cp.tau, cp.nu, p); end
            r = tc(r, sprintf('expected A,B == old (ici=%d, P=%d)', ici, P), max(relerr(Aeo, Aold), relerr(Beo, Bold)), 1e-12);
            r = tc(r, sprintf('OFDM MC A == E{A} (ici=%d, P=%d)', ici, P), relerr(Ao, Aeo), 3e-2);
            r = tc(r, sprintf('OTFS MC A == E{A} (ici=%d, P=%d)', ici, P), relerr(At, Aet), 3e-2);
            r = tc(r, sprintf('E{A} OFDM == OTFS (ici=%d, P=%d)', ici, P), relerr(Aeo, Aet), 1e-12);
            r = tc(r, sprintf('OFDM MC B == E{B} incl. Cov(b) (ici=%d, P=%d)', ici, P), relerr(Bo, Beo), 3e-2);
            r = tc(r, sprintf('OTFS MC B == E{B} incl. Cov(b) (ici=%d, P=%d)', ici, P), relerr(Bt, Bet), 3e-2);
            fprintf('    (info) E{B} OFDM vs OTFS gap (Cov(b) differs): %.3e; MC B OFDM vs OTFS: %.3e\n', relerr(Beo, Bet), relerr(Bo, Bt));
            if ~ici && P == 1
                % no-ICI single target: Cov(b) = 4p(1-p)|alpha|^2 J / sigma2 for OFDM
                Cb = covb_cp_ofdm(g0, cp, false, sigma2, p);
                J  = fim_expected('cp-ofdm', g0, cp, false, sigma2);
                r = tc(r, 'OFDM no-ICI Cov(b) == 4p(1-p)|a|^2 J/sigma2', relerr(Cb, 4*p*(1-p)*abs(cp.alpha)^2/sigma2*J), 1e-10);
                Cbt = covb_cp_otfs(g0, cp, false, sigma2, p);
                idx = eta_index(1);
                fprintf('    (info) OTFS no-ICI Cov(b): gain block vs collapse %.2e; kinematic block vs collapse %.2e\n', ...
                    relerr(Cbt(idx.gain, idx.gain), 4*p*(1-p)*abs(cp.alpha)^2/sigma2*J(idx.gain, idx.gain)), ...
                    relerr(Cbt(idx.kin, idx.kin),  4*p*(1-p)*abs(cp.alpha)^2/sigma2*J(idx.kin, idx.kin)));
                % mean-field identity E{A} = -(1-2p) J on kinematic and gain-kinematic entries (no ICI)
                Ea = Aeo;  chk = -(1-2*p)*J;  mask = true(4);  mask(1:2, 1:2) = false;
                r = tc(r, 'E{A} = -(1-2p) J on kinematic/gain-kinematic entries (no ICI)', relerr(Ea(mask), chk(mask)), 1e-10);
                r = tc(r, 'E{A} = -J on the gain-gain block', relerr(Ea(1:2,1:2), -J(1:2,1:2)), 1e-10);
            end
        end
    end
end
