function r = test_rx()
% Phase-2 receiver tests: DD tap extraction vs the Dirichlet formula, MP
% detector on AWGN / one tap, posterior mode vs the legacy detector, frame
% masks and pilot fairness, OFDM CE vs the legacy copy, BSC moments.
    r = [];
    oracle_paths('add', 'old');  c = onCleanup(@() oracle_paths('remove', 'old'));
    N = 64; M = 16; df = 15e3;  Ts = 1/(N*df);
    g = make_grid(N, M, df, 0, 'cp');            % Tcp = 0 so TS = T (legacy Dirichlet formula uses T)
    mdm = modem('4qam');

    % ---- dd_response_io: integer taps -> single entries; fractional Doppler row
    cr.P = 2;  cr.alpha = [0.9*exp(1j*0.2); 0.5*exp(-1j*1.0)];  cr.tau = [1; 3]*Ts;  cr.nu = [2; -1]/(M*g.TS);
    cr.eta_bar = eta_pack(cr.alpha, cr.tau, cr.nu);  cr.k_real = cr.nu*M*g.TS;
    taps = dd_response_io(g, cr, false, 10, 1e-9, 'origin');
    r = tc(r, 'dd_response_io: integer (l,k) taps -> 2 taps',        abs(taps.S - 2), 0);
    [~, i1] = min(abs(taps.gain - cr.alpha(1)));
    r = tc(r, 'dd_response_io: gain/shift of tap 1',                  abs(taps.gain(i1) - cr.alpha(1)) + abs(taps.dl(i1) - 1) + abs(mod(taps.dk(i1), M) - 2), 1e-12);
    crf = cr;  crf.nu = [2.3; -0.6]/(M*g.TS);  crf.k_real = crf.nu*M*g.TS;
    tf = dd_response_io(g, crf, false, Inf, 0, 'origin');
    Gold = build_effective_dd_eva(crf.alpha, round(crf.tau/Ts), crf.nu, N, M, df, Inf);
    % Relation to the legacy Raviteja eq. (13) row (TS = T here): remove the
    % twist e^{-j2pi k_real l/(NM)}, and conjugate the Dirichlet leakage factor
    % (the io layer, derived from cp_otfs_io.tex and checked against the
    % explicit double sums, gives sum_m e^{j2pi m (k_real - k')/M}; the printed
    % eq. (13) carries the conjugate argument).  Magnitudes are identical.
    Gt = Gold;
    for p = 1:crf.P
        lp = round(crf.tau(p)/Ts);
        row = Gt(lp+1, :) * exp(1j*2*pi*crf.k_real(p)*lp/(N*M));      % twist removed
        Gt(lp+1, :) = crf.alpha(p) * conj(row / crf.alpha(p));          % conjugate leakage factor
    end
    r = tc(r, 'dd_response_io |h| == legacy Dirichlet row |G| (TS=T)',       relerr(abs(tf.h), abs(Gold)), 1e-10);
    r = tc(r, 'dd_response_io == legacy row (twist removed, leakage phase conjugated)', relerr(tf.h, Gt), 1e-10);
    te = dd_response_io(g, crf, false, 10, 0, 'origin');
    r = tc(r, 'E = 10 keeps > 99 % of the tap energy',                 1 - te.energy_kept, 0.01);
    gi = make_grid(N, M, df, 3*Ts, 'cp');  cri = crf;  cri.nu = [0.3; -0.4]*444/(M*gi.TS)*M*gi.TS/M;  % small Doppler
    cri.nu = [40; -60];  cri.k_real = cri.nu*M*gi.TS;
    tp = dd_response_io(gi, cri, true, 10, 0, 'per_row');
    fprintf('    (info) ICI shift-variance residual of the DD response at nu T ~ 0.004: %.2e (O(nu T) leakage the MP ignores)\n', tp.residual);
    r = tc(r, 'per_row residual O(nu T) with ICI (nu T ~ 0.004)',        tp.residual, 5e-2);

    % ---- MP detector: one tap == slicer; AWGN BER vs theory (small grid)
    cr1.P = 1;  cr1.alpha = 0.8*exp(1j*0.7);  cr1.tau = 2*Ts;  cr1.nu = 3/(M*g.TS);  cr1.eta_bar = eta_pack(cr1.alpha, cr1.tau, cr1.nu);  cr1.k_real = 3;
    rs = crn_stream(77, 0, 1);
    x = gen_symbols(N, M, '4qam', rs);
    y = io_cp_otfs(x, g, cr1, false) + sqrt(0.05/2)*(randn(rs, N, M) + 1j*randn(rs, N, M));
    t1 = dd_response_io(g, cr1, false, 10, 1e-9, 'origin');
    xh = mp_detect_dd(y, t1, mdm.alphabet, 0.05, false(N, M), zeros(N, M), 10, 0.7, 'edge');
    xs = reshape(mdm.map(mdm.demap(circshift(y, [-2, -3])/cr1.alpha)), N, M);
    r = tc(r, 'MP one tap == slicer',                                  double(~isequal(xh, xs)), 0);
    % AWGN BER at 6 dB, 200 frames 64 x 16, both MP modes
    ch = chan_controls(g, 'awgn');  nb = 0; ne = 0; nep = 0;
    for t = 1:200
        rs = crn_stream(78, 0, t);  cra = chan_realize(ch, g, rs);  cra.k_real = 0;
        [x, bits] = gen_symbols(N, M, '4qam', rs);
        N0 = 10^(-0.6);
        y = cra.alpha*x + sqrt(N0/2)*(randn(rs, N, M) + 1j*randn(rs, N, M));
        ta = struct('dl', 0, 'dk', 0, 'gain', cra.alpha);
        xh = mp_detect_dd(y, ta, mdm.alphabet, N0, false(N, M), zeros(N, M), 10, 0.7, 'edge');
        xp = mp_detect_dd(y, ta, mdm.alphabet, N0, false(N, M), zeros(N, M), 10, 0.3, 'posterior');
        ne = ne + nnz(mdm.demap(xh(:)) ~= bits);  nep = nep + nnz(mdm.demap(xp(:)) ~= bits);  nb = nb + numel(bits);
    end
    th = 0.5*erfc(sqrt(10^0.6/2));
    r = tc(r, 'MP edge-mode AWGN BER at 6 dB vs Q(sqrt(SNR))',           abs(ne/nb/th - 1), 0.10);
    r = tc(r, 'MP posterior-mode AWGN BER at 6 dB vs Q(sqrt(SNR))',      abs(nep/nb/th - 1), 0.10);
    % posterior mode bit-exact vs legacy mp_detector on a 2-tap channel with a pilot frame
    chE = chan_eva(g, 'reduced', true);  chE.l_tau = 3;
    rs = crn_stream(79, 0, 1);  crE = chan_realize(chE, g, rs);
    [xd] = gen_symbols(N, M, '4qam', rs);
    [xf, mk] = otfs_pilot_frame(xd, g, chE, 5);
    tE = dd_response_io(g, crE, false, 4, 0.02, 'origin');
    yE = io_cp_otfs(xf, g, crE, false) + sqrt(0.1/2)*(randn(rs, N, M) + 1j*randn(rs, N, M));
    xn = mp_detect_dd(yE, tE, mdm.alphabet, 0.1, mk.known, mk.known_val, 8, 0.3, 'posterior');
    xo = mp_detector(yE, tE.dl, tE.dk, tE.gain, 0.1, mdm.alphabet, 8, mk.data, mk.known_val, 0.3);
    xo(mk.known) = mk.known_val(mk.known);
    r = tc(r, 'posterior mode == legacy mp_detector (data cells)',       double(~isequal(xn(mk.data), xo(mk.data))), 0);

    % ---- frames and fairness on the headline EVA grid
    N2 = 512; M2 = 128;  g2 = make_grid(N2, M2, df, 20/(N2*df), 'cp');  ch2 = chan_eva(g2);
    [xf, mk] = otfs_pilot_frame(zeros(N2, M2), g2, ch2, 3);
    r = tc(r, 'OTFS full guard 5248 cells (8.0 %)',                      abs(mk.n_nondata - 5248), 0);
    r = tc(r, 'OTFS pilot at (l_tau, M/2)',                             abs(xf(21, 65) - 3), 0);
    [~, mo, a_p, info] = ofdm_pilot_lattice(zeros(N2, M2), g2, ch2, mk.n_nondata, mk.n_nondata);
    fprintf('    (info) OFDM lattice (dN,dM) = (%d,%d), %d pilots vs target %d: mismatch %.1f%%\n', mo.dN, mo.dM, mo.n_pilot, mk.n_nondata, 100*info.mismatch);
    r = tc(r, 'OFDM lattice count within 5 % of the OTFS budget',        info.mismatch, 0.05);
    r = tc(r, 'OFDM pilot energy == OTFS pilot energy',                 abs(mo.n_pilot*a_p^2 - mk.n_nondata)/mk.n_nondata, 1e-12);
    r = tc(r, 'dft CE Nyquist l_tau < N/dN',                            double(~(ch2.l_tau < N2/mo.dN)), 0);

    % ---- pinned OFDM lattice (speed sweep: the search result depends on k_nu)
    chR = chan_eva(g2, 'reduced', true);  ch60 = chan_eva(g2, 'reduced', true, 'v_kmh', 60);
    nR = (2*chR.l_tau + 1)*M2;  Xr = gen_symbols(N2, M2, '4qam', crn_stream(82, 0, 1));
    [dA, mA] = ofdm_pilot_lattice(Xr, g2, chR, nR, nR);
    [dB, mB] = ofdm_pilot_lattice(Xr, g2, chR, nR, nR, [], [37 2]);
    r = tc(r, 'EVA-reduced 120 km/h: search gives lattice (37, 2)',     abs(mA.dN - 37) + abs(mA.dM - 2), 0);
    r = tc(r, 'pinned (37, 2) == searched lattice at 120 km/h',          double(~(isequal(dA, dB) && isequal(mA, mB))), 0);
    [~, m6] = ofdm_pilot_lattice(Xr, g2, ch60, nR, nR);
    [~, m6p] = ofdm_pilot_lattice(Xr, g2, ch60, nR, nR, [], [37 2]);
    r = tc(r, 'EVA-reduced 60 km/h: search gives lattice (4, 19)',       abs(m6.dN - 4) + abs(m6.dM - 19), 0);
    r = tc(r, 'pinned (37, 2) at 60 km/h: same pilot mask as 120 km/h',  double(~(isequal(m6p.pilot, mA.pilot) && m6p.n_pilot == m6.n_pilot)), 0);

    % ---- OFDM CE == legacy ce_ofdm (via mcrb_ce_ici local copy re-implemented here) on a mid grid
    N3 = 128; M3 = 16;  g3 = make_grid(N3, M3, df, 0, 'cp');  ch3 = chan_eva(g3, 'reduced', true);
    rs = crn_stream(80, 0, 1);  cr3 = chan_realize(ch3, g3, rs);
    Xd = gen_symbols(N3, M3, '4qam', rs);
    [d3, m3, ~] = ofdm_pilot_lattice(Xd, g3, ch3, 2*(ch3.l_tau+1)*M3, 2*(ch3.l_tau+1)*M3);
    N0 = 0.01;  Y3 = io_cp_ofdm(d3, g3, cr3, true) + sqrt(N0/2)*(randn(rs, N3, M3) + 1j*randn(rs, N3, M3));
    Hg = io_cp_ofdm(ones(N3, M3), g3, cr3, false);                      % genie diagonal channel
    for mode = {'interp', 'dft', 'dftdd'}
        H = ofdm_ce(Y3, d3, m3, ch3.l_tau, N0, mdm, mode{1}, 2);
        Hl = ce_ofdm_legacy(Y3, d3, m3.pilot, m3.n_pilots, m3.m_pilots, N3, M3, ch3.l_tau, N0, mode{1}, 2, mdm);
        r = tc(r, sprintf('ofdm_ce %-6s == legacy ce_ofdm', mode{1}), relerr(H, Hl), 1e-12);
        fprintf('    (info) %-6s CE MSE vs genie: %.2e\n', mode{1}, mean(abs(H(:) - Hg(:)).^2)/mean(abs(Hg(:)).^2));
    end
    Hd = ofdm_ce(Y3, d3, m3, ch3.l_tau, N0, mdm, 'dftdd', 2);
    r = tc(r, 'dftdd CE relative MSE vs genie < 1 % at 20 dB',       mean(abs(Hd(:) - Hg(:)).^2)/mean(abs(Hg(:)).^2), 0.01);

    % ---- BSC moments
    rs = crn_stream(81, 0, 1);  X = gen_symbols(200, 200, '4qam', rs);  p = 0.1;
    Xh = bsc_symbols(X, p, mdm, rs, 'bit');  E = X - Xh;
    r = tc(r, 'BSC: E[Xhat|X] = (1-2p) X',       abs(mean(Xh(:).*conj(X(:))) - (1-2*p)), 0.01);
    r = tc(r, 'BSC: E[E* Xhat] = -2p',            abs(mean(conj(E(:)).*Xh(:)) + 2*p), 0.01);
    r = tc(r, 'BSC: E|E|^2 = 4p',                 abs(mean(abs(E(:)).^2) - 4*p), 0.02);
    Xb = bsc_match(X, Xh, true(size(X)), mdm, rs);
    r = tc(r, 'bsc_match preserves the error count', abs(nnz(Xb ~= X) - nnz(Xh ~= X)), 0);
    r = tc(r, 'mp memory estimate for full EVA < 1 GB (S=189, Q=4)', double(8*512*128*189*4 > 1e9), 0);
end

function H_hat = ce_ofdm_legacy(Y, d, pilot_mask, n_pilots, m_pilots, N, M, l_tau, N0, ce_mode, n_dd, mdm)
% verbatim logic of the old ce_ofdm local (qammod/qamdemod replaced by the modem)
    H_ls = Y(n_pilots, m_pilots) ./ d(n_pilots, m_pilots);
    [Nq, Mq] = ndgrid(1:N, 1:M);
    Afull = exp(-1j*2*pi * (0:N-1).' * (0:l_tau) / N);
    switch ce_mode
        case 'interp'
            F_H   = griddedInterpolant({n_pilots, m_pilots}, H_ls, 'linear','nearest');
            H_hat = F_H(Nq, Mq);  return
        case {'dft','dftdd'}
            Abas  = exp(-1j*2*pi * (n_pilots(:)-1) * (0:l_tau) / N);
            Hp    = Afull * (Abas \ H_ls);
            F_H   = griddedInterpolant({(1:N).', m_pilots}, Hp, 'linear','nearest');
            H_hat = F_H(Nq, Mq);
    end
    if strcmp(ce_mode,'dftdd')
        for it = 1:n_dd
            Xs   = conj(H_hat).*Y ./ (abs(H_hat).^2 + N0);
            dhat = reshape(mdm.map(mdm.demap(Xs(:))), N, M);
            dhat(pilot_mask) = d(pilot_mask);
            H_hat = Afull * (Afull \ (Y ./ dhat));
        end
    end
end
