function r = test_channels()
% Channel specifications (Phase 0.5, channels.md).
    r = [];
    N = 512; M = 128; df = 15e3;  Ts = 1/(N*df);
    g = make_grid(N, M, df, 20*Ts, 'cp');
    ch = chan_eva(g);
    r = tc(r, 'EVA taps [0 0 1 2 3 5 8 13 19]', double(~isequal(ch.l_tap.', [0 0 1 2 3 5 8 13 19])), 0);
    r = tc(r, 'EVA l_tau = 20',                 abs(ch.l_tau - 20), 0);
    r = tc(r, 'EVA k_nu = 4',                   abs(ch.k_nu - 4), 0);
    r = tc(r, 'EVA sum P_l = 1',                abs(sum(ch.P_lin) - 1), 1e-15);
    r = tc(r, 'EVA linear power sum 4.14594',   abs(sum(10.^(ch.pow_db/10)) - 4.14594)/4.14594, 1e-5);
    r = tc(r, 'EVA nu_max = 444.44 Hz',         abs(ch.nu_max - 444.444)/444.444, 1e-5);
    r = tc(r, 'EVA collisions {[1 2]} (and [4 5]? no: 310->2, 370->3)', ...
        double(~(numel(ch.meta.collisions) == 1 && isequal(ch.meta.collisions{1}, [1 2]))), 0);
    r = tc(r, 'EVA Tcp = 2.604 us',             abs(g.Tcp - 2.6042e-6)/2.6042e-6, 1e-4);

    % draws
    rs = crn_stream(9000, 0, 1);
    cr = chan_realize(ch, g, rs);
    r = tc(r, 'EVA realization size P = 9',     abs(cr.P - 9), 0);
    r = tc(r, 'EVA snapped tau = l_tap Ts',     relerr(cr.tau, ch.l_tap*Ts), 1e-15);
    r = tc(r, 'EVA |nu| <= nu_max',             double(any(abs(cr.nu) > ch.nu_max + 1e-9)), 0);
    % Rayleigh gain power over many draws ~ P_l
    K = 4000;  acc = zeros(9, 1);
    for t = 1:K, c = chan_realize(ch, g, crn_stream(9000, 0, t)); acc = acc + abs(c.alpha).^2; end
    r = tc(r, 'EVA Rayleigh E|alpha|^2 == P_l (4000 draws)', relerr(acc/K, ch.P_lin), 0.08);
    chc = chan_eva(g, 'gain_model', 'constmag');
    c = chan_realize(chc, g, rs);
    r = tc(r, 'EVA constmag |alpha|^2 == P_l exactly', relerr(abs(c.alpha).^2, chc.P_lin), 1e-14);
    chr = chan_eva(g, 'reduced', true);
    r = tc(r, 'EVA-reduced taps [0 1 3]',        double(~isequal(chr.l_tap.', [0 1 3])), 0);
    chm = chan_eva(g, 'merge_colliding', true);
    cm = chan_realize(chm, g, rs);
    r = tc(r, 'merge_colliding -> P = 8',        abs(cm.P - 8), 0);

    % legacy preset pairs (tau, theta, phi) by index like old/eva/verify_crb.m
    gl = make_grid(16, 8, df, 0, 'cp');
    chl = chan_eva(gl, 'legacy_preset', true);
    cl  = chan_realize(chl, gl, []);
    numax = (120/3.6)*2e9/3e8;
    alpha_ref = 10.^([0; -0.6; -1.4]/20) .* exp(1j*[0.7; 1.1; -0.5]);
    nu_ref    = numax*cos([0.5; 2.0; 3.5]);
    r = tc(r, 'legacy preset alpha',  relerr(cl.alpha, alpha_ref), 1e-15);
    r = tc(r, 'legacy preset nu',     relerr(cl.nu, nu_ref), 1e-15);
    r = tc(r, 'legacy preset tau [0 370 150] ns', relerr(cl.tau, [0; 370; 150]*1e-9), 1e-15);

    % MT40 geometry
    gm = make_grid(128, 1, 160e6/128, 0, 'cp');
    cmt = chan_mt40(gm);
    rng(20260514);  u = sort(88*rand(40,1));  tau_ref = u + (0:39).';
    r = tc(r, 'MT40 geometry == rng(20260514) draw', relerr(cmt.tau_samp, tau_ref), 0);
    r = tc(r, 'MT40 min separation >= 1 sample',   double(min(diff(cmt.tau_samp)) < 1 - 1e-12), 0);
    r = tc(r, 'MT40 max delay <= 127 samples',     double(max(cmt.tau_samp) > 127), 0);
    r = tc(r, 'MT40 range unit 0.9375 m',          abs(cmt.range_unit_m - 0.9375), 1e-12);
    crm = chan_realize(cmt, gm, []);
    r = tc(r, 'MT40 unit gains, nu = 0',           double(any(crm.alpha ~= 1) || any(crm.nu ~= 0)), 0);
    gm2 = make_grid(128, 16, 160e6/128, 127/160e6, 'cp');
    cmd = chan_mt40(gm2, 'doppler', 'mt40d');
    r = tc(r, 'MT40-D nu_max = 0.5/(M TS)',        abs(cmd.nu_max - 0.5/(16*gm2.TS))/cmd.nu_max, 1e-12);
    r = tc(r, 'MT40-D |nu| <= nu_max',             double(any(abs(cmd.nu) > cmd.nu_max)), 0);
end
