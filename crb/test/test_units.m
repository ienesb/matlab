function r = test_units()
% Grid constants and the exact bin-unit change of variables (Phase 0.1).
    r = [];
    N = 512; M = 128; df = 15e3;
    Ts = 1/(N*df);
    g  = make_grid(N, M, df, 20*Ts, 'cp');
    r = tc(r, 'cp: TS = T + Tcp',        abs(g.TS - (1/df + 20*Ts)), 1e-18);
    r = tc(r, 'cp: bin_nu = 1/(M TS)',   abs(g.bin_nu - 1/(M*g.TS)), 1e-12);
    r = tc(r, 'cp: Doppler bin 112.8 Hz', abs(g.bin_nu - 112.79)/112.79, 2e-4);
    gr = make_grid(N, M, df, 20*Ts, 'rcp');
    r = tc(r, 'rcp: TS = T (frame CP)',   abs(gr.TS - 1/df), 1e-18);
    r = tc(r, 'rcp: Tframe = M T + Tcp',  abs(gr.Tframe - (M/df + 20*Ts)), 1e-15);

    % D*(D\J/D)*D == J exactly for a random SPD J
    P = 3;  D = eta_scale(g, P);
    A = randn(4*P) + 1j*randn(4*P);  J = real(A*A');
    r = tc(r, 'D*(D\J/D)*D == J', relerr(D*(D\J/D)*D, J), 1e-13);

    % pack / unpack / index
    alpha = randn(P,1) + 1j*randn(P,1);  tau = rand(P,1)*1e-6;  nu = randn(P,1)*100;
    [a2, t2, n2] = eta_unpack(eta_pack(alpha, tau, nu));
    r = tc(r, 'eta_pack/unpack round trip', relerr([a2; t2; n2], [alpha; tau; nu]), 1e-15);
    idx = eta_index(P);
    r = tc(r, 'eta_index tau slots', double(~isequal(idx.tau, [3 7 11])), 0);
end
