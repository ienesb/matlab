function Gb = getGb(nu, tau, phi, b, B)
    Na = 64;
    Nrf = 4;
    Ns = 1;
    M = 64;          % number of subcarriers
    N = 30;          % number of subsymbols/frame
    df = 15e3;       % make this the frequency bin spacing of LTE
    T = 1/df;

    U = getU(Na, Nrf, B, b);
    F = getF(Na, Ns);
    
    a = array_response(phi, Na);

    dd = dd_crosstalk_coefficientsv2(nu, tau, T, N, M);
    ddm = dd2m(dd, M, N);

    C = U' * (a * a') * F;
    Gb = kron(C, ddm);
end