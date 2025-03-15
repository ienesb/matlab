function Gb = getGbv2(nu, tau, phi, b)
    Na = 64;
    Nrf = 4;
    Ns = 1;
    B = 6;
    M = 64;          % number of subcarriers
    N = 30;          % number of subsymbols/frame
    df = 15e3;       % make this the frequency bin spacing of LTE
    T = 1/df;

    U = getU(Na, Nrf, B, b);
    F = getF(Na, Ns);
    
    a = array_responsev2(phi, Na);
    a = permute(a, [1,3,2]);
    aH = conj(permute(a, [2,1,3]));
    
    aaH = pagemtimes(a, aH);
    
    C = pagemtimes(U', pagemtimes(aaH, F));

    ddm = dd_crosstalk_coefficientsv3(nu, tau, T, N, M);
    ddm = permute(ddm, [1,2,5,3,4]);
    Gb = my_kron_nd(C, ddm);
    
end