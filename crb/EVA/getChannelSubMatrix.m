function Hb = getChannelSubMatrix(N, m, tau, nu, alpha, deltaf, Ts) % ICI
    T = 1 / deltaf;

    ns = 0:(N-1);
    nps = (0:(N-1))';

    Hb = exp(-1j*2*pi*ns*deltaf*tau) .* get_ici_coeffs(nps - ns, nu, T);
    
    Hb = alpha * exp(1j*2*pi*m*Ts*nu) * Hb;
end