function H = getChannelMatrix(N, M, tau, nu, alpha, deltaf, Ts) % ICI
    T = 1 / deltaf;

    ns = 0:(N-1);
    nps = (0:(N-1))';
    ms = 0:(M-1);

    Hb = exp(-1j*2*pi*nps*deltaf*tau) .* get_ici_coeffs(nps - ns, nu, T);

    H = alpha * kron(diag(exp(1j*2*pi*ms*Ts*nu)), Hb);
end