function CRB = getCRB_OFDM(X, deltaf, T, sigma2, alpha, nu, tau, mask) % monostatic [alpha_R, alpha_I, nu, tau]
    [N, M] = size(X);

    X(~mask) = 0;

    ns = (0:(N-1)).';
    ms = 0:(M-1);

    c = (exp(-1j*2*pi*deltaf*ns*tau) * exp(1j*2*pi*T*ms*nu)) .* X;
    c = c(:);

    m = ones(N, 1) * ms;
    m = m(:);
    n = ns * ones(1, M);
    n = n(:);
    
    g0 = c;
    g1 = 1j * c;
    g2 = 1j * 2 * pi * T * alpha .* m .* c;
    g3 = -1j * 2 * pi * deltaf * alpha.* n .* c;
    
    G  = [g0, g1, g2, g3];             % NM x 4
    J  = (2/sigma2) * real(G' * G);    % exact J_like
    CRB = inv(J);
end