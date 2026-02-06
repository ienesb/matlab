function CRB = getCRB_OFDM_ICI(X, deltaf, T, fc, sigma2, alpha, nu, tau, mask) % monostatic [alpha_R, alpha_I, nu, tau]
    [N, M] = size(X);

    X(~mask) = 0;

    ns = (0:(N-1)).';
    ms = 0:(M-1);

    b = exp(-1j*2*pi*deltaf*ns*tau);
    c = exp(-1j*2*pi*T*ms*nu);
    D = diag(exp(1j*2*pi*T/N*ns));
    F = dftmtx(N) / sqrt(N);
    % D = eye(N);
    % F = eye(N);

    C = D * F' * ((b * c) .* X);
    C = C(:);

    m = ones(N, 1) * ms;
    m = m(:);
    n = ns * ones(1, M);
    
    C2 = D * F' * ((-1j*2*pi*deltaf*n) .* (b * c) .* X);
    C2 = C2(:);

    n = n(:);

    g0 = C;
    g1 = 1j * C;
    g2 = 1j * 2 * pi * alpha * T * (m + n/N) .* C;
    g3 = alpha * C2;
    
    G  = [g0, g1, g2, g3];             % NM x 4
    J  = (2/sigma2) * real(G' * G);    % exact J_like
    CRB = inv(J);
end