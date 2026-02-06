function CRB = getCRB_OTFS (X_DD, deltaf, T, sigma2, alpha, nu, tau, mask) % monostatic [alpha_R, alpha_I, nu, tau]
    [N, M] = size(X_DD);

    X_DD(~mask) = 0;

    ns = (0:(N-1)).';
    ms = 0:(M-1);
    
    C_FT = exp(-1j*2*pi*deltaf*ns*tau) * exp(1j*2*pi*T*ms*nu);
    C_tau_FT = (-1j*2*pi*deltaf*ns) .* C_FT;
    C_nu_FT = (1j*2*pi*T*ms) .* C_FT;
    
    % H_FT = alpha * C_FT;

    C_DD = sfft(C_FT);
    C_tau_DD = sfft(C_tau_FT);
    C_nu_DD = sfft(C_nu_FT);
    % H_DD = sfft(H_FT);

    g0 = circular_conv2d(C_DD, X_DD);
    g1 = 1j * circular_conv2d(C_DD, X_DD);
    g2 = alpha * circular_conv2d(C_nu_DD, X_DD);
    g3 = alpha * circular_conv2d(C_tau_DD, X_DD);

    g0 = g0(:);
    g1 = g1(:);
    g2 = g2(:);
    g3 = g3(:);

    G = [g0 g1 g2 g3];
    FIM = 2/sigma2 * real(G' * G);
    CRB = inv(FIM);
end