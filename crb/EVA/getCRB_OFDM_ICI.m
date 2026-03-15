function CRB = getCRB_OFDM_ICI(X, deltaf, T, sigma2, alpha, nu, tau, mask) % monostatic [alpha_R, alpha_I, nu, tau]
    [N, M] = size(X);

    X(~mask) = 0;

    ns = (0:(N-1)).';
    ms = 0:(M-1);

    m_array = ones(N, 1) * ms;
    m_array = m_array(:);
    n_array = ns * ones(1, M);
    n_array = n_array(:);

    % c = (exp(-1j*2*pi*deltaf*ns*tau) * exp(1j*2*pi*T*ms*nu)) .* X;
    % c = channel_noiseless(X, 1, tau, nu, deltaf, T);
    g0 = zeros(N, M);
    g2 = zeros(N, M);
    g3 = zeros(N, M);
    for n = 0:(N-1)
        I = get_ici_coeffs(ns-n, nu, T);
        dI = get_ici_coeffs_derivative(ns-n, nu, T);

        temp = (exp(-1j*2*pi*deltaf*ns*tau) * exp(1j*2*pi*T*ms*nu));
        temp = X .* temp;
        temp2 = temp .* ((1j*2*pi*ms*T) .* I + dI);
        temp = temp .* I;
        temp3 = temp .* (-1j*2*pi*deltaf*ns);

        g0(n+1, :) = sum(temp, 1);
        g2(n+1, :) = sum(temp2, 1);
        g3(n+1, :) = sum(temp3, 1);
    end
    g0 = g0(:);
    g2 = alpha * g2(:);
    g3 = alpha * g3(:);

    
    % g0 = g0;
    g1 = 1j * g0;
    % g2 = 1j * 2 * pi * T * alpha .* m_array .* g0; % Incorrect !!!!!!!!!
    % g3 = -1j * 2 * pi * deltaf * alpha.* n_array .* g0;
    
    G  = [g0, g1, g2, g3];             % NM x 4
    J  = (2/sigma2) * real(G' * G);    % exact J_like
    CRB = inv(J);
end