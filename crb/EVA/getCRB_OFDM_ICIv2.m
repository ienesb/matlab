function CRB = getCRB_OFDM_ICIv2(X, deltaf, Ts, sigma2, alphas, nus, taus, mask) % only for delays, monostatic [alpha_R, alpha_I, nu, tau]
    [N, M] = size(X);

    T = 1 / deltaf;

    X(~mask) = 0;

    ns = (0:(N-1)).';
    ms = 0:(M-1);

    % m_array = ones(N, 1) * ms;
    % m_array = m_array(:);
    % n_array = ns * ones(1, M);
    % n_array = n_array(:);

    P = length(alphas);
    % CRB = zeros(4*P, 4*P);
    G = zeros(N*M, 4*P);
    for p = 1:P
        alpha = alphas(p);
        nu = nus(p);
        tau = taus(p);

        % c = (exp(-1j*2*pi*deltaf*ns*tau) * exp(1j*2*pi*T*ms*nu)) .* X;
        % c = channel_noiseless(X, 1, tau, nu, deltaf, T);
        g0 = zeros(N, M);
        g2 = zeros(N, M);
        g3 = zeros(N, M);
        for n = 0:(N-1)
            I = get_ici_coeffs(ns-n, nu, T);
            dI = get_ici_coeffs_derivative(ns-n, nu, T);
    
            temp = (exp(-1j*2*pi*deltaf*ns*tau) * exp(1j*2*pi*Ts*ms*nu));
            temp = X .* temp;
            temp2 = temp .* ((1j*2*pi*ms*Ts) .* I + dI);
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
        
        % G  = [g0, g1, g2, g3];             % NM x 4
        % J  = (2/sigma2) * real(G' * G);    % exact J_like
        % CRB((4*p-3):(4*p), (4*p-3):(4*p)) = inv(J);
        G(:, 4*p-3) = g0; % alphaR
        G(:, 4*p-2) = g1; % alphaI
        G(:, 4*p-1) = g2; % nu
        G(:, 4*p) = g3; % tau
    end

    G = G(:, [(1:P)*4-3, (1:P)*4-2, (1:P)*4-1, (1:P)*4]);

    J  = (2/sigma2) * real(G' * G);    % exact J_like

    A = J(1:27, 1:27);
    B = J(1:27, 28:36);
    C = J(28:36, 1:27);
    D = J(28:36, 28:36);
    J_tau = D - C * inv(A) * B;
    
    CRB = inv(J_tau);
end