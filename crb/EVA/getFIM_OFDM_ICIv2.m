function MCRB = getFIM_OFDM_ICIv2(X_hat, mu_gt, deltaf, Ts, sigma2, alpha_gt, nu_gt, tau_gt, P) % deltaf, T, X_hat, mu_gt, nu_gt, tau_gt
    T = 1 / deltaf;
    [N, M] = size(X_hat);

    ns = (0:(N-1)).'; % 1
    ms = 0:(M-1); % 2

    nps = 0:(N-1);
    nps = permute(nps, [1, 3, 2]); % 3

    X_permuted = permute(X_hat, [3 , 2, 1]);
    
    %% Pseudo-True Parameter Calculation
    % For now use GT parameters for PT parameters
    tau_0 = tau_gt;
    nu_0 = nu_gt;
    alpha_0 = alpha_gt;
    mu_tilde = channel(X_hat, alpha_0, tau_0, nu_0, deltaf, Ts, 0); % Assumed mean function evaluated at pseudo-true parameter.
    mu_tilde = mu_tilde(:);


    G = zeros(N*M, 4*P);
    H = zeros(N*M, 4*P, 4*P);

    for p = 1:P
        alpha = alpha_gt(p);
        nu = nu_gt(p);
        tau = tau_gt(p);
        I = get_ici_coeffs(nps-ns, nu, T);
        dI = get_ici_coeffs_derivative(nps-ns, nu, T);
        ddI = get_ici_coeffs_derivative(nps-ns, nu, T); % replace this with second derivative

        temp = (exp(-1j*2*pi*deltaf*nps*tau) .* exp(1j*2*pi*Ts*ms*nu));
        temp = X_permuted .* temp;
        temp2 = temp .* ((1j*2*pi*ms*Ts) .* I + dI);
        temp = temp .* I;
        temp3 = temp .* (-1j*2*pi*deltaf*nps);

        g0 = sum(temp, 3);
        g1 = 1j * g0;
        g2 = alpha * sum(temp2, 3);
        g3 = alpha * sum(temp3, 3);

        G(:, 4*p-3) = g0(:); % alphaR % Change the order here !!!!!!!!!!!!!!!!!!!
        G(:, 4*p-2) = g1(:); % alphaI
        G(:, 4*p-1) = g2(:); % nu
        G(:, 4*p) = g3(:); % tau 
    end

    
    % G: NM x 4P matrix
    % H: NM x 4P^2 matrix
    % Delta: NM x 1 vector
    
    J = 2 / sigma2 * real(G' * G);

    A = J(1:27, 1:27);
    B = J(1:27, 28:36);
    C = J(28:36, 1:27);
    D = J(28:36, 28:36);
    J_tau = D - C * inv(A) * B;

end