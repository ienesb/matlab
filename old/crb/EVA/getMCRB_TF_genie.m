function MCRB = getMCRB_TF_genie(X_true, deltaf, Ts, sigma2, alpha_gt, nu_gt, tau_gt, P)
    T = 1 / deltaf;
    [N, M] = size(X_true);

    ns = (0:(N-1)).'; % 1
    ms = 0:(M-1); % 2

    nps = 0:(N-1);
    nps = permute(nps, [1, 3, 2]); % 3

    X_permuted = permute(X_true, [3 , 2, 1]);

    G = zeros(N*M, 4*P);

    for p = 1:P
        alpha = alpha_gt(p);
        nu = nu_gt(p);
        tau = tau_gt(p);
        I = get_ici_coeffs(nps-ns, nu, T);
        dI = get_ici_coeffs_derivative(nps-ns, nu, T);

        temp = (exp(-1j*2*pi*deltaf*nps*tau) .* exp(1j*2*pi*Ts*ms*nu));
        temp = X_permuted .* temp;
        temp2 = temp .* ((1j*2*pi*ms*Ts) .* I + dI);
        temp = temp .* I;
        temp3 = temp .* (-1j*2*pi*deltaf*nps);

        g0 = sum(temp, 3);
        g1 = 1j * g0;
        g2 = alpha * sum(temp2, 3);
        g3 = alpha * sum(temp3, 3);

        G(:, 4*p-3) = g0(:);
        G(:, 4*p-2) = g1(:);
        G(:, 4*p-1) = g2(:);
        G(:, 4*p) = g3(:);
    end

    J = 2 / sigma2 * real(G' * G);

    MCRB = inv(J);
end
