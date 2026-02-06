function MCRB = getMCRB_OTFS(X_DD, X_DD_hat, deltaf, T, sigma2, alpha, nu, tau) % monostatic [alpha_R, alpha_I, nu, tau]
    [N, M] = size(X_DD);

    ns = (0:(N-1)).';
    ms = 0:(M-1);

    X_TF_hat = isfft(X_DD_hat);
    
    C_TF = exp(-1j*2*pi*deltaf*ns*tau) * exp(1j*2*pi*T*ms*nu); % convert these to pseudo-true parameters
    C_tau_TF = (-1j*2*pi*deltaf*ns) .* C_TF;
    C_nu_TF = (1j*2*pi*T*ms) .* C_TF;
    C_tau2_TF = -(2*pi*deltaf*ns).^2 .* C_TF;
    C_nu2_TF = -(2*pi*T*ms).^2 .* C_TF;
    C_taunu_TF = (2*pi*T*ms) .* ((2*pi*deltaf*ns) .* C_TF);
    
    
    H_TF = alpha * C_TF;

    % C_DD = sfft(C_TF);
    % C_tau_DD = sfft(C_tau_TF);
    % C_nu_DD = sfft(C_nu_TF);
    % C_tau2_DD = sfft(C_tau2_TF);
    % C_nu2_DD = sfft(C_nu2_TF);
    % C_taunu_DD = sfft(C_taunu_TF);
    H_DD = sfft(H_TF);

    % g_tilde_0 = circular_conv2d(C_DD, X_DD_hat);
    % g_tilde_1 = 1j * circular_conv2d(C_DD, X_DD_hat);
    % g_tilde_2 = alpha * circular_conv2d(C_nu_DD, X_DD_hat);
    % g_tilde_3 = alpha * circular_conv2d(C_tau_DD, X_DD_hat);
    g_tilde_0 = sfft(C_TF .* X_TF_hat);
    g_tilde_1 = 1j * sfft(C_TF .* X_TF_hat);
    g_tilde_2 = alpha * sfft(C_nu_TF .* X_TF_hat);
    g_tilde_3 = alpha * sfft(C_tau_TF .* X_TF_hat);

    g_tilde_0 = g_tilde_0(:);
    g_tilde_1 = g_tilde_1(:);
    g_tilde_2 = g_tilde_2(:);
    g_tilde_3 = g_tilde_3(:);

    G = [g_tilde_0 g_tilde_1 g_tilde_2 g_tilde_3];
    J = 2/sigma2 * real(G' * G);

    mu_gt = circular_conv2d(H_DD, X_DD);
    mu_tilde = circular_conv2d(H_DD, X_DD_hat);

    mu_gt = mu_gt(:);
    mu_tilde = mu_tilde(:);

    Delta = mu_gt - mu_tilde;
    b = zeros(4, 1);
    b(1) = Delta' * g_tilde_0;
    b(2) = Delta' * g_tilde_1;
    b(3) = Delta' * g_tilde_2;
    b(4) = Delta' * g_tilde_3;
    b = 2 / sigma2 * real(b);

    % h_03 = circular_conv2d(C_tau_DD, X_DD_hat);
    % h_02 = circular_conv2d(C_nu_DD, X_DD_hat);
    h_03 = sfft(C_tau_TF .* X_TF_hat);
    h_02 = sfft(C_nu_TF .* X_TF_hat);
    h_03 = h_03(:);
    h_02 = h_02(:);
    h_13 = 1j * h_03;
    h_12 = 1j * h_02;

    % h_23 = alpha * circular_conv2d(C_taunu_DD, X_DD_hat);
    % h_22 = alpha * circular_conv2d(C_nu2_DD, X_DD_hat);
    % h_33 = alpha * circular_conv2d(C_tau2_DD, X_DD_hat);
    h_23 = alpha * sfft(C_taunu_TF .* X_TF_hat);
    h_22 = alpha * sfft(C_nu2_TF .* X_TF_hat);
    h_33 = alpha * sfft(C_tau2_TF .* X_TF_hat);
    h_23 = h_23(:);
    h_22 = h_22(:);
    h_33 = h_33(:);

    A_delta = zeros(4, 4);
    A_delta(1, 1) = 0;
    A_delta(1, 2) = 0;
    A_delta(2, 1) = 0;
    A_delta(2, 2) = 0;

    A_delta(1, 3) = Delta' * h_02;
    A_delta(1, 4) = Delta' * h_03;
    A_delta(2, 3) = Delta' * h_12;
    A_delta(2, 4) = Delta' * h_13;
    A_delta(3, 3) = Delta' * h_22;
    A_delta(3, 4) = Delta' * h_23;
    A_delta(4, 4) = Delta' * h_33;

    A_delta(3, 1) = A_delta(1, 3);
    A_delta(4, 1) = A_delta(1, 4);
    A_delta(3, 2) = A_delta(2, 3);
    A_delta(4, 2) = A_delta(2, 4); 
    A_delta(4, 3) = A_delta(3, 4);

    A_delta = 2 / sigma2 * real(A_delta);

    % B Matrix
    B = J + b * b.';

    % A matrix
    A = -J + A_delta;
    
    % MCRB
    % invA = inv(A); % !!!!!!!!!!!!
    MCRB = A \ B / A;

    % eta_gt = [real(alpha_gt), imag(alpha_gt), nu_gt, tau_gt].';
    % eta_0 = [real(alpha_0), imag(alpha_0), nu_0, tau_0].';
    % 
    % MCRB = MCRB + (eta_gt - eta_0) * (eta_gt - eta_0).';

end