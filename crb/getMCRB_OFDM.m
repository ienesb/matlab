function MCRB = getMCRB_OFDM(X_hat, mu_gt, deltaf, T, sigma2, alpha_gt, nu_gt, tau_gt) % deltaf, T, X_hat, mu_gt, nu_gt, tau_gt

    options = optimoptions('fmincon', ...
                  'Algorithm','interior-point', ...
                  'Display','off', ...
                  'MaxFunctionEvaluations',1e5, ...
                  'FiniteDifferenceType','central', ...
                  'FiniteDifferenceStepSize',[1e-2; 1e-10], ...  % (Hz step, seconds step)
                  'StepTolerance',1e-14, ...
                  'OptimalityTolerance',1e-12);

    eta0 = [nu_gt; tau_gt];   % initial guess
        
    lb = [nu_gt-10; tau_gt-1e-7];
    ub = [nu_gt+10; tau_gt+1e-7];

    [eta_0, ~] = fmincon(@(eta) kl_distance(eta, deltaf, T, X_hat, mu_gt), ...
              eta0, [], [], [], [], lb, ub, [], options);

    nu_0  = eta_0(1);
    tau_0 = eta_0(2);

    [N, M] = size(X_hat);

    ns = (0:(N-1)).';
    ms = 0:(M-1);


    c_tilde = (exp(-1j*2*pi*deltaf*ns*tau_0) * exp(1j*2*pi*T*ms*nu_0)) .* X_hat;
    c_tilde = c_tilde(:);

    alpha_0 = pinv(c_tilde) * mu_gt;

    mu_tilde = alpha_0 * c_tilde;


    m = ones(N, 1) * ms;
    m = m(:);
    n = ns * ones(1, M);
    n = n(:);
    
    %% MCRB calculation
    g0 = c_tilde;
    g1 = 1j * c_tilde;
    g2 = 1j * 2 * pi * T * alpha_0 .* m .* c_tilde;
    g3 = -1j * 2 * pi * deltaf * alpha_0 .* n .* c_tilde;
    
    G  = [g0, g1, g2, g3];             % NM x 4
    J  = (2/sigma2) * real(G' * G);    % exact J_like

    % b vector

    Delta = mu_gt - mu_tilde;
    g_tilde_0 = c_tilde;
    g_tilde_1 = 1j * c_tilde;
    g_tilde_2 = 1j * 2 * pi * T * alpha_0 * m .* c_tilde;
    g_tilde_3 = -1j * 2 * pi * deltaf * alpha_0 * n .* c_tilde;
    b = zeros(4, 1);
    b(1) = Delta' * g_tilde_0;
    b(2) = Delta' * g_tilde_1;
    b(3) = Delta' * g_tilde_2;
    b(4) = Delta' * g_tilde_3;
    b = 2 / sigma2 * real(b);

    % A_delta matrix
    h_03 = -1j * 2 * pi * deltaf * n .* c_tilde;
    h_02 = 1j * 2 * pi * T * m .* c_tilde;
    h_13 = 2 * pi * deltaf * n .* c_tilde;
    h_12 = -2 * pi * T * m .* c_tilde;
    h_23 = alpha_0 * 4 * pi^2 * deltaf * T * m .* n .* c_tilde;
    h_22 = -alpha_0 * 4 * pi^2 * T^2 * m.^2 .* c_tilde;
    h_33 = -alpha_0 * 4 * pi^2 * deltaf^2 * n.^2 .* c_tilde;

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

    eta_gt = [real(alpha_gt), imag(alpha_gt), nu_gt, tau_gt].';
    eta_0 = [real(alpha_0), imag(alpha_0), nu_0, tau_0].';

    MCRB = MCRB + (eta_gt - eta_0) * (eta_gt - eta_0).';
end