function MCRB = getMCRB_OFDM(X_hat, mu_gt, deltaf, Ts, sigma2, alpha_gt, nu_gt, tau_gt, P) % deltaf, T, X_hat, mu_gt, nu_gt, tau_gt
    T = 1 / deltaf;
    [N, M] = size(X_hat);

    ns = (0:(N-1)).';
    ms = 0:(M-1);
    
    %% Pseudo-True Parameter Calculation
    % For now use GT parameters for PT parameters
    tau_0 = tau_gt;
    nu_0 = nu_gt;
    alpha_0 = alpha_gt;
    mu_tilde = channel(X_hat, alpha_0, tau_0, nu_0, deltaf, Ts, 0); % Assumed mean function evaluated at pseudo-true parameter.
    mu_tilde = mu_tilde(:);
    % options = optimoptions('fmincon', ...
    %               'Algorithm','interior-point', ...
    %               'Display','off', ...
    %               'MaxFunctionEvaluations',1e5, ...
    %               'FiniteDifferenceType','central', ...
    %               'FiniteDifferenceStepSize',[1e-2; 1e-10], ...  % (Hz step, seconds step)
    %               'StepTolerance',1e-14, ...
    %               'OptimalityTolerance',1e-12);
    % 
    % eta0 = [nu_gt; tau_gt];   % initial guess
    % 
    % lb = [nu_gt-10; tau_gt-1e-7];
    % ub = [nu_gt+10; tau_gt+1e-7];
    % 
    % [eta_0, ~] = fmincon(@(eta) kl_distance(eta, deltaf, T, X_hat, mu_gt), ...
    %           eta0, [], [], [], [], lb, ub, [], options);
    % 
    % nu_0  = eta_0(1);
    % tau_0 = eta_0(2);
    %
    % c_tilde = (exp(-1j*2*pi*deltaf*ns*tau_0) * exp(1j*2*pi*T*ms*nu_0)) .* X_hat;
    % c_tilde = c_tilde(:);
    % 
    % alpha_0 = pinv(c_tilde) * mu_gt;
    % 
    % mu_tilde = alpha_0 * c_tilde;

    
    % m = ones(N, 1) * ms;
    % m = m(:);
    % n = ns * ones(1, M);
    % n = n(:);

    G = zeros(N*M, 4*P);
    H = zeros(4*P, 4*P, N*M);

    for p = 1:P
        alpha = alpha_gt(p);
        nu = nu_gt(p);
        tau = tau_gt(p);

        g0 = zeros(N, M);
        g2 = zeros(N, M);
        g3 = zeros(N, M);
        for n = 0:(N-1) % n -> n, ns -> np, ms -> m
            I = get_ici_coeffs(ns-n, nu, T);
            dI = get_ici_coeffs_derivative(ns-n, nu, T);
            ddI = get_ici_coeffs_derivative2(ns-n, nu, T);
    
            temp = (exp(-1j*2*pi*deltaf*ns*tau) * exp(1j*2*pi*Ts*ms*nu));
            temp = X .* temp;
            temp2 = temp .* ((1j*2*pi*ms*Ts) .* I + dI);
            h22 = temp .* (-(2*pi*ms*Ts).^2 .* I + (1j*4*pi*ms*Ts) .* dI + ddI);
            temp = temp .* I;
            temp3 = temp .* (-1j*2*pi*deltaf*ns);

            % h00 = 0; % 0: alpha_r, 1: alpha_i, 2: nu, 3: tau
            % h01 = 0;
            % h10 = 0;
            % h11 = 0;
            h02 = temp2;
            h03 = temp3;
            h12 = 1j * h02;
            h13 = 1j * h03;
            h33 = temp3 .* (-1j*2*pi*deltaf*ns);
            h23 = temp2 .* (-1j*2*pi*deltaf*ns);

            h02 = sum(h02, 1);
            h03 = sum(h03, 1);
            h12 = sum(h12, 1);
            h13 = sum(h13, 1);
            h33 = sum(h33, 1);
            h23 = sum(h23, 1);
            
    
            g0(n+1, :) = sum(temp, 1);
            g2(n+1, :) = sum(temp2, 1);
            g3(n+1, :) = sum(temp3, 1);
        end
        g0 = g0(:);
        g1 = 1j * g0;
        g2 = alpha * g2(:);
        g3 = alpha * g3(:);

        G(:, 4*p-3) = g0; % alphaR % Change the order here !!!!!!!!!!!!!!!!!!!
        G(:, 4*p-2) = g1; % alphaI
        G(:, 4*p-1) = g2; % nu
        G(:, 4*p) = g3; % tau
        
        % H(4*p-3, 4*p-3) = h00;
        % H(4*p-3, 4*p-2) = h01;
        H(4*p-3, 4*p-1) = h02;
        H(4*p-3, 4*p) = h03;

        % H(4*p-2, 4*p-3) = h10;
        % H(4*p-2, 4*p-2) = h11;
        H(4*p-2, 4*p-1) = h12;
        H(4*p-2, 4*p) = h13;

        H(4*p-1, 4*p-3) = h02;
        H(4*p-1, 4*p-2) = h12;
        H(4*p-1, 4*p-1) = h22; % ------------------
        H(4*p-1, 4*p) = h23;

        H(4*p, 4*p-3) = h03;
        H(4*p, 4*p-2) = h13;
        H(4*p, 4*p-1) = h23;
        H(4*p, 4*p) = h33;
        
    end
    H = reshape(H, [(4*P)^2, N*M]);
    H = H.';
    % G: NM x 4P matrix
    % H: NM x 4P^2 matrix
    % Delta: NM x 1 vector

    b = 2 / sigma2 * real(G' * Delta);
    
    A_delta = real(H' * Delta);
    A_delta = 2 / sigma2 * reshape(A_delta, [4*P, 4*P]);
    
    J = 2 / sigma2 * real(G' * G);

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