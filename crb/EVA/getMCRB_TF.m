function MCRB = getMCRB_TF(X_hat, X_true, deltaf, Ts, sigma2, alpha_gt, nu_gt, tau_gt, P) % deltaf, T, X_hat, mu_gt, nu_gt, tau_gt
    T = 1 / deltaf;
    [N, M] = size(X_hat);

    ns = (0:(N-1)).'; % 1
    ms = 0:(M-1); % 2

    nps = 0:(N-1);
    nps = permute(nps, [1, 3, 2]); % 3

    X_permuted = permute(X_hat, [3 , 2, 1]);

    mu_gt = channel(X_true, alpha_gt, tau_gt, nu_gt, deltaf, Ts, 0);
    mu_gt = mu_gt(:);
    
    %% Pseudo-True Parameter Calculation
    % eta = [real(alpha), imag(alpha), nu, tau]
    % eta_0 = argmin_{eta} || mu_tilde(eta) - mu_gt ||^2
    % where mu_tilde is calculated with X_hat and mu_gt is calculated with
    % true X and true eta

    % Pack initial guess from GT parameters
    eta_gt = zeros(4*P, 1);
    for p = 1:P
        eta_gt(4*p-3) = real(alpha_gt(p));
        eta_gt(4*p-2) = imag(alpha_gt(p));
        eta_gt(4*p-1) = nu_gt(p);
        eta_gt(4*p) = tau_gt(p);
    end

    % Optimize: minimize ||mu_tilde(eta) - mu_gt||^2
    residual_fn = @(eta) pt_residual(eta, X_hat, mu_gt, P, deltaf, Ts);
    options = optimoptions('lsqnonlin', 'Display', 'off', 'Algorithm', 'levenberg-marquardt');
    eta_0 = lsqnonlin(residual_fn, eta_gt, [], [], options);

    % Unpack optimized pseudo-true parameters
    alpha_0 = zeros(P, 1);
    nu_0 = zeros(P, 1);
    tau_0 = zeros(P, 1);
    for p = 1:P
        alpha_0(p) = eta_0(4*p-3) + 1j*eta_0(4*p-2);
        nu_0(p) = eta_0(4*p-1);
        tau_0(p) = eta_0(4*p);
    end

    mu_tilde = channel(X_hat, alpha_0, tau_0, nu_0, deltaf, Ts, 0);
    mu_tilde = mu_tilde(:);


    G = zeros(N*M, 4*P);
    H = zeros(N*M, 4*P, 4*P);

    for p = 1:P
        alpha = alpha_gt(p);
        nu = nu_gt(p);
        tau = tau_gt(p);
        I = get_ici_coeffs(nps-ns, nu, T);
        dI = get_ici_coeffs_derivative(nps-ns, nu, T);
        ddI = get_ici_coeffs_derivative2(nps-ns, nu, T); % replace this with second derivative

        temp = (exp(-1j*2*pi*deltaf*nps*tau) .* exp(1j*2*pi*Ts*ms*nu));
        temp = X_permuted .* temp;
        temp2 = temp .* ((1j*2*pi*ms*Ts) .* I + dI);
        h22 = temp .* (-(2*pi*ms*Ts).^2 .* I + (1j*4*pi*ms*Ts) .* dI + ddI);
        temp = temp .* I;
        temp3 = temp .* (-1j*2*pi*deltaf*nps);

        g0 = sum(temp, 3);
        g1 = 1j * g0;
        g2 = alpha * sum(temp2, 3);
        g3 = alpha * sum(temp3, 3);

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

        h02 = sum(h02, 3);
        h03 = sum(h03, 3);
        h12 = sum(h12, 3);
        h13 = sum(h13, 3);
        h33 = sum(h33, 3);
        h23 = sum(h23, 3);
        h22 = sum(h22, 3);


        % g0 = zeros(N, M);
        % g2 = zeros(N, M);
        % g3 = zeros(N, M);
        % for n = 0:(N-1) % n -> n, ns -> np, ms -> m
        %     I = get_ici_coeffs(ns-n, nu, T);
        %     dI = get_ici_coeffs_derivative(ns-n, nu, T);
        %     ddI = get_ici_coeffs_derivative2(ns-n, nu, T);
        % 
        %     temp = (exp(-1j*2*pi*deltaf*ns*tau) * exp(1j*2*pi*Ts*ms*nu));
        %     temp = X .* temp;
        %     temp2 = temp .* ((1j*2*pi*ms*Ts) .* I + dI);
        %     h22 = temp .* (-(2*pi*ms*Ts).^2 .* I + (1j*4*pi*ms*Ts) .* dI + ddI);
        %     temp = temp .* I;
        %     temp3 = temp .* (-1j*2*pi*deltaf*ns);
        % 
        %     % h00 = 0; % 0: alpha_r, 1: alpha_i, 2: nu, 3: tau
        %     % h01 = 0;
        %     % h10 = 0;
        %     % h11 = 0;
        %     h02 = temp2;
        %     h03 = temp3;
        %     h12 = 1j * h02;
        %     h13 = 1j * h03;
        %     h33 = temp3 .* (-1j*2*pi*deltaf*ns);
        %     h23 = temp2 .* (-1j*2*pi*deltaf*ns);
        % 
        %     h02 = sum(h02, 1);
        %     h03 = sum(h03, 1);
        %     h12 = sum(h12, 1);
        %     h13 = sum(h13, 1);
        %     h33 = sum(h33, 1);
        %     h23 = sum(h23, 1);
        % 
        % 
        %     g0(n+1, :) = sum(temp, 1);
        %     g2(n+1, :) = sum(temp2, 1);
        %     g3(n+1, :) = sum(temp3, 1);
        % end
        % g0 = g0(:);
        % g1 = 1j * g0;
        % g2 = alpha * g2(:);
        % g3 = alpha * g3(:);

        G(:, 4*p-3) = g0(:); % alphaR % Change the order here !!!!!!!!!!!!!!!!!!!
        G(:, 4*p-2) = g1(:); % alphaI
        G(:, 4*p-1) = g2(:); % nu
        G(:, 4*p) = g3(:); % tau
        
        % H(4*p-3, 4*p-3) = h00;
        % H(4*p-3, 4*p-2) = h01;
        H(:, 4*p-3, 4*p-1) = h02(:);
        H(:, 4*p-3, 4*p) = h03(:);

        % H(4*p-2, 4*p-3) = h10;
        % H(4*p-2, 4*p-2) = h11;
        H(:, 4*p-2, 4*p-1) = h12(:);
        H(:, 4*p-2, 4*p) = h13(:);

        H(:, 4*p-1, 4*p-3) = h02(:);
        H(:, 4*p-1, 4*p-2) = h12(:);
        H(:, 4*p-1, 4*p-1) = h22(:);
        H(:, 4*p-1, 4*p) = h23(:);

        H(:, 4*p, 4*p-3) = h03(:);
        H(:, 4*p, 4*p-2) = h13(:);
        H(:, 4*p, 4*p-1) = h23(:);
        H(:, 4*p, 4*p) = h33(:);     
    end
    H = reshape(H, [N*M, (4*P)^2]);

    Delta = mu_gt - mu_tilde;
    
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

    MCRB = MCRB + (eta_gt - eta_0) * (eta_gt - eta_0).';
end

function r = pt_residual(eta, X_hat, mu_gt, P, deltaf, Ts)
    alpha = zeros(P, 1);
    nu = zeros(P, 1);
    tau = zeros(P, 1);
    for p = 1:P
        alpha(p) = eta(4*p-3) + 1j*eta(4*p-2);
        nu(p) = eta(4*p-1);
        tau(p) = eta(4*p);
    end
    mu_tilde = channel(X_hat, alpha, tau, nu, deltaf, Ts, 0);
    r_complex = mu_tilde(:) - mu_gt;
    r = [real(r_complex); imag(r_complex)];
end