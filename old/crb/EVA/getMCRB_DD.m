function [MCRB, A, B, bias_mat] = getMCRB_DD(x_hat_dd, x_true_dd, deltaf, Ts, sigma2, alpha_gt, nu_gt, tau_gt, P)
% getMCRB_DD - MCRB for OTFS channel parameters in the DD domain.
%
% Returns:
%   MCRB     [4P x 4P] full MCRB = A^{-1} B A^{-1} + (eta_gt-eta_0)(eta_gt-eta_0)'
%   A        [4P x 4P] = -J + A_Delta
%   B        [4P x 4P] = J + b b'
%   bias_mat [4P x 4P] = (eta_gt - eta_0) * (eta_gt - eta_0)'
% Variance-only MCRB (no bias term):  MCRB - bias_mat   (== A \ B / A)
%
% Implements the MCRB matrices A, B, J defined in otfs_io_derivatives.tex:
%   J = (2/sigma2) * Re{ G' G }
%   b = (2/sigma2) * Re{ G' Delta }
%   A_Delta = (2/sigma2) * Re{ H' Delta }     (per-element, then reshape)
%   A = -J + A_Delta
%   B =  J + b b'
%   MCRB = A^{-1} B A^{-1} + (eta_gt - eta_0)(eta_gt - eta_0)'
%
% with G, H built from the per-path partial derivatives of the DD mean
%
%   mu_p[k,l] = SFFT( TF_channel( ISFFT(x); alpha=1, tau_p, nu_p ) )[k,l]
%
% (last equation of otfs_io_equations.tex, factored as ISFFT -> TF ICI
% channel -> SFFT). Since SFFT is unitary, the inner products
% G' G, G' Delta, H' Delta are domain-invariant and the returned MCRB
% matches getMCRB_TF up to numerical error.

    T = 1 / deltaf;
    [N, M] = size(x_hat_dd);

    ns  = (0:(N-1)).';                      % dim 1: n  (output subcarrier)
    ms  = 0:(M-1);                          % dim 2: m
    nps = permute(0:(N-1), [1, 3, 2]);      % dim 3: n' (summed over)

    % --- DD inputs -> TF for using the TF channel kernel directly ---
    X_hat_TF  = isfft(x_hat_dd);
    X_true_TF = isfft(x_true_dd);
    X_hat_perm = permute(X_hat_TF, [3, 2, 1]);   % [1, M, N]

    % --- Ground-truth DD mean ---
    mu_gt = sfft(channel(X_true_TF, alpha_gt, tau_gt, nu_gt, deltaf, Ts, 0));
    mu_gt = mu_gt(:);

    % --- Pseudo-true parameter: argmin_eta || mu_tilde(eta; X_hat) - mu_gt ||^2 ---
    eta_gt = zeros(4*P, 1);
    for p = 1:P
        eta_gt(4*p-3) = real(alpha_gt(p));
        eta_gt(4*p-2) = imag(alpha_gt(p));
        eta_gt(4*p-1) = nu_gt(p);
        eta_gt(4*p)   = tau_gt(p);
    end
    residual_fn = @(eta) pt_residual_dd(eta, X_hat_TF, mu_gt, P, deltaf, Ts);
    options = optimoptions('lsqnonlin', 'Display', 'off', ...
                           'Algorithm', 'levenberg-marquardt');
    eta_0 = lsqnonlin(residual_fn, eta_gt, [], [], options);

    alpha_0 = zeros(P, 1); nu_0 = zeros(P, 1); tau_0 = zeros(P, 1);
    for p = 1:P
        alpha_0(p) = eta_0(4*p-3) + 1j*eta_0(4*p-2);
        nu_0(p)    = eta_0(4*p-1);
        tau_0(p)   = eta_0(4*p);
    end

    mu_tilde = sfft(channel(X_hat_TF, alpha_0, tau_0, nu_0, deltaf, Ts, 0));
    mu_tilde = mu_tilde(:);

    % --- Jacobian (G) and Hessian (H) blocks, evaluated at GT params ---
    G = zeros(N*M, 4*P);
    H = zeros(N*M, 4*P, 4*P);

    for p = 1:P
        alpha = alpha_gt(p);
        nu    = nu_gt(p);
        tau   = tau_gt(p);

        I   = get_ici_coeffs(nps - ns, nu, T);
        dI  = get_ici_coeffs_derivative(nps - ns, nu, T);
        ddI = get_ici_coeffs_derivative2(nps - ns, nu, T);

        phase  = exp(-1j*2*pi*deltaf*nps*tau) .* exp(1j*2*pi*Ts*ms*nu);
        base   = X_hat_perm .* phase;            % [1, M, N]: X̂ * phase

        kern_n = -1j*2*pi*deltaf*nps;            % d/dtau pulls down -j2π Δf n'
        a_nu   = base .* ((1j*2*pi*ms*Ts) .* I + dI);                                    % d/dnu kernel * X
        a_nn   = base .* (-(2*pi*ms*Ts).^2 .* I + (1j*4*pi*ms*Ts) .* dI + ddI);          % d^2/dnu^2 kernel * X
        b_I    = base .* I;
        a_tau  = b_I  .* kern_n;                                                         % d/dtau kernel * X

        % Reductions over n' (dim 3)  -> N x M arrays in TF domain
        g0_TF  = sum(b_I,    3);
        g2_TF  = alpha * sum(a_nu,  3);
        g3_TF  = alpha * sum(a_tau, 3);

        h02_TF =         sum(a_nu,  3);                          % d^2 / d a_R d nu  (no alpha)
        h03_TF =         sum(a_tau, 3);                          % d^2 / d a_R d tau (no alpha)
        h22_TF = alpha * sum(a_nn,  3);                          % d^2 / d nu^2
        h23_TF = alpha * sum(a_nu  .* kern_n, 3);                % d^2 / d tau d nu
        h33_TF = alpha * sum(a_tau .* kern_n, 3);                % d^2 / d tau^2

        % --- TF -> DD via SFFT (linear, applied per column) ---
        g0  = sfft(g0_TF);   g1  = 1j * g0;
        g2  = sfft(g2_TF);
        g3  = sfft(g3_TF);

        h02 = sfft(h02_TF);  h12 = 1j * h02;
        h03 = sfft(h03_TF);  h13 = 1j * h03;
        h22 = sfft(h22_TF);
        h23 = sfft(h23_TF);
        h33 = sfft(h33_TF);

        G(:, 4*p-3) = g0(:);    % alpha_R
        G(:, 4*p-2) = g1(:);    % alpha_I
        G(:, 4*p-1) = g2(:);    % nu
        G(:, 4*p)   = g3(:);    % tau

        H(:, 4*p-3, 4*p-1) = h02(:);
        H(:, 4*p-3, 4*p)   = h03(:);
        H(:, 4*p-2, 4*p-1) = h12(:);
        H(:, 4*p-2, 4*p)   = h13(:);
        H(:, 4*p-1, 4*p-3) = h02(:);
        H(:, 4*p-1, 4*p-2) = h12(:);
        H(:, 4*p-1, 4*p-1) = h22(:);
        H(:, 4*p-1, 4*p)   = h23(:);
        H(:, 4*p, 4*p-3)   = h03(:);
        H(:, 4*p, 4*p-2)   = h13(:);
        H(:, 4*p, 4*p-1)   = h23(:);
        H(:, 4*p, 4*p)     = h33(:);
    end
    H = reshape(H, [N*M, (4*P)^2]);

    Delta = mu_gt - mu_tilde;     % bar(mu) - tilde(mu)

    b       = 2/sigma2 * real(G' * Delta);
    A_delta = real(H' * Delta);
    A_delta = 2/sigma2 * reshape(A_delta, [4*P, 4*P]);
    J       = 2/sigma2 * real(G' * G);

    B = J + b * b.';
    A = -J + A_delta;

    bias_mat = (eta_gt - eta_0) * (eta_gt - eta_0).';
    MCRB     = A \ B / A + bias_mat;
end

function r = pt_residual_dd(eta, X_hat_TF, mu_gt_dd, P, deltaf, Ts)
    alpha = zeros(P, 1); nu = zeros(P, 1); tau = zeros(P, 1);
    for p = 1:P
        alpha(p) = eta(4*p-3) + 1j*eta(4*p-2);
        nu(p)    = eta(4*p-1);
        tau(p)   = eta(4*p);
    end
    mu_tilde_dd = sfft(channel(X_hat_TF, alpha, tau, nu, deltaf, Ts, 0));
    r_complex = mu_tilde_dd(:) - mu_gt_dd;
    r = [real(r_complex); imag(r_complex)];
end
