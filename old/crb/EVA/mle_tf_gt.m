function [alpha_hat, tau_hat, nu_hat] = mle_tf_gt(Y_TF, X_hat, deltaf, Ts, tau_gt, nu_gt)
% MLE of multipath (alpha, tau, nu) parameters in the TF domain,
% initialized with ground-truth delays and Dopplers.
%
% Signal model (ICI neglected, valid for nu*T << 1):
%     g_{n,m}(tau,nu) = X_hat_{n,m} * exp(-j 2*pi*deltaf*n*tau)
%                                   * exp( j 2*pi*Ts*m*nu)
%
% Joint 2P-dimensional optimization over (tau_1..tau_P, nu_1..nu_P)
% via lsqnonlin + Levenberg-Marquardt, with alpha concentrated out
% (variable projection). Initial point = (tau_gt, nu_gt).

    [N, M] = size(Y_TF);
    ns = (0:N-1).';
    ms = 0:M-1;

    P = length(tau_gt);
    tau_gt = tau_gt(:);
    nu_gt  = nu_gt(:);

    y = Y_TF(:);
    theta0 = [tau_gt; nu_gt];

    residual_fn = @(theta) ml_residual(theta, X_hat, y, P, ...
                                       deltaf, Ts, ns, ms);

    opts_lm = optimoptions('lsqnonlin', ...
                           'Display', 'off', ...
                           'Algorithm', 'levenberg-marquardt', ...
                           'MaxIterations', 200, ...
                           'FunctionTolerance', 1e-14, ...
                           'StepTolerance', 1e-14);
    theta_hat = lsqnonlin(residual_fn, theta0, [], [], opts_lm);

    tau_hat = theta_hat(1:P);
    nu_hat  = theta_hat(P+1:2*P);

    G = build_basis(tau_hat, nu_hat, X_hat, deltaf, Ts, ns, ms);
    alpha_hat = G \ y;
end

function r = ml_residual(theta, X_hat, y, P, deltaf, Ts, ns, ms)
    tau = theta(1:P);
    nu  = theta(P+1:2*P);
    G   = build_basis(tau, nu, X_hat, deltaf, Ts, ns, ms);
    alpha     = G \ y;
    r_complex = G * alpha - y;
    r = [real(r_complex); imag(r_complex)];
end

function G = build_basis(tau, nu, X_hat, deltaf, Ts, ns, ms)
    [N, M] = size(X_hat);
    P = length(tau);
    G = zeros(N*M, P);
    for p = 1:P
        gp = X_hat .* (exp(-1j*2*pi*deltaf*ns*tau(p)) * ...
                       exp( 1j*2*pi*Ts*ms*nu(p)));
        G(:, p) = gp(:);
    end
end
