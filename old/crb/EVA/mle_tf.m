function [alpha_hat, tau_hat, nu_hat] = mle_tf(Y_TF, X_hat, deltaf, Ts, P)
% MLE of multipath (alpha, tau, nu) parameters in the TF domain.
%
% Signal model (ICI neglected, valid for nu*T << 1):
%     g_{n,m}(tau,nu) = X_hat_{n,m} * exp(-j 2*pi*deltaf*n*tau)
%                                   * exp( j 2*pi*Ts*m*nu)
%
% Stage 1 (init, RELAX/SIC):
%   For each path: coarse 2D-FFT matched filter -> fminsearch refine ->
%   closed-form alpha -> residual update.
%
% Stage 2 (joint refinement):
%   Joint 2P-dimensional optimization over (tau_1..tau_P, nu_1..nu_P)
%   via lsqnonlin + Levenberg-Marquardt, with alpha concentrated out
%   (variable projection). Mirrors the pseudo-true parameter solver
%   used in getMCRB_TF.m / pt_residual.

    [N, M] = size(Y_TF);
    ns = (0:N-1).';
    ms = 0:M-1;

    %% ---- Stage 1: RELAX / SIC initialization ----
    osr_n = 4;
    osr_m = 4;
    Nf = N * osr_n;
    Mf = M * osr_m;

    r = Y_TF;
    tau_init = zeros(P, 1);
    nu_init  = zeros(P, 1);

    opts_nm = optimset('TolX', 1e-14, 'TolFun', 1e-14, ...
                       'MaxFunEvals', 400, 'Display', 'off');

    for p = 1:P
        A = r .* conj(X_hat);

        Ap = zeros(Nf, Mf);
        Ap(1:N, 1:M) = A;
        MF = Nf * ifft(Ap, Nf, 1);
        MF = fft(MF, Mf, 2);

        [~, idx] = max(abs(MF(:)));
        [li, ki] = ind2sub([Nf, Mf], idx);

        tau0 = (li - 1) / (Nf * deltaf);
        if (ki - 1) <= Mf/2
            nu0 = (ki - 1) / (Mf * Ts);
        else
            nu0 = (ki - 1 - Mf) / (Mf * Ts);
        end

        negll = @(x) -abs( sum(sum( A .* ( ...
            exp(1j*2*pi*deltaf*ns*x(1)) * ...
            exp(-1j*2*pi*Ts*ms*x(2)) ) )) )^2;
        x_opt = fminsearch(negll, [tau0; nu0], opts_nm);
        tau_p = x_opt(1);
        nu_p  = x_opt(2);

        g = X_hat .* (exp(-1j*2*pi*deltaf*ns*tau_p) * ...
                      exp( 1j*2*pi*Ts*ms*nu_p));
        alpha_p = (g(:)' * r(:)) / (g(:)' * g(:));

        tau_init(p) = tau_p;
        nu_init(p)  = nu_p;
        r = r - alpha_p * g;
    end

    %% ---- Stage 2: joint 2P-dim refinement (variable projection) ----
    y = Y_TF(:);
    theta0 = [tau_init; nu_init];

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
