% plot_l2_bsc.m
% Two diagnostics for the OFDM (ICI) reduced-EVA channel, both as a function
% of the BSC bit error probability rho:
%   (1) L2 error of the received mean   : ||mu_bar - mu_tilde||_2
%   (2) Bias of the pseudo-true parameter:
%       eta_0  = argmin_eta ||mu_bar(x; eta_bar) - mu_tilde(x_hat; eta)||_2^2
%       and the L2 norm ||eta_0 - eta_bar||_2 in scaled coordinates.
% mu follows the ICI form of eq (4) in ofdm_fim.tex.

clear; clc; close all;
addpath('../crb/EVA');   % for get_ici_coeffs

% ---- Grid / channel parameters ----
N      = 512;
M      = 128;
deltaf = 15e3;
T      = 1/deltaf;
Tcp    = 0;
Ts     = T + Tcp;

c        = 3e8;
fc       = 4e9;
velocity = 120 * 1e3 / 3600;
nu_max   = velocity * fc / c;

% ---- Reduced EVA channel (P = 3) ----
time_delays_ns  = [0,  370,  150];
relative_powers = [0,  -0.6, -1.4];
P     = numel(time_delays_ns);
taus  = time_delays_ns * 1e-9;
P_lin = 10.^(relative_powers/10);
P_lin = P_lin / sum(P_lin);
amps  = sqrt(P_lin);

rng(0);
phases = 2*pi*rand(P,1);
thetas = 2*pi*rand(P,1);
alphas = amps(:) .* exp(1j*phases);
nus    = nu_max * cos(thetas);

% ---- True parameter vector + per-parameter scaling ----
% Scaling: alpha unitless (s=1), tau in ns (s=1e9), nu in Hz (s=1).
% After scaling, all components live in a comparable range (~1 to ~400),
% which keeps lsqnonlin's finite-difference Jacobian well-conditioned.
eta_bar = zeros(4*P, 1);
for p = 1:P
    eta_bar(4*(p-1) + (1:4)) = [real(alphas(p)); imag(alphas(p)); taus(p); nus(p)];
end
scales         = repmat([1; 1; 1e9; 1], P, 1);
eta_bar_scaled = eta_bar .* scales;

box   = repmat([2; 2; 1000; 2*nu_max], P, 1);   % bound box (scaled units)
lb    = eta_bar_scaled - box;
ub    = eta_bar_scaled + box;

opts  = optimoptions('lsqnonlin','Display','off', ...
                     'Algorithm','trust-region-reflective', ...
                     'MaxIterations',60, 'MaxFunctionEvaluations',4000, ...
                     'FunctionTolerance',1e-12, 'StepTolerance',1e-12);

% ---- Sweep rho ----
rho_vals  = 0:0.05:0.5;
nMC       = 10;
l2_mu_avg  = zeros(size(rho_vals));
l2_eta_avg = zeros(size(rho_vals));

for ir = 1:numel(rho_vals)
    rho = rho_vals(ir);
    acc_mu  = 0;
    acc_eta = 0;
    for mc = 1:nMC
        bits = randi([0 3], N, M);
        x    = (1/sqrt(2)) * (2*bitand(bits,1) - 1 + 1j*(2*floor(bits/2) - 1));
        x_hat = bsc_qpsk(x, rho);

        mu_bar          = mu_ofdm_ici(x    , deltaf, Ts, alphas, taus, nus);
        mu_tilde_at_bar = mu_ofdm_ici(x_hat, deltaf, Ts, alphas, taus, nus);
        acc_mu = acc_mu + norm(mu_bar(:) - mu_tilde_at_bar(:));

        % Pseudo-true eta_0
        res = @(es) residual_ofdm(es, scales, mu_bar, x_hat, deltaf, Ts);
        eta0_scaled = lsqnonlin(res, eta_bar_scaled, lb, ub, opts);
        acc_eta = acc_eta + norm(eta0_scaled - eta_bar_scaled);
    end
    l2_mu_avg(ir)  = acc_mu  / nMC;
    l2_eta_avg(ir) = acc_eta / nMC;
    fprintf('rho = %.2f   ||mu_bar-mu_tilde||_2 = %.4e   ||eta_0-eta_bar||_2 (scaled) = %.4e\n', ...
            rho, l2_mu_avg(ir), l2_eta_avg(ir));
end

% ---- Plot ----
figure;
plot(rho_vals, l2_mu_avg, '-o', 'LineWidth', 1.6);
grid on;
xlabel('\rho (BSC bit error probability)');
ylabel('Average  ||\mu_{bar} - \mu_{tilde}||_2');
title('OFDM (ICI): \mu L2 error vs BSC \rho — reduced EVA, 3 paths');

figure;
plot(rho_vals, l2_eta_avg, '-o', 'LineWidth', 1.6);
grid on;
xlabel('\rho (BSC bit error probability)');
ylabel('Average  ||\eta_0 - \eta_{bar}||_2  (scaled units)');
title('OFDM (ICI): pseudo-true parameter bias vs BSC \rho — reduced EVA, 3 paths');

% ===================== local helpers =====================
function mu = mu_ofdm_ici(X, deltaf, Ts, alphas, taus, nus)
% Mean per eq (4) of ofdm_fim.tex, evaluated via per-path NxN Toeplitz I_mat.
    [N, M] = size(X);
    T = 1/deltaf;
    P = numel(alphas);

    ns = (0:N-1).';
    ms = (0:M-1);
    ds = (-(N-1):(N-1)).';

    mu = zeros(N, M);
    for p = 1:P
        A      = X .* exp(-1j*2*pi*deltaf*ns*taus(p)) .* exp(1j*2*pi*Ts*ms*nus(p));
        I_full = get_ici_coeffs(ds, nus(p), T);
        I_mat  = toeplitz(I_full(N:-1:1), I_full(N:end));
        mu     = mu + alphas(p) * (I_mat * A);
    end
end

function r = residual_ofdm(eta_s, scales, mu_bar, x_hat, deltaf, Ts)
% Returns the real-valued residual r = [Re(mu_bar-mu_tilde); Im(...)] of size 2NM.
    eta    = eta_s ./ scales;
    alphas = complex(eta(1:4:end), eta(2:4:end));
    taus   = eta(3:4:end);
    nus    = eta(4:4:end);
    mu     = mu_ofdm_ici(x_hat, deltaf, Ts, alphas, taus, nus);
    d      = mu_bar(:) - mu(:);
    r      = [real(d); imag(d)];
end
