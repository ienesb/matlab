% plot_crb_multi.m
% CRLB(tau) and CRLB(nu) vs SNR for the reduced EVA channel (3 paths),
% comparing no-ICI vs ICI for each path. Two figures, 6 curves per figure.
% Per-path operating point (alpha, tau, nu) is fixed; QPSK data is averaged.

clear; clc; close all;
addpath('../crb/EVA');

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

% ---- SNR sweep ----
SNR_dbs  = 0:2:20;
SNR_lins = 10.^(SNR_dbs/10);

% ---- Monte Carlo over QPSK data (operating point fixed) ----
nMC = 20;

CRB_tau_noICI_unit = zeros(P,1);   CRB_nu_noICI_unit = zeros(P,1);
CRB_tau_ICI_unit   = zeros(P,1);   CRB_nu_ICI_unit   = zeros(P,1);

for mc = 1:nMC
    bits = randi([0 3], N, M);
    X    = (1/sqrt(2)) * (2*bitand(bits,1) - 1 + 1j*(2*floor(bits/2) - 1));

    J_no = getFIM_OFDM_multi_noICI(X, deltaf, Ts, 1, alphas, taus, nus);
    J_ic = getFIM_OFDM_multi_ICI  (X, deltaf, Ts, 1, alphas, taus, nus);

    C_no = invFIM_scaled(J_no);
    C_ic = invFIM_scaled(J_ic);

    for p = 1:P
        i_tau = 4*(p-1)+3;
        i_nu  = 4*(p-1)+4;
        CRB_tau_noICI_unit(p) = CRB_tau_noICI_unit(p) + C_no(i_tau,i_tau);
        CRB_nu_noICI_unit (p) = CRB_nu_noICI_unit (p) + C_no(i_nu, i_nu );
        CRB_tau_ICI_unit  (p) = CRB_tau_ICI_unit  (p) + C_ic(i_tau,i_tau);
        CRB_nu_ICI_unit   (p) = CRB_nu_ICI_unit   (p) + C_ic(i_nu, i_nu );
    end
end
CRB_tau_noICI_unit = CRB_tau_noICI_unit / nMC;
CRB_nu_noICI_unit  = CRB_nu_noICI_unit  / nMC;
CRB_tau_ICI_unit   = CRB_tau_ICI_unit   / nMC;
CRB_nu_ICI_unit    = CRB_nu_ICI_unit    / nMC;

% Scale by 1/SNR_lin
CRB_tau_noICI = CRB_tau_noICI_unit ./ SNR_lins;   % P x length(SNR_dbs)
CRB_tau_ICI   = CRB_tau_ICI_unit   ./ SNR_lins;
CRB_nu_noICI  = CRB_nu_noICI_unit  ./ SNR_lins;
CRB_nu_ICI    = CRB_nu_ICI_unit    ./ SNR_lins;

% ---- Plot ----
colors = lines(P);

% Figure 1: CRLB(tau)
figure;
for p = 1:P
    semilogy(SNR_dbs, CRB_tau_noICI(p,:), '-o', 'Color', colors(p,:), ...
        'LineWidth', 1.5, 'DisplayName', sprintf('Path %d, no ICI', p)); hold on;
    semilogy(SNR_dbs, CRB_tau_ICI  (p,:), '--s', 'Color', colors(p,:), ...
        'LineWidth', 1.5, 'DisplayName', sprintf('Path %d, ICI',    p));
end
grid on;
xlabel('SNR [dB]');
ylabel('CRLB(\tau_p) [s^2]');
title('OFDM CRLB(\tau) — reduced EVA, 3 paths');
legend('Location','northeast');

% Figure 2: CRLB(nu)
figure;
for p = 1:P
    semilogy(SNR_dbs, CRB_nu_noICI(p,:), '-o', 'Color', colors(p,:), ...
        'LineWidth', 1.5, 'DisplayName', sprintf('Path %d, no ICI', p)); hold on;
    semilogy(SNR_dbs, CRB_nu_ICI  (p,:), '--s', 'Color', colors(p,:), ...
        'LineWidth', 1.5, 'DisplayName', sprintf('Path %d, ICI',    p));
end
grid on;
xlabel('SNR [dB]');
ylabel('CRLB(\nu_p) [Hz^2]');
title('OFDM CRLB(\nu) — reduced EVA, 3 paths');
legend('Location','northeast');

% ---------- local helper ----------
function C = invFIM_scaled(J)
    d  = sqrt(diag(J));
    Js = J ./ (d * d.');
    Cs = inv(Js);
    C  = Cs ./ (d * d.');
end
