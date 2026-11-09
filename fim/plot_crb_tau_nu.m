% plot_crb_tau_nu.m
% Plot CRLB(tau) and CRLB(nu) vs SNR for OFDM (no-ICI vs ICI), single first-path target.
% SNR is defined per data symbol with |alpha|^2 = 1 and |X[n,m]|^2 = 1, so
% sigma2 = 1/SNR_lin.

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

% First EVA path
tau   = 0;
nu    = nu_max;     % worst-case Doppler for the first path
alpha = 1;

% ---- SNR sweep ----
SNR_dbs  = 0:2:20;
SNR_lins = 10.^(SNR_dbs/10);

% ---- Monte Carlo over QPSK data (FIM is data-dependent in ICI case) ----
nMC = 20;
rng(0);

% Compute FIM at sigma2 = 1; CRB scales as 1/SNR_lin since J ∝ 1/sigma2.
CRB_tau_noICI_unit = 0;  CRB_nu_noICI_unit = 0;
CRB_tau_ICI_unit   = 0;  CRB_nu_ICI_unit   = 0;
for mc = 1:nMC
    bits = randi([0 3], N, M);
    X    = (1/sqrt(2)) * (2*bitand(bits,1) - 1 + 1j*(2*floor(bits/2) - 1));

    J_no = getFIM_OFDM_noICI(X, deltaf, Ts, 1, alpha, tau, nu);
    J_ic = getFIM_OFDM_ICI (X, deltaf, Ts, 1, alpha, tau, nu);

    % FIM diagonal spans ~15 decades (alpha ~ O(1), tau ~ O(1e-20)),
    % so a direct inv() is numerically singular. Symmetric diagonal
    % preconditioning: D = diag(sqrt(diag(J))), invert J./(D*D'), unscale.
    C_no = invFIM_scaled(J_no);
    C_ic = invFIM_scaled(J_ic);

    CRB_tau_noICI_unit = CRB_tau_noICI_unit + C_no(3,3);
    CRB_nu_noICI_unit  = CRB_nu_noICI_unit  + C_no(4,4);
    CRB_tau_ICI_unit   = CRB_tau_ICI_unit   + C_ic(3,3);
    CRB_nu_ICI_unit    = CRB_nu_ICI_unit    + C_ic(4,4);
end
CRB_tau_noICI_unit = CRB_tau_noICI_unit / nMC;
CRB_nu_noICI_unit  = CRB_nu_noICI_unit  / nMC;
CRB_tau_ICI_unit   = CRB_tau_ICI_unit   / nMC;
CRB_nu_ICI_unit    = CRB_nu_ICI_unit    / nMC;

% Scale by 1/SNR_lin
CRB_tau_noICI = CRB_tau_noICI_unit ./ SNR_lins;
CRB_tau_ICI   = CRB_tau_ICI_unit   ./ SNR_lins;
CRB_nu_noICI  = CRB_nu_noICI_unit  ./ SNR_lins;
CRB_nu_ICI    = CRB_nu_ICI_unit    ./ SNR_lins;

% ---- Plot ----
figure;

subplot(1,2,1);
semilogy(SNR_dbs, CRB_tau_noICI, '-o', 'LineWidth', 1.6, 'DisplayName', 'No ICI'); hold on;
semilogy(SNR_dbs, CRB_tau_ICI,   '-s', 'LineWidth', 1.6, 'DisplayName', 'ICI');
grid on;
xlabel('SNR [dB]');
ylabel('CRLB(\tau) [s^2]');
title('CRLB(\tau)');
legend('Location','northeast');

subplot(1,2,2);
semilogy(SNR_dbs, CRB_nu_noICI, '-o', 'LineWidth', 1.6, 'DisplayName', 'No ICI'); hold on;
semilogy(SNR_dbs, CRB_nu_ICI,   '-s', 'LineWidth', 1.6, 'DisplayName', 'ICI');
grid on;
xlabel('SNR [dB]');
ylabel('CRLB(\nu) [Hz^2]');
title('CRLB(\nu)');
legend('Location','northeast');

sgtitle(sprintf('OFDM single-path CRLB, \\nu = \\nu_{max} = %.1f Hz', nu));

fprintf('Ratio CRB_ICI/CRB_noICI @ tau: %.6f\n', CRB_tau_ICI_unit / CRB_tau_noICI_unit);
fprintf('Ratio CRB_ICI/CRB_noICI @ nu : %.6f\n', CRB_nu_ICI_unit  / CRB_nu_noICI_unit );

function C = invFIM_scaled(J)
% Invert a poorly-scaled symmetric positive-definite FIM by symmetric
% diagonal preconditioning. Returns the CRB matrix in original units.
    d = sqrt(diag(J));
    Js = J ./ (d * d.');
    Cs = inv(Js);
    C  = Cs ./ (d * d.');
end

