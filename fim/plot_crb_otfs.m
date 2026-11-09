% plot_crb_otfs.m
% Plot CRLB(tau) and CRLB(nu) vs SNR for OTFS (single first-path target).
% OFDM no-ICI is overlaid as a reference -- otfs_fim.tex Section 3.1 proves
% the expected OTFS FIM equals the OFDM no-ICI FIM under i.i.d. QPSK.

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
nu    = nu_max;
alpha = 1;

% ---- SNR sweep ----
SNR_dbs  = 0:2:20;
SNR_lins = 10.^(SNR_dbs/10);

% ---- Monte Carlo over QPSK realizations ----
nMC = 20;
rng(0);

CRB_tau_OTFS_unit = 0;  CRB_nu_OTFS_unit = 0;
CRB_tau_OFDM_unit = 0;  CRB_nu_OFDM_unit = 0;
for mc = 1:nMC
    bits = randi([0 3], N, M);
    x    = (1/sqrt(2)) * (2*bitand(bits,1) - 1 + 1j*(2*floor(bits/2) - 1));

    J_otfs = getFIM_OTFS      (x, deltaf, Ts, 1, alpha, tau, nu);
    J_ofdm = getFIM_OFDM_noICI(x, deltaf, Ts, 1, alpha, tau, nu);

    C_otfs = invFIM_scaled(J_otfs);
    C_ofdm = invFIM_scaled(J_ofdm);

    CRB_tau_OTFS_unit = CRB_tau_OTFS_unit + C_otfs(3,3);
    CRB_nu_OTFS_unit  = CRB_nu_OTFS_unit  + C_otfs(4,4);
    CRB_tau_OFDM_unit = CRB_tau_OFDM_unit + C_ofdm(3,3);
    CRB_nu_OFDM_unit  = CRB_nu_OFDM_unit  + C_ofdm(4,4);
end
CRB_tau_OTFS_unit = CRB_tau_OTFS_unit / nMC;
CRB_nu_OTFS_unit  = CRB_nu_OTFS_unit  / nMC;
CRB_tau_OFDM_unit = CRB_tau_OFDM_unit / nMC;
CRB_nu_OFDM_unit  = CRB_nu_OFDM_unit  / nMC;

% Scale by 1/SNR_lin (FIM ∝ 1/sigma2)
CRB_tau_OTFS = CRB_tau_OTFS_unit ./ SNR_lins;
CRB_nu_OTFS  = CRB_nu_OTFS_unit  ./ SNR_lins;
CRB_tau_OFDM = CRB_tau_OFDM_unit ./ SNR_lins;
CRB_nu_OFDM  = CRB_nu_OFDM_unit  ./ SNR_lins;

% ---- Plot ----
figure;

subplot(1,2,1);
semilogy(SNR_dbs, CRB_tau_OTFS, '-o',  'LineWidth', 1.6, 'DisplayName', 'OTFS'); hold on;
semilogy(SNR_dbs, CRB_tau_OFDM, '--s', 'LineWidth', 1.6, 'DisplayName', 'OFDM no-ICI');
grid on;
xlabel('SNR [dB]');
ylabel('CRLB(\tau) [s^2]');
title('CRLB(\tau)');
legend('Location','northeast');

subplot(1,2,2);
semilogy(SNR_dbs, CRB_nu_OTFS, '-o',  'LineWidth', 1.6, 'DisplayName', 'OTFS'); hold on;
semilogy(SNR_dbs, CRB_nu_OFDM, '--s', 'LineWidth', 1.6, 'DisplayName', 'OFDM no-ICI');
grid on;
xlabel('SNR [dB]');
ylabel('CRLB(\nu) [Hz^2]');
title('CRLB(\nu)');
legend('Location','northeast');

sgtitle(sprintf('OTFS single-path CRLB, first EVA path, \\nu = \\nu_{max} = %.1f Hz', nu));

fprintf('Ratio CRB_OTFS/CRB_OFDM @ tau: %.6f\n', CRB_tau_OTFS_unit / CRB_tau_OFDM_unit);
fprintf('Ratio CRB_OTFS/CRB_OFDM @ nu : %.6f\n', CRB_nu_OTFS_unit  / CRB_nu_OFDM_unit );

% ---------- local helper ----------
function C = invFIM_scaled(J)
    d  = sqrt(diag(J));
    Js = J ./ (d * d.');
    Cs = inv(Js);
    C  = Cs ./ (d * d.');
end
