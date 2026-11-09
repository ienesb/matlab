% plot_mcrb.m
% Load mcrb_results.mat and plot MC-averaged MCRB(tau) and MCRB(nu) vs BSC rho
% for all four cases: {OFDM,OTFS} x {single,multi}.
% Single-target case in Figure 1 (2 subplots), multi-target in Figure 2 (2 subplots).

clear; clc; close all;
S = load('mcrb_results.mat');
rho_vals = S.rho_vals;

% =================== Figure 1: single target ===================
figure('Name','MCRB single target');

subplot(1,2,1);
semilogy(rho_vals, S.results.ofdm_single.MCRB_tau_avg(:,1), '-o', 'LineWidth', 1.6, ...
    'DisplayName', 'OFDM (ICI)'); hold on;
semilogy(rho_vals, S.results.otfs_single.MCRB_tau_avg(:,1), '--s', 'LineWidth', 1.6, ...
    'DisplayName', 'OTFS');
grid on;
xlabel('\rho (BSC bit error probability)');
ylabel('MCRB(\tau) [s^2]');
title('MCRB(\tau) — single target (first EVA path)');
legend('Location','northwest');

subplot(1,2,2);
semilogy(rho_vals, S.results.ofdm_single.MCRB_nu_avg(:,1), '-o', 'LineWidth', 1.6, ...
    'DisplayName', 'OFDM (ICI)'); hold on;
semilogy(rho_vals, S.results.otfs_single.MCRB_nu_avg(:,1), '--s', 'LineWidth', 1.6, ...
    'DisplayName', 'OTFS');
grid on;
xlabel('\rho (BSC bit error probability)');
ylabel('MCRB(\nu) [Hz^2]');
title('MCRB(\nu) — single target (first EVA path)');
legend('Location','northwest');

% =================== Figure 2: multi target (reduced EVA, 3 paths) ===================
P = 3;
colors = lines(P);

figure('Name','MCRB multi target');

subplot(1,2,1);
for p = 1:P
    semilogy(rho_vals, S.results.ofdm_multi.MCRB_tau_avg(:,p), '-o', 'Color', colors(p,:), ...
        'LineWidth', 1.5, 'DisplayName', sprintf('Path %d, OFDM (ICI)', p)); hold on;
    semilogy(rho_vals, S.results.otfs_multi.MCRB_tau_avg(:,p), '--s', 'Color', colors(p,:), ...
        'LineWidth', 1.5, 'DisplayName', sprintf('Path %d, OTFS', p));
end
grid on;
xlabel('\rho (BSC bit error probability)');
ylabel('MCRB(\tau_p) [s^2]');
title('MCRB(\tau) — reduced EVA multi-target');
legend('Location','northwest','NumColumns',2);

subplot(1,2,2);
for p = 1:P
    semilogy(rho_vals, S.results.ofdm_multi.MCRB_nu_avg(:,p), '-o', 'Color', colors(p,:), ...
        'LineWidth', 1.5, 'DisplayName', sprintf('Path %d, OFDM (ICI)', p)); hold on;
    semilogy(rho_vals, S.results.otfs_multi.MCRB_nu_avg(:,p), '--s', 'Color', colors(p,:), ...
        'LineWidth', 1.5, 'DisplayName', sprintf('Path %d, OTFS', p));
end
grid on;
xlabel('\rho (BSC bit error probability)');
ylabel('MCRB(\nu_p) [Hz^2]');
title('MCRB(\nu) — reduced EVA multi-target');
legend('Location','northwest','NumColumns',2);
