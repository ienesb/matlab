% plot_mmle.m
% --------------------------------------------------------------------------
% Empirical misspecified-MLE (MMLE) delay error on the reduced EVA channel
% (multi-target), trace over the delay (tau) block, for the four configs
% {OFDM, CP-OTFS} x {no-ICI, ICI}.  Reads compute_mmle_mc.m output (the
% per-(config,p,SNR) empirical MSE about the true parameter, MSE_bar):
%
%   tr{MSE_tau} = sum_p MSE_bar(tau_p, tau_p)    [s^2],
%   SNR := |alpha_1|^2 / sigma^2.
%
% Mirrors plot_mcrb.m, so the MMLE curves can be overlaid on the pseudo-true
% lower bound (both measured about eta_bar).  Two figures, four curves each:
%   Fig 1: fixed p = 0.1,    x-axis = SNR [dB]
%   Fig 2: fixed SNR = 0 dB, x-axis = p
% --------------------------------------------------------------------------
clear; clc; close all;
thisdir = fileparts(mfilename('fullpath'));
resdir  = fullfile(thisdir, 'results');

S = load(fullfile(resdir, 'compute_mmle_mc_results.mat'), ...
         'MSE_bar','configs','p_list','snr_db_list','idx_tau','P','nTrials');
nC = size(S.configs,1);  it = S.idx_tau;
np = numel(S.p_list);    ns = numel(S.snr_db_list);

% Scalar tau-trace of the empirical MSE, per (config, p, SNR).
mseT = zeros(nC, np, ns);          % tr{MSE_tau}  [s^2]
for ci = 1:nC
    for ip = 1:np
        for is = 1:ns
            mseT(ci,ip,is) = sum(real(diag(S.MSE_bar(it,it,ip,is,ci))));
        end
    end
end

styles = {'-o','-s','--^','--d'};
colors = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.47 0.67 0.19; 0.49 0.18 0.56];

% ===== Figure 1: tr{MSE_tau} vs SNR at fixed p = 0.1 =====================
p_fix = 0.1;  ip = find(abs(S.p_list - p_fix) < 1e-9, 1);
assert(~isempty(ip), 'p = %.2f not in p_list', p_fix);
figure('Name','MMLE MSE vs SNR (fixed p)','Position',[200 220 780 520]);
for ci = 1:nC
    semilogy(S.snr_db_list, squeeze(mseT(ci,ip,:)), styles{ci}, ...
             'Color', colors(ci,:), 'LineWidth', 1.8, 'MarkerSize', 6); hold on;
end
hold off; grid on;
xlabel('SNR = |\alpha_1|^2 / \sigma^2  [dB]');
ylabel('tr\{MSE_\tau\} = \Sigma_p MSE(\tau_p)  [s^2]');
title(sprintf('Empirical MMLE delay MSE vs SNR  (reduced EVA, P=%d, p=%.1f, %d trials)', ...
              S.P, p_fix, S.nTrials));
legend(S.configs(:,1), 'Location','southwest');  set(gca,'FontSize',11);
png1 = fullfile(resdir,'mmle_mse_vs_snr.png');  exportgraphics(gcf, png1, 'Resolution',150);

% ===== Figure 2: tr{MSE_tau} vs p at fixed SNR = 0 dB ===================
snr_fix = 0;  is = find(abs(S.snr_db_list - snr_fix) < 1e-9, 1);
assert(~isempty(is), 'SNR = %d dB not in snr_db_list', snr_fix);
figure('Name','MMLE MSE vs p (fixed SNR)','Position',[280 180 780 520]);
for ci = 1:nC
    semilogy(S.p_list, squeeze(mseT(ci,:,is)), styles{ci}, ...
             'Color', colors(ci,:), 'LineWidth', 1.8, 'MarkerSize', 6); hold on;
end
hold off; grid on;
xlabel('BSC error rate  p');
ylabel('tr\{MSE_\tau\} = \Sigma_p MSE(\tau_p)  [s^2]');
title(sprintf('Empirical MMLE delay MSE vs p  (reduced EVA, P=%d, SNR=%d dB, %d trials)', ...
              S.P, snr_fix, S.nTrials));
legend(S.configs(:,1), 'Location','northwest');  set(gca,'FontSize',11);
png2 = fullfile(resdir,'mmle_mse_vs_p.png');  exportgraphics(gcf, png2, 'Resolution',150);

fprintf('Saved %s\n       %s\n', png1, png2);
