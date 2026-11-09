% plot_mle_vs_mcrb.m
% Per-path delay MLE RMSE overlaid with the corresponding MCRB curves.
% One figure per path, separately for OFDM and OTFS -> 2*P figures total.
clear;

mc  = load('eva_ber_mp_mcrb_results.mat');
mle = load('eva_ber_mp_mle_results.mat');

SNR_dbs = mc.SNR_dbs;
P       = mc.P;
nSNR    = length(SNR_dbs);

% ---- MLE RMSE (nSNR x P) ----
rmse_tau_ofdm_csi   = sqrt(squeeze(mean(mle.mse_tau_csi,        2)));
rmse_tau_ofdm_ce    = sqrt(squeeze(mean(mle.mse_tau_ce,         2)));
rmse_tau_ofdm_genie = sqrt(squeeze(mean(mle.mse_tau_genie,      2)));

rmse_tau_otfs_csi   = sqrt(squeeze(mean(mle.mse_tau_otfs_csi,   2)));
rmse_tau_otfs_genie = sqrt(squeeze(mean(mle.mse_tau_otfs_genie, 2)));
rmse_tau_otfs_ce    = sqrt(squeeze(mean(mle.mse_tau_otfs_ce,    2)));
if ndims(rmse_tau_otfs_ce) == 3
    rmse_tau_otfs_ce = rmse_tau_otfs_ce(:, :, 1);   % first pilot SNR
end

% ---- MCRB delay RMSE (nSNR x P) ----

MCRB_OFDM       = mean(mc.MCRBs_OFDM,    4);
MCRB_OFDM_CE    = mean(mc.MCRBs_OFDM_CE,    4);
MCRB_OFDM_genie = mean(mc.MCRBs_OFDM_genie,    4);
MCRB_OTFS       = mean(mc.MCRBs_OTFS,    4);
MCRB_OTFS_CE    = mean(mc.MCRBs_OTFS_CE,    4);
MCRB_OTFS_genie = mean(mc.MCRBs_OTFS_genie,    4);

crb_tau_ofdm_csi   = zeros(nSNR, P);
crb_tau_ofdm_ce    = zeros(nSNR, P);
crb_tau_ofdm_genie = zeros(nSNR, P);
crb_tau_otfs_csi   = zeros(nSNR, P);
crb_tau_otfs_ce    = zeros(nSNR, P);
crb_tau_otfs_genie = zeros(nSNR, P);
for p = 1:P
    crb_tau_ofdm_csi  (:, p) = sqrt(squeeze(MCRB_OFDM      (4*p, 4*p, :)));
    crb_tau_ofdm_ce   (:, p) = sqrt(squeeze(MCRB_OFDM_CE   (4*p, 4*p, :)));
    crb_tau_ofdm_genie(:, p) = sqrt(squeeze(MCRB_OFDM_genie(4*p, 4*p, :)));
    crb_tau_otfs_csi  (:, p) = sqrt(squeeze(MCRB_OTFS      (4*p, 4*p, :)));
    crb_tau_otfs_ce   (:, p) = sqrt(squeeze(MCRB_OTFS_CE   (4*p, 4*p, :, 1)));
    crb_tau_otfs_genie(:, p) = sqrt(squeeze(MCRB_OTFS_genie(4*p, 4*p, :)));
end

legend_entries = {'MLE CSI', 'MLE CE', 'MLE genie', ...
                  'MCRB CSI', 'MCRB CE', 'CRB genie'};

colors = lines(P);

%% ---- OFDM per-path figures ----
for p = 1:P
    figure;
    semilogy(SNR_dbs, rmse_tau_ofdm_csi(:,p),   '-o',  'LineWidth', 1.8, "Color", colors(1, :)); hold on;
    semilogy(SNR_dbs, rmse_tau_ofdm_ce(:,p),    '-s',  'LineWidth', 1.8, "Color", colors(2, :));
    semilogy(SNR_dbs, rmse_tau_ofdm_genie(:,p), '-^',  'LineWidth', 1.8, "Color", colors(3, :));
    semilogy(SNR_dbs, crb_tau_ofdm_csi(:,p),    '--',  'LineWidth', 1.4, "Color", colors(1, :));
    semilogy(SNR_dbs, crb_tau_ofdm_ce(:,p),     ':',   'LineWidth', 1.4, "Color", colors(2, :));
    semilogy(SNR_dbs, crb_tau_ofdm_genie(:,p),  '-.',  'LineWidth', 1.4, "Color", colors(3, :));
    grid on;
    xlabel('SNR [dB]');
    ylabel('Delay RMSE [s]');
    title(sprintf('OFDM — Path %d: Delay MLE vs MCRB', p));
    legend(legend_entries, 'Location', 'southwest');
    theme(gcf, "light");
    saveas(gcf, sprintf('mcrb_ofdm_path%d.png', p));
end

%% ---- OTFS per-path figures ----
for p = 1:P
    figure;
    semilogy(SNR_dbs, rmse_tau_otfs_csi(:,p),   '-o',  'LineWidth', 1.8, "Color", colors(1, :)); hold on;
    semilogy(SNR_dbs, rmse_tau_otfs_ce(:,p),    '-s',  'LineWidth', 1.8, "Color", colors(2, :));
    semilogy(SNR_dbs, rmse_tau_otfs_genie(:,p), '-^',  'LineWidth', 1.8, "Color", colors(3, :));
    semilogy(SNR_dbs, crb_tau_otfs_csi(:,p),    '--',  'LineWidth', 1.4, "Color", colors(1, :));
    semilogy(SNR_dbs, crb_tau_otfs_ce(:,p),     ':',   'LineWidth', 1.4, "Color", colors(2, :));
    semilogy(SNR_dbs, crb_tau_otfs_genie(:,p),  '-.',  'LineWidth', 1.4, "Color", colors(3, :));
    grid on;
    xlabel('SNR [dB]');
    ylabel('Delay RMSE [s]');
    title(sprintf('OTFS — Path %d: Delay MLE vs MCRB', p));
    legend(legend_entries, 'Location', 'southwest');
    theme(gcf, "light");
    saveas(gcf, sprintf('mcrb_otfs_path%d.png', p));
end
