% plot_mcrb_dd.m
% Plot BER and DD-domain MCRB results saved by eva_ber_mp_mcrb_dd.m.
% Loads eva_ber_mp_mcrb_dd_results.mat and produces:
%   - OTFS BER (CSI, CE)
%   - Per-path delay RMSE (CSI vs CE overlay)
%   - Per-path Doppler RMSE (CSI vs CE overlay)
%   - Per-path |alpha| RMSE (CSI vs CE overlay)
clear;

r_dd = load('eva_ber_mp_mcrb_dd_results.mat');
r_tf = load('eva_ber_mp_mcrb_results.mat');

SNR_dbs        = r_dd.SNR_dbs;
P              = r_dd.P;
SNR_p_dBs      = r_dd.SNR_p_dBs;
n_pilot        = r_dd.n_pilot;
bers_otfs_csi  = r_dd.bers_otfs_csi;
bers_otfs_ce   = r_dd.bers_otfs_ce;
MCRBs_OTFS_CSI_DD = r_dd.MCRBs_OTFS_CSI;
MCRBs_OTFS_CE_DD  = r_dd.MCRBs_OTFS_CE;

MCRBs_OTFS_CSI_TF = r_tf.MCRBs_OTFS;
MCRBs_OTFS_CE_TF  = r_tf.MCRBs_OTFS_CE;

% Parameter ordering per path: [alpha_R, alpha_I, nu, tau]
%   alpha_R: 4p-3   alpha_I: 4p-2   nu: 4p-1   tau: 4p
IDX_AR  = @(p) 4*p - 3;
IDX_AI  = @(p) 4*p - 2;
IDX_NU  = @(p) 4*p - 1;
IDX_TAU = @(p) 4*p;

% Reduce MCRB across MC dimension (median is more robust than mean for
% MCRB which can have heavy tails from poor pseudo-true fits).
MCRB_CSI_DD = mean(MCRBs_OTFS_CSI_DD, 4);                   % [4P x 4P x nSNR]
MCRB_CE_DD  = squeeze(mean(MCRBs_OTFS_CE_DD,  4));          % [4P x 4P x nSNR x n_pilot]
MCRB_CSI_TF = mean(MCRBs_OTFS_CSI_TF, 4);                   % [4P x 4P x nSNR]
MCRB_CE_TF  = squeeze(mean(MCRBs_OTFS_CE_TF,  4));          % [4P x 4P x nSNR x n_pilot]

% Average BER across MC
ber_csi = mean(bers_otfs_csi, 2);                       % [nSNR x 1]
ber_ce  = squeeze(mean(bers_otfs_ce,  2));              % [nSNR x n_pilot]

colors = lines(P);

%% ---- BER ----
figure;
semilogy(SNR_dbs, ber_csi, '-o', 'LineWidth', 2); hold on;
markers = {'--s', '--^', '--v', '--d'};
for p_idx = 1:n_pilot
    semilogy(SNR_dbs, ber_ce(:, p_idx), markers{min(p_idx, end)}, 'LineWidth', 2);
end
grid on;
xlabel('SNR [dB]'); ylabel('BER');
ce_labels = arrayfun(@(x) sprintf('OTFS (CE + MP, SNR_p=%d dB)', x), ...
                     SNR_p_dBs, 'UniformOutput', false);
legend(['OTFS (CSI + MP)', ce_labels], 'Location', 'southwest');
title('OTFS BER — DD-domain MCRB run');
theme(gcf, "light");

%% ---- Delay RMSE (per-path, CSI vs CE) ----
for p = 1:P
    figure;
    h_lines = gobjects(1, (1 + n_pilot)*2);
    rmse_csi_dd = sqrt(squeeze(MCRB_CSI_DD(IDX_TAU(p), IDX_TAU(p), :)));
    rmse_csi_tf = sqrt(squeeze(MCRB_CSI_TF(IDX_TAU(p), IDX_TAU(p), :)));
    semilogy(SNR_dbs, rmse_csi_tf, '-', ...
        'Color', colors(1, :), 'LineWidth', 2); hold on;
    semilogy(SNR_dbs, rmse_csi_dd, '--', ...
        'Color', colors(1, :), 'LineWidth', 2); hold on;
    for p_idx = 1:n_pilot
        rmse_ce_dd = sqrt(squeeze(MCRB_CE_DD(IDX_TAU(p), IDX_TAU(p), :, p_idx)));
        rmse_ce_tf = sqrt(squeeze(MCRB_CE_TF(IDX_TAU(p), IDX_TAU(p), :, p_idx)));
        semilogy(SNR_dbs, rmse_ce_tf, '-', ...
            'Color', colors(2, :), 'LineWidth', 2);
        semilogy(SNR_dbs, rmse_ce_dd, '--', ...
            'Color', colors(2, :), 'LineWidth', 2);
    end
    grid on;
    xlabel('SNR [dB]'); ylabel('Delay RMSE [s]');
    title(sprintf('OTFS Delay MCRB (DD) Path %d', p));
    path_labels = sprintf('Path %d', p);
    ce_path_labels = cell(1, P*n_pilot);
    for p_idx = 1:n_pilot
        ce_path_labels{(p-1)*n_pilot + p_idx} = ...
            sprintf('Path %d (CE %d dB)', p, SNR_p_dBs(p_idx));
    end
    legend(["CSI (TF)", "CSI (DD)", "CE (TF)", "CE (DD)"], 'Location', 'southwest', 'NumColumns', 2);
    theme(gcf, "light");
end


return;

%% ---- Doppler RMSE (per-path, CSI vs CE) ----
figure;
for p = 1:P
    rmse_csi_dd = sqrt(squeeze(MCRB_CSI_DD(IDX_NU(p), IDX_NU(p), :)));
    semilogy(SNR_dbs, rmse_csi_dd, '-', 'Color', colors(p, :), 'LineWidth', 2); hold on;
    for p_idx = 1:n_pilot
        rmse_ce_dd = sqrt(squeeze(MCRB_CE_DD(IDX_NU(p), IDX_NU(p), :, p_idx)));
        semilogy(SNR_dbs, rmse_ce_dd, '--', 'Color', colors(p, :), 'LineWidth', 2);
    end
end
grid on;
xlabel('SNR [dB]'); ylabel('Doppler RMSE [Hz]');
title('OTFS Doppler MCRB (DD): CSI (solid) vs CE (dashed)');

%% ---- |alpha| RMSE (per-path, CSI vs CE) ----
% Total amplitude variance ≈ var(alpha_R) + var(alpha_I) for small bias.
figure;
for p = 1:P
    var_a_csi = squeeze(MCRB_CSI_DD(IDX_AR(p), IDX_AR(p), :)) + ...
                squeeze(MCRB_CSI_DD(IDX_AI(p), IDX_AI(p), :));
    semilogy(SNR_dbs, sqrt(var_a_csi), '-', ...
        'Color', colors(p, :), 'LineWidth', 2); hold on;
    for p_idx = 1:n_pilot
        var_a_ce = squeeze(MCRB_CE_DD(IDX_AR(p), IDX_AR(p), :, p_idx)) + ...
                   squeeze(MCRB_CE_DD(IDX_AI(p), IDX_AI(p), :, p_idx));
        semilogy(SNR_dbs, sqrt(var_a_ce), '--', ...
            'Color', colors(p, :), 'LineWidth', 2);
    end
end
grid on;
xlabel('SNR [dB]'); ylabel('|\alpha| RMSE');
title('OTFS Path-Gain MCRB (DD): CSI (solid) vs CE (dashed)');
