% plot_merged_mcrb.m
% Load and merge results from two MCRB runs, then plot.
clear;

r1 = load('eva_ber_mp_mcrb_results1.mat');
r2 = load('eva_ber_mp_mcrb_results2.mat');

SNR_dbs   = r1.SNR_dbs;
P         = r1.P;
SNR_p_dBs = r1.SNR_p_dBs;
n_pilot   = r1.n_pilot;

% Merge BERs (both are already averaged over 5 MC trials — take mean)
bers_ofdm_csi = (r1.bers_ofdm_csi + r2.bers_ofdm_csi) / 2;
bers_otfs_csi = (r1.bers_otfs_csi + r2.bers_otfs_csi) / 2;
bers_otfs_ce  = (r1.bers_otfs_ce  + r2.bers_otfs_ce)  / 2;

% Merge MCRBs (concatenate along MC dimension, then take median)
MCRBs_OFDM    = cat(4, r1.MCRBs_OFDM,    r2.MCRBs_OFDM);
MCRBs_OTFS    = cat(4, r1.MCRBs_OTFS,    r2.MCRBs_OTFS);
MCRBs_OTFS_CE = cat(4, r1.MCRBs_OTFS_CE, r2.MCRBs_OTFS_CE);

%% ---- BER Plot ----
figure;
semilogy(SNR_dbs, bers_ofdm_csi, '-d',  'LineWidth', 2); hold on;
semilogy(SNR_dbs, bers_otfs_csi, '-o',  'LineWidth', 2);
markers = {'--s', '--^', '--v'};
for p_idx = 1:n_pilot
    semilogy(SNR_dbs, bers_otfs_ce(:,:,p_idx), markers{p_idx}, 'LineWidth', 2);
end
grid on;
xlabel('SNR_d [dB]');
ylabel('BER');
ce_labels = arrayfun(@(x) sprintf('OTFS (CE + MP, SNR_p=%d dB)', x), SNR_p_dBs, 'UniformOutput', false);
legend(['OFDM (CSI)', 'OTFS (CSI + MP)', ce_labels], 'Location', 'southwest');
title('EVA Channel — Raviteja 2019 Fig. 9');
theme(gcf, "light");

%% ---- MCRB Plots ----
MCRB_OFDM = median(MCRBs_OFDM, 4);
colors = lines(P);

figure;
for p = [1, 5, 3]
    MCRB_p = squeeze(MCRB_OFDM(4*p, 4*p, :));
    MCRB_p = sqrt(MCRB_p);
    semilogy(SNR_dbs, MCRB_p, "Color", colors(p, :), "LineStyle", '--', 'LineWidth', 2); hold on;
end
legend("1", "2", "3", "4 ", "5", "6", "7", "8", "9");
grid on;
xlabel("SNR [dB]");
ylabel("Delay RMSE [s]");
title("OFDM CSI");
theme(gcf, "light");
ylim([5e-11, 1e-8]);

MCRB_OTFS = median(MCRBs_OTFS, 4);

figure;
for p = [1, 5, 3]
    MCRB_p = squeeze(MCRB_OTFS(4*p, 4*p, :));
    MCRB_p = sqrt(MCRB_p);
    semilogy(SNR_dbs, MCRB_p, "Color", colors(p, :), "LineStyle", '--', 'LineWidth', 2); hold on;
end
legend("1", "2", "3", "4 ", "5", "6", "7", "8", "9");
grid on;
xlabel("SNR [dB]");
ylabel("Delay RMSE [s]");
title("OTFS CSI");
theme(gcf, "light");
ylim([5e-11, 1e-8]);

MCRB_OTFS_CE = median(MCRBs_OTFS_CE, 4);

for p_idx = 1:n_pilot
    figure;
    for p = [1, 5, 3]
        MCRB_p = squeeze(MCRB_OTFS_CE(4*p, 4*p, :, p_idx));
        MCRB_p = sqrt(MCRB_p);
        semilogy(SNR_dbs, MCRB_p, "Color", colors(p, :), "LineStyle", '--', 'LineWidth', 2); hold on;
    end
    legend("1", "2", "3", "4 ", "5", "6", "7", "8", "9");
    grid on;
    xlabel("SNR [dB]");
    ylabel("Delay RMSE [s]");
    title(sprintf("OTFS CE (SNR_p = %d dB)", SNR_p_dBs(p_idx)));
    theme(gcf, "light");
    ylim([5e-11, 1e-8]);
end

%% ---- Per-path MCRB Comparison ----
line_styles = {'-', '--', ':'};
scheme_names = {'OFDM (CSI)', 'OTFS (CSI)', 'OTFS (CE)'};
MCRBs_all = {MCRB_OFDM, MCRB_OTFS, MCRB_OTFS_CE};
path_indices = [1, 5, 3];
path_labels  = {'Path 1 (0 dB)', 'Path 5 (-0.6 dB)', 'Path 3 (-1.4 dB)'};

for pi = 1:length(path_indices)
    p = path_indices(pi);
    figure;
    for s = 1:3
        M_cur = MCRBs_all{s};
        MCRB_p = squeeze(M_cur(4*p, 4*p, :));
        MCRB_p = sqrt(MCRB_p);
        semilogy(SNR_dbs, MCRB_p, "LineStyle", line_styles{s}, 'LineWidth', 2); hold on;
    end
    legend(scheme_names, 'Location', 'southwest');
    grid on;
    xlabel("SNR [dB]");
    ylabel("Delay RMSE [s]");
    title(sprintf("MCRB Comparison — %s", path_labels{pi}));
    theme(gcf, "light");
    % ylim([5e-11, 1e-8]);
end

%% ---- Save merged results ----
save('eva_ber_mp_mcrb_merged.mat', 'SNR_dbs', 'P', 'SNR_p_dBs', 'n_pilot', ...
    'bers_ofdm_csi', 'bers_otfs_csi', 'bers_otfs_ce', ...
    'MCRBs_OFDM', 'MCRBs_OTFS', 'MCRBs_OTFS_CE');
