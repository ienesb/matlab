% plot_merged_mcrb.m
% Load results from two MCRB runs, then plot.
clear;

r1 = load('eva_ber_mp_mcrb_results.mat');        % OTFS CSI/CE + OFDM CSI
r2 = load('eva_mcrb_genie_results.mat');         % Genie CRBs
r3 = load('eva_ber_mp_mcrb_ofdm_results.mat');   % OFDM CSI + OFDM CE (new)

SNR_dbs   = r1.SNR_dbs;
P         = r1.P;
SNR_p_dBs = r1.SNR_p_dBs;
n_pilot   = r1.n_pilot;

nMonteCarlo = r1.nMonteCarlo;

% Merge BERs
bers_ofdm_csi = r3.bers_ofdm_csi;
bers_ofdm_ce  = r3.bers_ofdm_ce;
bers_otfs_csi = r1.bers_otfs_csi;
bers_otfs_ce  = r1.bers_otfs_ce;

% Merge MCRBs
MCRBs_OFDM    = r3.MCRBs_OFDM;      % prefer r3 (has matching CE run)
MCRBs_OFDM_CE = r3.MCRBs_OFDM_CE;
MCRBs_OTFS    = r1.MCRBs_OTFS;
MCRBs_OTFS_CE = r1.MCRBs_OTFS_CE;

% Genie CRBs
MCRBs_OFDM_genie = r2.MCRBs_OFDM_genie;
MCRBs_OTFS_genie = r2.MCRBs_OTFS_genie;

%% ---- BER Plot ----
figure;
semilogy(SNR_dbs, bers_ofdm_csi, '-d',  'LineWidth', 2); hold on;
semilogy(SNR_dbs, bers_ofdm_ce,  '-x',  'LineWidth', 2);
semilogy(SNR_dbs, bers_otfs_csi, '-o',  'LineWidth', 2);
markers = {'--s', '--^', '--v'};
for p_idx = 1:n_pilot
    semilogy(SNR_dbs, bers_otfs_ce(:,:,p_idx), markers{p_idx}, 'LineWidth', 2);
end
grid on;
xlabel('SNR_d [dB]');
ylabel('BER');
ce_labels = arrayfun(@(x) sprintf('OTFS (CE + MP, SNR_p=%d dB)', x), SNR_p_dBs, 'UniformOutput', false);
legend(['OFDM (CSI)', 'OFDM (CE, scattered)', 'OTFS (CSI + MP)', ce_labels], 'Location', 'southwest');
title('EVA Channel — Raviteja 2019 Fig. 9');
theme(gcf, "light");

%% ---- MCRB Plots ----
MCRB_OFDM    = median(MCRBs_OFDM,    4);
MCRB_OFDM_CE = median(MCRBs_OFDM_CE, 4);
colors = lines(P);

% OFDM CSI vs OFDM CE overlaid on one figure (solid = CSI, dashed = CE)
figure;
h_lines = gobjects(P, 2);
for p = 1:P
    MCRB_p_csi = sqrt(squeeze(MCRB_OFDM(4*p, 4*p, :)));
    MCRB_p_ce  = sqrt(squeeze(MCRB_OFDM_CE(4*p, 4*p, :)));
    h_lines(p, 1) = semilogy(SNR_dbs, MCRB_p_csi, "Color", colors(p, :), "LineStyle", '-',  'LineWidth', 2); hold on;
    h_lines(p, 2) = semilogy(SNR_dbs, MCRB_p_ce,  "Color", colors(p, :), "LineStyle", '--', 'LineWidth', 2);
end
grid on;
xlabel("SNR [dB]");
ylabel("Delay RMSE [s]");
title("OFDM: CSI (solid) vs CE (dashed)");
path_legend = arrayfun(@(p) sprintf('Path %d', p), 1:P, 'UniformOutput', false);
legend([h_lines(:, 1); h_lines(:, 2)], ...
       [strcat(path_legend, ' (CSI)'), strcat(path_legend, ' (CE)')], ...
       'Location', 'southwest', 'NumColumns', 2);
theme(gcf, "light");
ylim([5e-11, 1e-8]);

MCRB_OTFS = median(MCRBs_OTFS, 4);

figure;
for p = 1:P
    MCRB_p = squeeze(MCRB_OTFS(4*p, 4*p, :));
    MCRB_p = sqrt(MCRB_p);
    semilogy(SNR_dbs, MCRB_p, "Color", colors(p, :), "LineStyle", '--', 'LineWidth', 2); hold on;
end
legend("1", "2", "3");
grid on;
xlabel("SNR [dB]");
ylabel("Delay RMSE [s]");
title("OTFS CSI");
theme(gcf, "light");
ylim([5e-11, 1e-8]);

MCRB_OTFS_CE = median(MCRBs_OTFS_CE, 4);

for p_idx = 1:n_pilot
    figure;
    for p = 1:P
        MCRB_p = squeeze(MCRB_OTFS_CE(4*p, 4*p, :, p_idx));
        MCRB_p = sqrt(MCRB_p);
        semilogy(SNR_dbs, MCRB_p, "Color", colors(p, :), "LineStyle", '--', 'LineWidth', 2); hold on;
    end
    legend("1", "2", "3");
    grid on;
    xlabel("SNR [dB]");
    ylabel("Delay RMSE [s]");
    title(sprintf("OTFS CE (SNR_p = %d dB)", SNR_p_dBs(p_idx)));
    theme(gcf, "light");
    ylim([5e-11, 1e-8]);
end

%% ---- Per-path MCRB Comparison ----
MCRB_OFDM_genie = median(MCRBs_OFDM_genie, 4);
MCRB_OTFS_genie = median(MCRBs_OTFS_genie, 4);

line_styles  = {'-', '--', '-.', ':', '--', '-.'};
scheme_names = {'OFDM (CSI)', 'OFDM (CE)', 'OTFS (CSI)', 'OTFS (CE)', 'OFDM (Genie)', 'OTFS (Genie)'};
MCRBs_all    = {MCRB_OFDM, MCRB_OFDM_CE, MCRB_OTFS, MCRB_OTFS_CE, MCRB_OFDM_genie, MCRB_OTFS_genie};
path_indices = 1:P;
path_labels  = {'Path 1 (0 dB)', 'Path 5 (-0.6 dB)', 'Path 3 (-1.4 dB)'};

for pi = 1:length(path_indices)
    p = path_indices(pi);
    figure;
    for s = 1:length(MCRBs_all)
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
save('eva_ber_mp_mcrb_merged.mat', ...
    'SNR_dbs', 'P', 'SNR_p_dBs', 'n_pilot', 'nMonteCarlo', ...
    'bers_ofdm_csi', 'bers_ofdm_ce', 'bers_otfs_csi', 'bers_otfs_ce', ...
    'MCRBs_OFDM', 'MCRBs_OFDM_CE', 'MCRBs_OTFS', 'MCRBs_OTFS_CE', ...
    'MCRBs_OFDM_genie', 'MCRBs_OTFS_genie', ...
    'MCRB_OFDM', 'MCRB_OFDM_CE', 'MCRB_OTFS', 'MCRB_OTFS_CE', ...
    'MCRB_OFDM_genie', 'MCRB_OTFS_genie');
