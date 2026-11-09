% eva_ber_mp_mcrb_dd.m
% BER + MCRB for OTFS only (CSI + MP, CE + MP) using getMCRB_DD,
% which computes the MCRB directly in the DD domain from
% otfs_io_derivatives.tex.
clear;

initialize_parameters;

n_iter_mp     = 30;
damping       = 0.5;
constellation = qammod((0:M_order-1).', M_order, UnitAveragePower=true);

% ---- Pilot design (fractional Doppler, full guard — TVT 2019 Eq. 8) ----
SNR_p_dBs = [45];
n_pilot   = length(SNR_p_dBs);

l_p = l_tau;
k_p = floor(M/2);

guard_row_lo = l_p - l_tau;
guard_row_hi = l_p + l_tau;

data_mask = true(N, M);
data_mask(guard_row_lo+1 : guard_row_hi+1, :) = false;

ce_thresholds = 3 ./ sqrt(10.^(SNR_p_dBs/10));
csi_rel_dB    = 30;
ce_rows       = (l_p + 1) : (guard_row_hi + 1);

% Pre-allocate
nSNR = length(SNR_dbs);
bers_otfs_csi = zeros(nSNR, nMonteCarlo);
bers_otfs_ce  = zeros(nSNR, nMonteCarlo, n_pilot);

MCRBs_OTFS_CSI = zeros(4*P, 4*P, nSNR, nMonteCarlo);
MCRBs_OTFS_CE  = zeros(4*P, 4*P, nSNR, nMonteCarlo, n_pilot);

total_duration = 0;

rng(2026);
Z_unit = zeros(N, M, nMonteCarlo);
for mc_idx = 1:nMonteCarlo
    Z_unit(:, :, mc_idx) = (randn(N, M) + 1j*randn(N, M)) / sqrt(2);
end

for SNR_idx = 1:nSNR
    SNR_db  = SNR_dbs(SNR_idx);
    SNR_lin = db2pow(SNR_db);
    N0      = 1 / SNR_lin;

    tic
    parfor mc_idx = 1:nMonteCarlo
        % ---- Channel realization + data ----
        current_doppler    = nu_max * cos(2*pi*rand(1, P) - pi);
        current_path_gains = path_gains_amp .* exp(1j*rand(1, P)*2*pi);
        symbols = randi(M_order, N, M) - 1;
        X       = qammod(symbols, M_order, UnitAveragePower=true);

        % Shared noise across CSI/CE variants
        Z_TF = Z_unit(:, :, mc_idx) * sqrt(N0);
        Z_DD = sfft(Z_TF);

        %% ======== Noiseless OTFS responses ========
        X_DD_data = zeros(N, M);
        X_DD_data(data_mask) = X(data_mask);
        Y_DD_data_clean = sfft(channel(isfft(X_DD_data), current_path_gains, ...
                               time_delays, current_doppler, deltaf, Ts, 0));

        X_impulse = zeros(N, M);
        X_impulse(l_p+1, k_p+1) = 1;
        h_dd_full = sfft(channel(isfft(X_impulse), current_path_gains, ...
                         time_delays, current_doppler, deltaf, Ts, 0));

        %% ======== OTFS CSI + MP ========
        x_p_max = sqrt(10^(max(SNR_p_dBs)/10) * N0);
        Y_DD = Y_DD_data_clean + x_p_max * h_dd_full + Z_DD;

        max_h      = max(abs(h_dd_full(:)));
        csi_thresh = max_h * 10^(-csi_rel_dB/20);
        h_dd_csi   = h_dd_full;
        h_dd_csi(abs(h_dd_csi) < csi_thresh) = 0;

        [r_csi, c_csi, g_csi] = find(h_dd_csi);
        dl_csi = mod(r_csi - (l_p + 1), N);
        dk_csi = mod(c_csi - (k_p + 1), M);

        known_sym = zeros(N, M);
        known_sym(l_p+1, k_p+1) = x_p_max;

        X_hat_csi = mp_detector(Y_DD, dl_csi, dk_csi, g_csi, ...
                                N0, constellation, n_iter_mp, ...
                                data_mask, known_sym, damping);

        % Hard decisions (DD domain) — replace soft data with constellation
        % symbols and reinsert pilot.  This is the X_hat fed to getMCRB_DD.
        X_hat_csi_hard = qammod(qamdemod(X_hat_csi, M_order, UnitAveragePower=true), ...
                                M_order, UnitAveragePower=true);
        X_hat_csi_hard(~data_mask) = 0;
        X_hat_csi_hard = X_hat_csi_hard + x_p_max*X_impulse;

        X_csi_DD = X_DD_data + x_p_max*X_impulse;

        MCRBs_OTFS_CSI(:, :, SNR_idx, mc_idx) = getMCRB_DD( ...
            X_hat_csi_hard, X_csi_DD, deltaf, Ts, N0, ...
            current_path_gains, current_doppler, time_delays, P);

        tx_bits     = qamdemod(X(data_mask),         M_order, OutputType="bit", UnitAveragePower=true);
        rx_bits_csi = qamdemod(X_hat_csi(data_mask), M_order, OutputType="bit", UnitAveragePower=true);
        bers_otfs_csi(SNR_idx, mc_idx) = mean(tx_bits ~= rx_bits_csi);

        %% ======== OTFS CE + MP for each pilot SNR ========
        ber_ce_mc   = zeros(1, n_pilot);
        mcrb_ce_mc  = zeros(4*P, 4*P, n_pilot);
        for p_idx = 1:n_pilot
            x_p = sqrt(10^(SNR_p_dBs(p_idx)/10) * N0);

            Y_DD_ce = Y_DD_data_clean + x_p * h_dd_full + Z_DD;

            h_dd_ce = zeros(N, M);
            h_dd_ce(ce_rows, :) = Y_DD_ce(ce_rows, :) / x_p;
            h_dd_ce(abs(h_dd_ce) < ce_thresholds(p_idx)) = 0;

            [r_ce, c_ce, g_ce] = find(h_dd_ce);
            dl_ce = mod(r_ce - (l_p + 1), N);
            dk_ce = mod(c_ce - (k_p + 1), M);

            known_sym_ce = zeros(N, M);
            known_sym_ce(l_p+1, k_p+1) = x_p;

            X_hat_ce = mp_detector(Y_DD_ce, dl_ce, dk_ce, g_ce, ...
                                   N0, constellation, n_iter_mp, ...
                                   data_mask, known_sym_ce, damping);

            X_hat_ce_hard = qammod(qamdemod(X_hat_ce, M_order, UnitAveragePower=true), ...
                                   M_order, UnitAveragePower=true);
            X_hat_ce_hard(~data_mask) = 0;
            X_hat_ce_hard = X_hat_ce_hard + x_p*X_impulse;

            X_ce_DD = X_DD_data + x_p*X_impulse;

            mcrb_ce_mc(:, :, p_idx) = getMCRB_DD( ...
                X_hat_ce_hard, X_ce_DD, deltaf, Ts, N0, ...
                current_path_gains, current_doppler, time_delays, P);

            rx_bits_ce       = qamdemod(X_hat_ce(data_mask), M_order, OutputType="bit", UnitAveragePower=true);
            ber_ce_mc(p_idx) = mean(tx_bits ~= rx_bits_ce);
        end
        bers_otfs_ce(SNR_idx, mc_idx, :)     = ber_ce_mc;
        MCRBs_OTFS_CE(:, :, SNR_idx, mc_idx, :) = mcrb_ce_mc;
    end
    duration       = toc;
    total_duration = total_duration + duration;
    remaining_time = (nSNR - SNR_idx) * duration;
    fprintf("%d/%d | CSI=%.3e", SNR_idx, nSNR, mean(bers_otfs_csi(SNR_idx,:)));
    for p_idx = 1:n_pilot
        fprintf(" CE%d=%.3e", SNR_p_dBs(p_idx), mean(bers_otfs_ce(SNR_idx,:,p_idx)));
    end
    fprintf(" | %.1fs (remain: %.1fs)\n", duration, remaining_time);
    save eva_ber_mp_mcrb_dd_results.mat
end

fprintf("Completed in %.1f seconds\n", total_duration);

save eva_ber_mp_mcrb_dd_resultsv2.mat

return;

%% ===== Plotting =====
bers_otfs_csi_avg = mean(bers_otfs_csi, 2);
bers_otfs_ce_avg  = mean(bers_otfs_ce,  2);

figure;
semilogy(SNR_dbs, bers_otfs_csi_avg, '-o', 'LineWidth', 2); hold on;
markers = {'--s', '--^', '--v'};
for p_idx = 1:n_pilot
    semilogy(SNR_dbs, squeeze(bers_otfs_ce_avg(:,:,p_idx)), markers{p_idx}, 'LineWidth', 2);
end
grid on;
xlabel('SNR_d [dB]');
ylabel('BER');
ce_labels = arrayfun(@(x) sprintf('OTFS (CE + MP, SNR_p=%d dB)', x), SNR_p_dBs, 'UniformOutput', false);
legend(['OTFS (CSI + MP)', ce_labels], 'Location', 'southwest');
title('OTFS BER — getMCRB\_DD (DD-domain MCRB)');

%% MCRB delay RMSE per path
colors = lines(P);

MCRB_OTFS_CSI = median(MCRBs_OTFS_CSI, 4);
figure;
for p = 1:P
    MCRB_p = squeeze(MCRB_OTFS_CSI(4*p, 4*p, :));
    semilogy(SNR_dbs, sqrt(MCRB_p), 'Color', colors(p,:), 'LineStyle', '--', 'LineWidth', 2); hold on;
end
legend(arrayfun(@(p) sprintf('Path %d', p), 1:P, 'UniformOutput', false));
grid on;
xlabel('SNR [dB]');
ylabel('Delay RMSE [s]');
title('OTFS CSI — MCRB (DD)');

MCRB_OTFS_CE = median(MCRBs_OTFS_CE, 4);
for p_idx = 1:n_pilot
    figure;
    for p = 1:P
        MCRB_p = squeeze(MCRB_OTFS_CE(4*p, 4*p, :, p_idx));
        semilogy(SNR_dbs, sqrt(MCRB_p), 'Color', colors(p,:), 'LineStyle', '--', 'LineWidth', 2); hold on;
    end
    legend(arrayfun(@(p) sprintf('Path %d', p), 1:P, 'UniformOutput', false));
    grid on;
    xlabel('SNR [dB]');
    ylabel('Delay RMSE [s]');
    title(sprintf('OTFS CE — MCRB (DD), SNR_p = %d dB', SNR_p_dBs(p_idx)));
end
