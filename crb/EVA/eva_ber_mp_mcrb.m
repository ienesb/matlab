% eva_ber_mp.m
% BER for EVA channel: OFDM (CSI), OTFS (CSI + MP), OTFS (CE + MP)
% Reproduces Raviteja et al. IEEE TVT 2019, Fig. 9
% Fractional Doppler with full guard symbols.
clear;

initialize_parameters;

n_iter_mp   = 30;
damping     = 0.5;
constellation = qammod((0:M_order-1).', M_order, UnitAveragePower=true);

% ---- Pilot design (fractional Doppler, full guard — TVT 2019 Eq. 8) ----
% SNR_p_dBs = [40, 45, 50];
SNR_p_dBs = [45];
n_pilot   = length(SNR_p_dBs);

l_p = l_tau;
k_p = floor(M/2);

guard_row_lo = l_p - l_tau;    % = 0
guard_row_hi = l_p + l_tau;    % = 40

data_mask = true(N, M);
data_mask(guard_row_lo+1 : guard_row_hi+1, :) = false;

ce_thresholds = 3 ./ sqrt(10.^(SNR_p_dBs/10));
csi_rel_dB = 30;
ce_rows = (l_p + 1) : (guard_row_hi + 1);  % MATLAB rows 21..41

% Pre-allocate
nSNR = length(SNR_dbs);
bers_ofdm_csi = zeros(nSNR, nMonteCarlo);
bers_otfs_csi = zeros(nSNR, nMonteCarlo);
bers_otfs_ce  = zeros(nSNR, nMonteCarlo, n_pilot);

MCRBs_OTFS = zeros(4*P, 4*P, nSNR, nMonteCarlo);
MCRBs_OFDM = zeros(4*P, 4*P, nSNR, nMonteCarlo);
MCRBs_OTFS_CE = zeros(4*P, 4*P, nSNR, nMonteCarlo, n_pilot);

total_duration = 0;

Z_unit = zeros(N, M, nMonteCarlo);
for mc_idx = 1:nMonteCarlo
    Z_unit(:, :, mc_idx) = (randn(N, M) + 1j*randn(N, M)) / sqrt(2);
end

for SNR_idx = 1:nSNR
    SNR_db  = SNR_dbs(SNR_idx);
    SNR_lin = db2pow(SNR_db);
    N0      = 1 / SNR_lin;

    % rng(2026);
    tic
    parfor mc_idx = 1:nMonteCarlo
        % ---- Common: channel realization + data ----
        current_doppler    = nu_max * cos(2*pi*rand(1, P) - pi);
        current_path_gains = path_gains_amp .* exp(1j*rand(1, P)*2*pi);
        symbols = randi(M_order, N, M) - 1;
        X       = qammod(symbols, M_order, UnitAveragePower=true);

        % ---- Scale pre-generated noise for this SNR ----
        Z_TF = Z_unit(:, :, mc_idx) * sqrt(N0);
        Z_DD = sfft(Z_TF);

        %% ======== OFDM (CSI) — MMSE equalization ========
        
        mu_gt_ofdm = channel(X, current_path_gains, time_delays, ...
                            current_doppler, deltaf, Ts, 0);
        Y_TF_ofdm = mu_gt_ofdm + Z_TF;
        % mu_gt_ofdm = mu_gt_ofdm(:);

        X_hat_ofdm = channel_inverse(Y_TF_ofdm, current_path_gains, ...
                                     time_delays, current_doppler, ...
                                     deltaf, Ts, N0);

        tx_bits_all  = qamdemod(X(:),          M_order, OutputType="bit", UnitAveragePower=true);
        rx_bits_ofdm = qamdemod(X_hat_ofdm(:), M_order, OutputType="bit", UnitAveragePower=true);
        
        X_hat_ofdm = qammod(qamdemod(X_hat_ofdm, M_order, UnitAveragePower=true), M_order, UnitAveragePower=true);
        
        MCRBs_OFDM(:, :, SNR_idx, mc_idx) = getMCRB_TF(X_hat_ofdm, X, deltaf, Ts, N0, current_path_gains, current_doppler, time_delays, P);
        
        bers_ofdm_csi(SNR_idx, mc_idx) = mean(tx_bits_all ~= rx_bits_ofdm);

        %% ======== Noiseless OTFS responses (computed once) ========
        % Data-only response
        X_DD_data = zeros(N, M);
        X_DD_data(data_mask) = X(data_mask);
        Y_DD_data_clean = sfft(channel(isfft(X_DD_data), current_path_gains, ...
                               time_delays, current_doppler, deltaf, Ts, 0));

        % DD impulse response (for CSI and pilot contribution)
        X_impulse = zeros(N, M);
        X_impulse(l_p+1, k_p+1) = 1;
        h_dd_full = sfft(channel(isfft(X_impulse), current_path_gains, ...
                         time_delays, current_doppler, deltaf, Ts, 0));

        %% ======== OTFS CSI + MP ========
        x_p_max = sqrt(10^(max(SNR_p_dBs)/10) * N0);
        Y_DD = Y_DD_data_clean + x_p_max * h_dd_full + Z_DD;

        % Threshold h_dd for CSI
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
        X_hat_csi_hard = qammod(qamdemod(X_hat_csi, M_order, UnitAveragePower=true), M_order, UnitAveragePower=true);
        X_hat_csi_hard(~data_mask) = 0;
        X_hat_csi_hard = X_hat_csi_hard + x_p_max*X_impulse;
        X_hat_csi_TF = isfft(X_hat_csi_hard);

        X_csi_TF = isfft(X_DD_data + x_p_max*X_impulse);
        % X_TF_gt = isfft(X_DD_data + x_p_max*X_impulse);
        % mu_gt_otfs = channel(X_TF_gt, current_path_gains, time_delays, ...
        %                      current_doppler, deltaf, Ts, 0);
        % mu_gt_otfs = mu_gt_otfs(:);
        
        MCRBs_OTFS(:, :, SNR_idx, mc_idx) = getMCRB_TF(X_hat_csi_TF, X_csi_TF, ...
            deltaf, Ts, N0, current_path_gains, current_doppler, time_delays, P);

        tx_bits     = qamdemod(X(data_mask),         M_order, OutputType="bit", UnitAveragePower=true);
        rx_bits_csi = qamdemod(X_hat_csi(data_mask), M_order, OutputType="bit", UnitAveragePower=true);
        bers_otfs_csi(SNR_idx, mc_idx) = mean(tx_bits ~= rx_bits_csi);

        %% ======== OTFS CE + MP for each pilot SNR ========
        ber_ce_mc = zeros(1, n_pilot);
        for p_idx = 1:n_pilot
            x_p = sqrt(10^(SNR_p_dBs(p_idx)/10) * N0);

            % Build Y_DD with this pilot power (same noise, same data)
            Y_DD_ce = Y_DD_data_clean + x_p * h_dd_full + Z_DD;

            % Channel estimation: threshold on guard region
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

            X_hat_ce_hard = qammod(qamdemod(X_hat_ce, M_order, UnitAveragePower=true), M_order, UnitAveragePower=true);
            X_hat_ce_hard(~data_mask) = 0;
            X_hat_ce_hard = X_hat_ce_hard + x_p*X_impulse;
            X_hat_ce_TF = isfft(X_hat_ce_hard);
            X_ce_TF = isfft(X_DD_data + x_p*X_impulse);
            MCRBs_OTFS_CE(:, :, SNR_idx, mc_idx, p_idx) = getMCRB_TF(X_hat_ce_TF, X_ce_TF, ...
                deltaf, Ts, N0, current_path_gains, current_doppler, time_delays, P);

            rx_bits_ce  = qamdemod(X_hat_ce(data_mask), M_order, OutputType="bit", UnitAveragePower=true);
            ber_ce_mc(p_idx) = mean(tx_bits ~= rx_bits_ce);
        end
        bers_otfs_ce(SNR_idx, mc_idx, :) = ber_ce_mc;
    end
    duration       = toc;
    total_duration = total_duration + duration;
    remaining_time = (nSNR - SNR_idx) * duration;
    fprintf("%d/%d | OFDM=%.3e CSI=%.3e", SNR_idx, nSNR, ...
        mean(bers_ofdm_csi(SNR_idx,:)), mean(bers_otfs_csi(SNR_idx,:)));
    for p_idx = 1:n_pilot
        fprintf(" CE%d=%.3e", SNR_p_dBs(p_idx), mean(bers_otfs_ce(SNR_idx,:,p_idx)));
    end
    fprintf(" | %.1fs (remain: %.1fs)\n", duration, remaining_time);
end

fprintf("Completed in %.1f seconds\n", total_duration);

bers_ofdm_csi = mean(bers_ofdm_csi, 2);
bers_otfs_csi = mean(bers_otfs_csi, 2);
bers_otfs_ce  = mean(bers_otfs_ce,  2);   % [nSNR x 1 x n_pilot]

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



MCRB_OFDM = median(MCRBs_OFDM, 4);

colors = lines(P);

figure;
for p = 1:7
    MCRB_p = MCRB_OFDM(4*p, 4*p, :);
    MCRB_p = squeeze(MCRB_p);
    MCRB_p = sqrt(MCRB_p);

    semilogy(SNR_dbs, MCRB_p,"Color", colors(p, :), "LineStyle", '--', 'LineWidth', 2); hold on;
end
legend("1", "2", "3", "4 " ,"5", "6", "7", "8", "9");
grid on;
xlabel("SNR [dB]");
ylabel("Delay RMSE [s]");
title("OFDM Full Data");
theme(gcf, "light");


MCRB_OTFS = median(MCRBs_OTFS, 4);

colors = lines(P);

figure;
for p = 1:7
    MCRB_p = MCRB_OTFS(4*p, 4*p, :);
    MCRB_p = squeeze(MCRB_p);
    MCRB_p = sqrt(MCRB_p);

    semilogy(SNR_dbs, MCRB_p,"Color", colors(p, :), "LineStyle", '--', 'LineWidth', 2); hold on;
end
legend("1", "2", "3", "4 " ,"5", "6", "7", "8", "9");
grid on;
xlabel("SNR [dB]");
ylabel("Delay RMSE [s]");
title("OTFS Full Data");
theme(gcf, "light");

MCRB_OTFS_CE = median(MCRBs_OTFS_CE, 4);

for p_idx = 1:n_pilot
    figure;
    for p = 1:7
        MCRB_p = MCRB_OTFS_CE(4*p, 4*p, :, p_idx);
        MCRB_p = squeeze(MCRB_p);
        MCRB_p = sqrt(MCRB_p);

        semilogy(SNR_dbs, MCRB_p, "Color", colors(p, :), "LineStyle", '--', 'LineWidth', 2); hold on;
    end
    legend("1", "2", "3", "4 ", "5", "6", "7", "8", "9");
    grid on;
    xlabel("SNR [dB]");
    ylabel("Delay RMSE [s]");
    title(sprintf("OTFS CE (SNR_p = %d dB)", SNR_p_dBs(p_idx)));
    theme(gcf, "light");
end