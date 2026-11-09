% eva_ber_mp_mcrb_mle.m
% EVA channel: single-pass script that computes BER, MCRB, and MLE MSE
% for OFDM and OTFS in one Monte Carlo loop. Merges the pipelines of
% eva_ber_mp_mcrb.m and eva_ber_mp_mle.m so every trial reuses the same
% channel realization, noise, data, and detected symbols.
clear;

initialize_parameters;

n_iter_mp     = 30;
damping       = 0.5;
constellation = qammod((0:M_order-1).', M_order, UnitAveragePower=true);

% ---- OTFS pilot design (fractional Doppler, full guard - TVT 2019 Eq. 8) ----
SNR_p_dBs = 45;
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

% ---- OFDM scattered pilot pattern ----
dN_p = 8;
dM_p = 8;
n_pilots_ofdm   = 1:dN_p:N;
m_pilots_ofdm   = 1:dM_p:M;
pilot_mask_ofdm = false(N, M);
pilot_mask_ofdm(n_pilots_ofdm, m_pilots_ofdm) = true;
data_mask_ofdm  = ~pilot_mask_ofdm;

nSNR = length(SNR_dbs);

% ---- BER storage ----
bers_ofdm_csi = zeros(nSNR, nMonteCarlo);
bers_ofdm_ce  = zeros(nSNR, nMonteCarlo);
bers_otfs_csi = zeros(nSNR, nMonteCarlo);
bers_otfs_ce  = zeros(nSNR, nMonteCarlo, n_pilot);

% ---- MCRB storage ----
MCRBs_OFDM       = zeros(4*P, 4*P, nSNR, nMonteCarlo);
MCRBs_OFDM_CE    = zeros(4*P, 4*P, nSNR, nMonteCarlo);
MCRBs_OFDM_genie = zeros(4*P, 4*P, nSNR, nMonteCarlo);
MCRBs_OTFS       = zeros(4*P, 4*P, nSNR, nMonteCarlo);
MCRBs_OTFS_CE    = zeros(4*P, 4*P, nSNR, nMonteCarlo, n_pilot);
MCRBs_OTFS_genie = zeros(4*P, 4*P, nSNR, nMonteCarlo);

% ---- MLE MSE storage (nSNR x nMC x P) ----
mse_tau_ofdm_csi   = zeros(nSNR, nMonteCarlo, P);
mse_nu_ofdm_csi    = zeros(nSNR, nMonteCarlo, P);
mse_tau_ofdm_ce    = zeros(nSNR, nMonteCarlo, P);
mse_nu_ofdm_ce     = zeros(nSNR, nMonteCarlo, P);
mse_tau_ofdm_genie = zeros(nSNR, nMonteCarlo, P);
mse_nu_ofdm_genie  = zeros(nSNR, nMonteCarlo, P);

mse_tau_otfs_csi   = zeros(nSNR, nMonteCarlo, P);
mse_nu_otfs_csi    = zeros(nSNR, nMonteCarlo, P);
mse_tau_otfs_genie = zeros(nSNR, nMonteCarlo, P);
mse_nu_otfs_genie  = zeros(nSNR, nMonteCarlo, P);
mse_tau_otfs_ce    = zeros(nSNR, nMonteCarlo, P, n_pilot);
mse_nu_otfs_ce     = zeros(nSNR, nMonteCarlo, P, n_pilot);

% ---- Pre-generated unit-variance TF noise, scaled per SNR point ----
Z_unit = zeros(N, M, nMonteCarlo);
for mc_idx = 1:nMonteCarlo
    Z_unit(:, :, mc_idx) = (randn(N, M) + 1j*randn(N, M)) / sqrt(2);
end

total_duration = 0;

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

        Z_TF = Z_unit(:, :, mc_idx) * sqrt(N0);
        Z_DD = sfft(Z_TF);

        %% ======== OFDM (CSI) ========
        Y_TF_ofdm = channel(X, current_path_gains, time_delays, ...
                            current_doppler, deltaf, Ts, 0) + Z_TF;

        X_eq_ofdm = channel_inverse(Y_TF_ofdm, current_path_gains, ...
                                    time_delays, current_doppler, ...
                                    deltaf, Ts, N0);

        tx_bits_all  = qamdemod(X(:),         M_order, OutputType="bit", UnitAveragePower=true);
        rx_bits_ofdm = qamdemod(X_eq_ofdm(:), M_order, OutputType="bit", UnitAveragePower=true);
        bers_ofdm_csi(SNR_idx, mc_idx) = mean(tx_bits_all ~= rx_bits_ofdm);

        X_hat_ofdm = qammod(qamdemod(X_eq_ofdm, M_order, UnitAveragePower=true), ...
                            M_order, UnitAveragePower=true);

        MCRBs_OFDM(:, :, SNR_idx, mc_idx) = getMCRB_TF( ...
            X_hat_ofdm, X, deltaf, Ts, N0, current_path_gains, ...
            current_doppler, time_delays, P);
        MCRBs_OFDM_genie(:, :, SNR_idx, mc_idx) = getMCRB_TF_genie( ...
            X, deltaf, Ts, N0, current_path_gains, current_doppler, ...
            time_delays, P);

        %% ======== OFDM (CE) ========
        H_ls = Y_TF_ofdm(n_pilots_ofdm, m_pilots_ofdm) ./ X(n_pilots_ofdm, m_pilots_ofdm);
        F_H  = griddedInterpolant({n_pilots_ofdm, m_pilots_ofdm}, H_ls, 'linear', 'nearest');
        [Nq, Mq] = ndgrid(1:N, 1:M);
        H_hat_ce = F_H(Nq, Mq);

        X_eq_ofdm_ce = conj(H_hat_ce) .* Y_TF_ofdm ./ (abs(H_hat_ce).^2 + N0);

        tx_bits_ce_ofdm = qamdemod(X(data_mask_ofdm),            M_order, OutputType="bit", UnitAveragePower=true);
        rx_bits_ce_ofdm = qamdemod(X_eq_ofdm_ce(data_mask_ofdm), M_order, OutputType="bit", UnitAveragePower=true);
        bers_ofdm_ce(SNR_idx, mc_idx) = mean(tx_bits_ce_ofdm ~= rx_bits_ce_ofdm);

        X_hat_ofdm_ce = qammod(qamdemod(X_eq_ofdm_ce, M_order, UnitAveragePower=true), ...
                               M_order, UnitAveragePower=true);
        X_hat_ofdm_ce(pilot_mask_ofdm) = X(pilot_mask_ofdm);

        MCRBs_OFDM_CE(:, :, SNR_idx, mc_idx) = getMCRB_TF( ...
            X_hat_ofdm_ce, X, deltaf, Ts, N0, current_path_gains, ...
            current_doppler, time_delays, P);

        %% ======== OFDM MLE ========
        [~, tau_csi, nu_csi] = mle_tf_gt(Y_TF, X_hat_csi, deltaf, Ts, time_delays, current_doppler);
        [~, tau_ce,  nu_ce ] = mle_tf_gt(Y_TF, X_hat_ce,  deltaf, Ts, time_delays, current_doppler);
        [~, tau_gen, nu_gen] = mle_tf_gt(Y_TF, X,         deltaf, Ts, time_delays, current_doppler);

        perm_csi = match_paths(time_delays, current_doppler, tau_csi, nu_csi);
        perm_ce  = match_paths(time_delays, current_doppler, tau_ce,  nu_ce );
        perm_gen = match_paths(time_delays, current_doppler, tau_gen, nu_gen);

        mse_tau_ofdm_csi  (SNR_idx, mc_idx, :) = (time_delays(:)     - tau_csi(perm_csi)).^2;
        mse_nu_ofdm_csi   (SNR_idx, mc_idx, :) = (current_doppler(:) - nu_csi (perm_csi)).^2;
        mse_tau_ofdm_ce   (SNR_idx, mc_idx, :) = (time_delays(:)     - tau_ce (perm_ce )).^2;
        mse_nu_ofdm_ce    (SNR_idx, mc_idx, :) = (current_doppler(:) - nu_ce  (perm_ce )).^2;
        mse_tau_ofdm_genie(SNR_idx, mc_idx, :) = (time_delays(:)     - tau_gen(perm_gen)).^2;
        mse_nu_ofdm_genie (SNR_idx, mc_idx, :) = (current_doppler(:) - nu_gen (perm_gen)).^2;

        %% ======== Noiseless OTFS responses (computed once) ========
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
        Y_DD    = Y_DD_data_clean + x_p_max * h_dd_full + Z_DD;

        max_h      = max(abs(h_dd_full(:)));
        csi_thresh = max_h * 10^(-csi_rel_dB/20);
        h_dd_csi   = h_dd_full;
        h_dd_csi(abs(h_dd_csi) < csi_thresh) = 0;

        [r_csi, c_csi, g_csi] = find(h_dd_csi);
        dl_csi = mod(r_csi - (l_p + 1), N);
        dk_csi = mod(c_csi - (k_p + 1), M);

        known_sym = zeros(N, M);
        known_sym(l_p+1, k_p+1) = x_p_max;

        X_hat_csi_dd = mp_detector(Y_DD, dl_csi, dk_csi, g_csi, ...
                                   N0, constellation, n_iter_mp, ...
                                   data_mask, known_sym, damping);

        tx_bits_otfs  = qamdemod(X(data_mask),            M_order, OutputType="bit", UnitAveragePower=true);
        rx_bits_csi   = qamdemod(X_hat_csi_dd(data_mask), M_order, OutputType="bit", UnitAveragePower=true);
        bers_otfs_csi(SNR_idx, mc_idx) = mean(tx_bits_otfs ~= rx_bits_csi);

        X_hat_csi_hard = qammod(qamdemod(X_hat_csi_dd, M_order, UnitAveragePower=true), ...
                                M_order, UnitAveragePower=true);
        X_hat_csi_hard(~data_mask) = 0;
        X_hat_csi_hard = X_hat_csi_hard + x_p_max * X_impulse;
        X_hat_csi_TF   = isfft(X_hat_csi_hard);

        X_csi_TF_genie = isfft(X_DD_data + x_p_max * X_impulse);

        MCRBs_OTFS(:, :, SNR_idx, mc_idx) = getMCRB_TF( ...
            X_hat_csi_TF, X_csi_TF_genie, deltaf, Ts, N0, ...
            current_path_gains, current_doppler, time_delays, P);
        MCRBs_OTFS_genie(:, :, SNR_idx, mc_idx) = getMCRB_TF_genie( ...
            X_csi_TF_genie, deltaf, Ts, N0, current_path_gains, ...
            current_doppler, time_delays, P);

        %% ======== OTFS MLE (CSI + genie) ========
        Y_TF_otfs = isfft(Y_DD);

        [~, tau_ocsi, nu_ocsi] = mle_tf_gt(Y_TF_otfs, X_hat_csi_TF,  deltaf, Ts, time_delays, current_doppler);
        [~, tau_ogen, nu_ogen] = mle_tf_gt(Y_TF_otfs, X_csi_TF_genie, deltaf, Ts, time_delays, current_doppler);

        perm_ocsi = match_paths(time_delays, current_doppler, tau_ocsi, nu_ocsi);
        perm_ogen = match_paths(time_delays, current_doppler, tau_ogen, nu_ogen);

        mse_tau_otfs_csi  (SNR_idx, mc_idx, :) = (time_delays(:)     - tau_ocsi(perm_ocsi)).^2;
        mse_nu_otfs_csi   (SNR_idx, mc_idx, :) = (current_doppler(:) - nu_ocsi (perm_ocsi)).^2;
        mse_tau_otfs_genie(SNR_idx, mc_idx, :) = (time_delays(:)     - tau_ogen(perm_ogen)).^2;
        mse_nu_otfs_genie (SNR_idx, mc_idx, :) = (current_doppler(:) - nu_ogen (perm_ogen)).^2;

        %% ======== OTFS CE + MP (per pilot SNR) ========
        ber_ce_mc      = zeros(1, n_pilot);
        mcrb_ce_mc     = zeros(4*P, 4*P, n_pilot);
        mse_tau_ce_p   = zeros(1, n_pilot, P);
        mse_nu_ce_p    = zeros(1, n_pilot, P);
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

            X_hat_ce_dd = mp_detector(Y_DD_ce, dl_ce, dk_ce, g_ce, ...
                                      N0, constellation, n_iter_mp, ...
                                      data_mask, known_sym_ce, damping);

            rx_bits_ce        = qamdemod(X_hat_ce_dd(data_mask), M_order, OutputType="bit", UnitAveragePower=true);
            ber_ce_mc(p_idx)  = mean(tx_bits_otfs ~= rx_bits_ce);

            X_hat_ce_hard = qammod(qamdemod(X_hat_ce_dd, M_order, UnitAveragePower=true), ...
                                   M_order, UnitAveragePower=true);
            X_hat_ce_hard(~data_mask) = 0;
            X_hat_ce_hard = X_hat_ce_hard + x_p * X_impulse;
            X_hat_ce_TF   = isfft(X_hat_ce_hard);

            X_ce_TF_genie = isfft(X_DD_data + x_p * X_impulse);

            mcrb_ce_mc(:, :, p_idx) = getMCRB_TF( ...
                X_hat_ce_TF, X_ce_TF_genie, deltaf, Ts, N0, ...
                current_path_gains, current_doppler, time_delays, P);

            Y_TF_ce = isfft(Y_DD_ce);
            [~, tau_oce, nu_oce] = mle_tf_gt(Y_TF_ce, X_hat_ce_TF, deltaf, Ts, time_delays, current_doppler);
            perm_oce = match_paths(time_delays, current_doppler, tau_oce, nu_oce);

            mse_tau_ce_p(1, p_idx, :) = (time_delays(:)     - tau_oce(perm_oce)).^2;
            mse_nu_ce_p (1, p_idx, :) = (current_doppler(:) - nu_oce (perm_oce)).^2;
        end
        bers_otfs_ce     (SNR_idx, mc_idx, :)    = ber_ce_mc;
        MCRBs_OTFS_CE    (:, :, SNR_idx, mc_idx, :) = mcrb_ce_mc;
        mse_tau_otfs_ce  (SNR_idx, mc_idx, :, :) = permute(mse_tau_ce_p, [1 3 2]);
        mse_nu_otfs_ce   (SNR_idx, mc_idx, :, :) = permute(mse_nu_ce_p,  [1 3 2]);
    end
    duration       = toc;
    total_duration = total_duration + duration;
    remaining      = (nSNR - SNR_idx) * duration;
    fprintf("%d/%d | OFDM csi=%.2e ce=%.2e | OTFS csi=%.2e", ...
        SNR_idx, nSNR, mean(bers_ofdm_csi(SNR_idx,:)), ...
        mean(bers_ofdm_ce(SNR_idx,:)), mean(bers_otfs_csi(SNR_idx,:)));
    for p_idx = 1:n_pilot
        fprintf(" ce%d=%.2e", SNR_p_dBs(p_idx), mean(bers_otfs_ce(SNR_idx,:,p_idx)));
    end
    fprintf(" | %.1fs (remain %.1fs)\n", duration, remaining);
end

fprintf("Completed in %.1f seconds\n", total_duration);

save eva_ber_mp_mcrb_mle_results.mat

return;

%% ================ Plotting ================
% Trial-mean MCRB diagonals
MCRB_OFDM_mean       = mean(MCRBs_OFDM,       4);
MCRB_OFDM_CE_mean    = mean(MCRBs_OFDM_CE,    4);
MCRB_OFDM_genie_mean = mean(MCRBs_OFDM_genie, 4);
MCRB_OTFS_mean       = mean(MCRBs_OTFS,       4);
MCRB_OTFS_CE_mean    = mean(MCRBs_OTFS_CE,    4);
MCRB_OTFS_genie_mean = mean(MCRBs_OTFS_genie, 4);

crb_tau_ofdm_csi   = zeros(nSNR, P);
crb_tau_ofdm_ce    = zeros(nSNR, P);
crb_tau_ofdm_genie = zeros(nSNR, P);
crb_tau_otfs_csi   = zeros(nSNR, P);
crb_tau_otfs_ce    = zeros(nSNR, P);
crb_tau_otfs_genie = zeros(nSNR, P);
for p = 1:P
    crb_tau_ofdm_csi  (:, p) = sqrt(squeeze(MCRB_OFDM_mean      (4*p, 4*p, :)));
    crb_tau_ofdm_ce   (:, p) = sqrt(squeeze(MCRB_OFDM_CE_mean   (4*p, 4*p, :)));
    crb_tau_ofdm_genie(:, p) = sqrt(squeeze(MCRB_OFDM_genie_mean(4*p, 4*p, :)));
    crb_tau_otfs_csi  (:, p) = sqrt(squeeze(MCRB_OTFS_mean      (4*p, 4*p, :)));
    crb_tau_otfs_ce   (:, p) = sqrt(squeeze(MCRB_OTFS_CE_mean   (4*p, 4*p, :, 1)));
    crb_tau_otfs_genie(:, p) = sqrt(squeeze(MCRB_OTFS_genie_mean(4*p, 4*p, :)));
end

% Trial-mean MLE RMSE
rmse_tau_ofdm_csi   = sqrt(squeeze(mean(mse_tau_ofdm_csi,   2)));
rmse_tau_ofdm_ce    = sqrt(squeeze(mean(mse_tau_ofdm_ce,    2)));
rmse_tau_ofdm_genie = sqrt(squeeze(mean(mse_tau_ofdm_genie, 2)));
rmse_tau_otfs_csi   = sqrt(squeeze(mean(mse_tau_otfs_csi,   2)));
rmse_tau_otfs_genie = sqrt(squeeze(mean(mse_tau_otfs_genie, 2)));
rmse_tau_otfs_ce    = sqrt(squeeze(mean(mse_tau_otfs_ce,    2)));
if ndims(rmse_tau_otfs_ce) == 3
    rmse_tau_otfs_ce = rmse_tau_otfs_ce(:, :, 1);
end

legend_entries = {'MLE CSI', 'MLE CE', 'MLE genie', ...
                  'MCRB CSI', 'MCRB CE', 'CRB genie'};
colors = lines(P);

%% ---- OFDM per-path figures ----
for p = 1:P
    figure;
    semilogy(SNR_dbs, rmse_tau_ofdm_csi(:,p),   '-o',  'LineWidth', 1.8, 'Color', colors(1,:)); hold on;
    semilogy(SNR_dbs, rmse_tau_ofdm_ce(:,p),    '-s',  'LineWidth', 1.8, 'Color', colors(2,:));
    semilogy(SNR_dbs, rmse_tau_ofdm_genie(:,p), '-^',  'LineWidth', 1.8, 'Color', colors(3,:));
    semilogy(SNR_dbs, crb_tau_ofdm_csi(:,p),    '--',  'LineWidth', 1.4, 'Color', colors(1,:));
    semilogy(SNR_dbs, crb_tau_ofdm_ce(:,p),     ':',   'LineWidth', 1.4, 'Color', colors(2,:));
    semilogy(SNR_dbs, crb_tau_ofdm_genie(:,p),  '-.',  'LineWidth', 1.4, 'Color', colors(3,:));
    grid on;
    xlabel('SNR [dB]');
    ylabel('Delay RMSE [s]');
    title(sprintf('OFDM - Path %d: Delay MLE vs MCRB', p));
    legend(legend_entries, 'Location', 'southwest');
end

%% ---- OTFS per-path figures ----
for p = 1:P
    figure;
    semilogy(SNR_dbs, rmse_tau_otfs_csi(:,p),   '-o',  'LineWidth', 1.8, 'Color', colors(1,:)); hold on;
    semilogy(SNR_dbs, rmse_tau_otfs_ce(:,p),    '-s',  'LineWidth', 1.8, 'Color', colors(2,:));
    semilogy(SNR_dbs, rmse_tau_otfs_genie(:,p), '-^',  'LineWidth', 1.8, 'Color', colors(3,:));
    semilogy(SNR_dbs, crb_tau_otfs_csi(:,p),    '--',  'LineWidth', 1.4, 'Color', colors(1,:));
    semilogy(SNR_dbs, crb_tau_otfs_ce(:,p),     ':',   'LineWidth', 1.4, 'Color', colors(2,:));
    semilogy(SNR_dbs, crb_tau_otfs_genie(:,p),  '-.',  'LineWidth', 1.4, 'Color', colors(3,:));
    grid on;
    xlabel('SNR [dB]');
    ylabel('Delay RMSE [s]');
    title(sprintf('OTFS - Path %d: Delay MLE vs MCRB', p));
    legend(legend_entries, 'Location', 'southwest');
end
