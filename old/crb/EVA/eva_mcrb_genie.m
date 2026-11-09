% eva_mcrb_genie.m
% Genie CRB for EVA channel: OFDM and OTFS with true data symbols.
clear;

initialize_parameters;

% ---- Pilot design (for OTFS signal construction) ----
SNR_p_dBs = [45];
n_pilot   = length(SNR_p_dBs);

l_p = l_tau;
k_p = floor(M/2);

guard_row_lo = l_p - l_tau;
guard_row_hi = l_p + l_tau;

data_mask = true(N, M);
data_mask(guard_row_lo+1 : guard_row_hi+1, :) = false;

% Pre-allocate
nSNR = length(SNR_dbs);
MCRBs_OFDM_genie = zeros(4*P, 4*P, nSNR, nMonteCarlo);
MCRBs_OTFS_genie = zeros(4*P, 4*P, nSNR, nMonteCarlo);

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

        %% ======== OFDM Genie ========
        MCRBs_OFDM_genie(:, :, SNR_idx, mc_idx) = getMCRB_TF_genie(X, ...
            deltaf, Ts, N0, current_path_gains, current_doppler, time_delays, P);

        %% ======== OTFS Genie ========
        x_p_max = sqrt(10^(max(SNR_p_dBs)/10) * N0);

        X_DD_data = zeros(N, M);
        X_DD_data(data_mask) = X(data_mask);

        X_impulse = zeros(N, M);
        X_impulse(l_p+1, k_p+1) = 1;

        X_csi_TF = isfft(X_DD_data + x_p_max*X_impulse);

        MCRBs_OTFS_genie(:, :, SNR_idx, mc_idx) = getMCRB_TF_genie(X_csi_TF, ...
            deltaf, Ts, N0, current_path_gains, current_doppler, time_delays, P);
    end
    duration       = toc;
    total_duration = total_duration + duration;
    remaining_time = (nSNR - SNR_idx) * duration;
    fprintf("%d/%d | %.1fs (remain: %.1fs)\n", SNR_idx, nSNR, duration, remaining_time);
end

fprintf("Completed in %.1f seconds\n", total_duration);
