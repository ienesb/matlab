% eva_mcrb_compare_tf_dd.m
% Compare OTFS (getMCRB_DD) vs OFDM (getMCRB_TF) MCRB at fixed bit-error rates.
% No detection or channel estimation — X_hat is built by flipping a controlled
% fraction of the data bits, so the BER is set exactly.
%
% Saves per-realization A, B, bias_mat for both schemes so the variance
% (A^{-1} B A^{-1}) and bias-squared components can be inspected post-hoc.
clear;

initialize_parameters;

% ---- Fixed BER points ----
% target_bers = [0, 0.01, 0.1];
target_bers = linspace(0, 4e-2, 20);
n_ber = length(target_bers);

% ---- OTFS pilot / guard layout (matches eva_ber_mp_mcrb_dd.m) ----
SNR_p_dB = 45;
l_p = l_tau;
k_p = floor(M/2);
guard_row_lo = l_p - l_tau;
guard_row_hi = l_p + l_tau;

data_mask = true(N, M);
data_mask(guard_row_lo+1 : guard_row_hi+1, :) = false;

X_impulse = zeros(N, M);
X_impulse(l_p+1, k_p+1) = 1;

bits_per_sym = log2(M_order);
n_bits_otfs  = nnz(data_mask) * bits_per_sym;   % OTFS: data positions only
n_bits_ofdm  = N * M           * bits_per_sym;  % OFDM: full TF grid is data

% ---- Pre-allocate ----
nSNR    = length(SNR_dbs);
shape5  = [4*P, 4*P, nSNR, nMonteCarlo, n_ber];

MCRBs_OTFS  = zeros(shape5);
As_OTFS     = zeros(shape5);
Bs_OTFS     = zeros(shape5);
biases_OTFS = zeros(shape5);

MCRBs_OFDM  = zeros(shape5);
As_OFDM     = zeros(shape5);
Bs_OFDM     = zeros(shape5);
biases_OFDM = zeros(shape5);

total_duration = 0;
rng(2028);

for SNR_idx = 1:nSNR
    SNR_db = SNR_dbs(SNR_idx);
    N0     = 1 / db2pow(SNR_db);
    x_p    = sqrt(10^(SNR_p_dB/10) * N0);

    tic
    parfor mc_idx = 1:nMonteCarlo
        % ---- Channel realization ----
        current_doppler    = nu_max * cos(2*pi*rand(1, P) - pi);
        current_path_gains = path_gains_amp .* exp(1j*rand(1, P)*2*pi);

        % ---- Random data ----
        symbols = randi(M_order, N, M) - 1;
        X       = qammod(symbols, M_order, UnitAveragePower=true);

        % ---- OTFS true symbols (DD, with embedded pilot) ----
        X_DD_data = zeros(N, M);
        X_DD_data(data_mask) = X(data_mask);
        X_true_DD = X_DD_data + x_p * X_impulse;

        true_bits_otfs = qamdemod(X(data_mask), M_order, ...
                                  OutputType="bit", UnitAveragePower=true);

        % ---- OFDM true symbols (TF, full grid is data, no pilot) ----
        X_true_TF      = X;
        true_bits_ofdm = qamdemod(X(:), M_order, ...
                                  OutputType="bit", UnitAveragePower=true);

        % Local accumulators (one slice per BER)
        mcrb_otfs_l = zeros(4*P, 4*P, n_ber);  A_otfs_l = zeros(4*P, 4*P, n_ber);
        B_otfs_l    = zeros(4*P, 4*P, n_ber);  bias_otfs_l = zeros(4*P, 4*P, n_ber);

        mcrb_ofdm_l = zeros(4*P, 4*P, n_ber);  A_ofdm_l = zeros(4*P, 4*P, n_ber);
        B_ofdm_l    = zeros(4*P, 4*P, n_ber);  bias_ofdm_l = zeros(4*P, 4*P, n_ber);

        for b_idx = 1:n_ber
            ber = target_bers(b_idx);

            % ---------- OTFS: flip OTFS-data bits, keep pilot ----------
            n_flips = round(ber * n_bits_otfs);
            flipped = true_bits_otfs;
            if n_flips > 0
                idx = randperm(n_bits_otfs, n_flips);
                flipped(idx) = ~flipped(idx);
            end
            X_hat_data = qammod(flipped, M_order, ...
                                InputType="bit", UnitAveragePower=true);
            X_hat_DD = zeros(N, M);
            X_hat_DD(data_mask) = X_hat_data;
            X_hat_DD = X_hat_DD + x_p * X_impulse;

            [mcrb_o, A_o, B_o, bias_o] = getMCRB_DD( ...
                X_hat_DD, X_true_DD, deltaf, Ts, N0, ...
                current_path_gains, current_doppler, time_delays, P);

            mcrb_otfs_l(:, :, b_idx) = mcrb_o;
            A_otfs_l   (:, :, b_idx) = A_o;
            B_otfs_l   (:, :, b_idx) = B_o;
            bias_otfs_l(:, :, b_idx) = bias_o;

            % ---------- OFDM: flip bits over full TF grid ----------
            n_flips = round(ber * n_bits_ofdm);
            flipped = true_bits_ofdm;
            if n_flips > 0
                idx = randperm(n_bits_ofdm, n_flips);
                flipped(idx) = ~flipped(idx);
            end
            X_hat_TF = reshape( ...
                qammod(flipped, M_order, InputType="bit", UnitAveragePower=true), ...
                [N, M]);

            [mcrb_f, A_f, B_f, bias_f] = getMCRB_TF( ...
                X_hat_TF, X_true_TF, deltaf, Ts, N0, ...
                current_path_gains, current_doppler, time_delays, P);

            mcrb_ofdm_l(:, :, b_idx) = mcrb_f;
            A_ofdm_l   (:, :, b_idx) = A_f;
            B_ofdm_l   (:, :, b_idx) = B_f;
            bias_ofdm_l(:, :, b_idx) = bias_f;
        end

        MCRBs_OTFS (:, :, SNR_idx, mc_idx, :) = mcrb_otfs_l;
        As_OTFS    (:, :, SNR_idx, mc_idx, :) = A_otfs_l;
        Bs_OTFS    (:, :, SNR_idx, mc_idx, :) = B_otfs_l;
        biases_OTFS(:, :, SNR_idx, mc_idx, :) = bias_otfs_l;

        MCRBs_OFDM (:, :, SNR_idx, mc_idx, :) = mcrb_ofdm_l;
        As_OFDM    (:, :, SNR_idx, mc_idx, :) = A_ofdm_l;
        Bs_OFDM    (:, :, SNR_idx, mc_idx, :) = B_ofdm_l;
        biases_OFDM(:, :, SNR_idx, mc_idx, :) = bias_ofdm_l;
    end
    duration       = toc;
    total_duration = total_duration + duration;
    remaining      = (nSNR - SNR_idx) * duration;
    fprintf("%d/%d | %.1fs (remain: %.1fs)\n", SNR_idx, nSNR, duration, remaining);
    save eva_mcrb_compare_tf_dd_results.mat
end

fprintf("Completed in %.1f seconds\n", total_duration);
save eva_mcrb_fixed_ber_resultsv5.mat
