% BER for EVA channel: OFDM (CSI), OTFS (CSI), OTFS (channel estimation)
% Channel estimation method: Raviteja et al., IEEE TVT 2019
%   "Embedded Pilot-Aided Channel Estimation for OTFS in Delay-Doppler Channels"
clear;
% close all;
% clc;

% len(SNR), nMC, time (s)
% 15      , 10  , 948

initialize_parameters;

% generate_indices;

bers_ofdm    = zeros(length(SNR_dbs), nMonteCarlo);
bers_otfs    = zeros(length(SNR_dbs), nMonteCarlo);
bers_otfs_ce = zeros(length(SNR_dbs), nMonteCarlo);

% --- Channel estimation (CE) design parameters ---
% Full guard symbols, fractional Doppler case
% (Raviteja 2019, Section III.B.1 -- corresponds to Fig. 9)
SNR_p_dB = 40;    % fixed pilot SNR (dB above noise floor)

% Pilot position (0-indexed)
l_p_ce = l_tau;          % delay  = 20  (= max delay tap)
k_p_ce = floor(M/2);     % Doppler = 64 (center of Doppler axis)

% Full guard: ALL M Doppler columns are zeroed at rows [l_p-l_tau, l_p+l_tau]
% Overhead = (2*l_tau+1)/N = 41/512 ≈ 8%  (matches paper Section III.B.1)
guard_row_lo = l_p_ce - l_tau;    % = 0
guard_row_hi = l_p_ce + l_tau;    % = 40

% Logical mask: true = data position, false = pilot/guard
data_mask_ce = true(N, M);
data_mask_ce(guard_row_lo+1:guard_row_hi+1, :) = false;  % all M cols are guard

total_duration = 0;

for SNR_idx = 1:length(SNR_dbs)
    SNR_db = SNR_dbs(SNR_idx);

    % Calculate noise variance N0 for a unit-power signal and unit-power channel
    SNR_lin = db2pow(SNR_db);
    N0 = 1 / SNR_lin;

    % Pilot amplitude scaled for fixed SNR_p = 40 dB above noise floor
    % |x_p|^2 / N0 = 10^(SNR_p_dB/10)  =>  x_p = sqrt(10^(SNR_p_dB/10) * N0)
    x_p_ce = sqrt(10^(SNR_p_dB/10) * N0);

    tic
    for mc_idx = 1:nMonteCarlo
        symbols = randi(M_order, N, M)-1;
        X = qammod(symbols, M_order, UnitAveragePower=true);
        X_TF = isfft(X);

        % Generate random phases and doppler shifts for THIS specific packet
        current_doppler = nu_max * cos(2*pi*rand(1, length(path_gains_amp)) - pi);
        current_path_gains = path_gains_amp .* exp(1j*rand(1, length(path_gains_amp))*2*pi);

        % Pass N0 as the noise variance, and use the randomized channel parameters
        Y = channel(X, current_path_gains, time_delays, current_doppler, deltaf, Ts, N0);
        Y_TF = channel(X_TF, current_path_gains, time_delays, current_doppler, deltaf, Ts, N0);

        X_hat = channel_inverse(Y, current_path_gains, time_delays, current_doppler, deltaf, Ts, N0);
        X_TF_hat = channel_inverse(Y_TF, current_path_gains, time_delays, current_doppler, deltaf, Ts, N0);
        X_DD_hat = sfft(X_TF_hat);

        transmitted_bits = qamdemod(X(:), M_order, OutputType="bit", UnitAveragePower=true);
        received_bits = qamdemod(X_hat(:), M_order, OutputType="bit", UnitAveragePower=true);
        bers_ofdm(SNR_idx, mc_idx) = mean(transmitted_bits ~= received_bits);

        received_bits = qamdemod(X_DD_hat(:), M_order, OutputType="bit", UnitAveragePower=true);
        bers_otfs(SNR_idx, mc_idx) = mean(transmitted_bits ~= received_bits);

        % ================================================================
        % OTFS with Embedded Pilot Channel Estimation
        % Full guard symbols, fractional Doppler (Fig. 9 -- Raviteja 2019 Sec. III.B.1)
        % ================================================================

        % 1. Build DD grid: data at non-guard rows, pilot at (l_p, k_p), zeros elsewhere
        X_DD_ce = zeros(N, M);
        X_DD_ce(data_mask_ce) = X(data_mask_ce);
        X_DD_ce(l_p_ce+1, k_p_ce+1) = x_p_ce;

        % 2. Transmit through channel
        Y_TF_ce = channel(isfft(X_DD_ce), current_path_gains, time_delays, ...
                          current_doppler, deltaf, Ts, N0);
        Y_DD_ce = sfft(Y_TF_ce);

        % 3. Estimate path parameters (tau, nu, alpha) from the pilot region.
        %    For each delay bin l' in [0, l_tau]:
        %      H_row[k''] = Y_DD[l_p+l', (k_p+k'') mod M] / x_p  for all k''
        %    With full guard, H_row is purely channel response + noise (no data ICI).
        %    Use oversampled FFT (4x) along the Doppler profile to estimate
        %    the fractional Doppler ν̂ at sub-bin resolution.
        %    Threshold per element: 3 * noise_std / x_p  (3σ criterion).
        K_os  = 4;                                 % Doppler oversampling factor
        est_alpha = zeros(1, P);
        est_tau   = zeros(1, P);
        est_nu    = zeros(1, P);

        % EVA delays are deterministic — iterate over the P known physical paths
        for i = 1:P
            l_prime_i = min(round(time_delays(i) * N * deltaf), l_tau);
            row   = l_p_ce + l_prime_i + 1;
            H_row = circshift(Y_DD_ce(row, :), -k_p_ce) / x_p_ce;

            % Find fractional Doppler via oversampled FFT of the Doppler profile
            H_row_os = fft(H_row, K_os * M);
            [~, k_os] = max(abs(H_row_os));
            k_frac = (k_os - 1) / K_os;          % fractional bin (0-indexed float)
            if k_frac >= M/2
                k_frac = k_frac - M;              % wrap to [-M/2, M/2)
            end

            % Gain: read H_row at integer bin nearest to estimated Doppler
            k_int_idx = mod(round(k_frac), M) + 1;

            est_tau(i)   = time_delays(i);        % exact known delay (not quantized)
            est_nu(i)    = k_frac * deltaf / M;
            est_alpha(i) = H_row(k_int_idx);
        end

        % 4. TF-domain MMSE equalization with estimated path parameters.
        X_TF_ce_hat = channel_inverse(Y_TF_ce, est_alpha, est_tau, est_nu, ...
                                      deltaf, Ts, N0);
        X_DD_ce_hat = sfft(X_TF_ce_hat);

        % 5. BER over data positions only
        tx_bits_ce = qamdemod(X(data_mask_ce), M_order, ...
                              OutputType="bit", UnitAveragePower=true);
        rx_bits_ce = qamdemod(X_DD_ce_hat(data_mask_ce), M_order, ...
                              OutputType="bit", UnitAveragePower=true);
        bers_otfs_ce(SNR_idx, mc_idx) = mean(tx_bits_ce ~= rx_bits_ce);
    end
    duration = toc;
    total_duration = total_duration + duration;
    remaining_time = (length(SNR_dbs) - SNR_idx)*duration;
    fprintf("%d/%d completed\n", SNR_idx, length(SNR_dbs));
    fprintf("Remaining Time: %f seconds\n", remaining_time);
    fprintf("----------------------------\n");
end

fprintf("Completed in %f seconds\n", total_duration);

bers_ofdm    = mean(bers_ofdm,    2);
bers_otfs    = mean(bers_otfs,    2);
bers_otfs_ce = mean(bers_otfs_ce, 2);

figure;
semilogy(SNR_dbs, bers_ofdm,    "LineStyle", '-',  'LineWidth', 2); hold on;
semilogy(SNR_dbs, bers_otfs,    "LineStyle", '-',  'LineWidth', 2);
semilogy(SNR_dbs, bers_otfs_ce, "LineStyle", '--', 'LineWidth', 2);
grid on;
xlabel("SNR [dB]");
ylabel("BER");
legend("OFDM (CSI)", "OTFS (CSI)", sprintf("OTFS (Ch. Est., SNR_p=%d dB)", SNR_p_dB));
theme(gcf, "light");
ylim([1e-5, 1]);
xlim([-10, 18]);
