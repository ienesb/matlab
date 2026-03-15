% BER for EVA channel for OFDM/OTFS (full CSI)
clear;
% close all;
% clc;

% len(SNR), nMC, time (s)
% 15      , 10  , 948

initialize_parameters;

% generate_indices;

bers_ofdm = zeros(length(SNR_dbs), nMonteCarlo);
bers_otfs = zeros(length(SNR_dbs), nMonteCarlo);

total_duration = 0;

for SNR_idx = 1:length(SNR_dbs)
    SNR_db = SNR_dbs(SNR_idx);
    
    % Calculate noise variance N0 for a unit-power signal and unit-power channel
    SNR_lin = db2pow(SNR_db);
    N0 = 1 / SNR_lin;

    tic
    parfor mc_idx = 1:nMonteCarlo     
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
    end
    duration = toc;
    total_duration = total_duration + duration;
    remaining_time = (length(SNR_dbs) - SNR_idx)*duration;
    fprintf("%d/%d completed\n", SNR_idx, length(SNR_dbs));
    fprintf("Remaining Time: %f seconds\n", remaining_time);
    fprintf("----------------------------\n");
end

fprintf("Completed in %f seconds\n", total_duration);

bers_ofdm = mean(bers_ofdm, 2);
bers_otfs = mean(bers_otfs, 2);

figure;
semilogy(SNR_dbs, bers_ofdm, "LineStyle", '-', 'LineWidth', 2); hold on;
semilogy(SNR_dbs, bers_otfs, "LineStyle", '-', 'LineWidth', 2);
grid on;
xlabel("SNR [dB]");
ylabel("BER");
legend("OFDM", "OTFS");
theme(gcf, "light");
ylim([1e-5, 1]);
xlim([-10, 18]);
