% CRB for EVA channel for OFDM (Pilot Only)
clear;
% close all;
% clc;

% nMC, time (s)
% 100 , 124

type = "OFDM";

initialize_parameters;

generate_indices;

FIMs_full_data_ici = zeros(nMonteCarlo, P, P);
FIMs_pilot_only_ici = zeros(nMonteCarlo, P, P);
FIMs_OTFS = zeros(nMonteCarlo, P, P);

N0 = 1 ./ SNR_lins;

tic
parfor mc_idx = 1:nMonteCarlo     
    symbols = randi(M_order, N, M)-1;
    X = qammod(symbols, M_order, UnitAveragePower=true);
    X_TF = isfft(X);

    % Generate random phases and doppler shifts for THIS specific packet
    current_doppler = nu_max * cos(2*pi*rand(1, length(path_gains_amp)) - pi);
    current_path_gains = path_gains_amp .* exp(1j*rand(1, length(path_gains_amp))*2*pi);

    FIMs_pilot_only_ici(mc_idx, :, :) = getFIM_OFDM_ICI(X, deltaf, T, current_path_gains, current_doppler, time_delays, pilot_indices);
    FIMs_full_data_ici(mc_idx, :, :) = getFIM_OFDM_ICI(X, deltaf, T, current_path_gains, current_doppler, time_delays, ones(N, M));
    FIMs_OTFS(mc_idx, :, :) = getFIM_OFDM_ICI(X_TF, deltaf, T, current_path_gains, current_doppler, time_delays, ones(N, M));
end
total_duration = toc;

fprintf("Completed in %f seconds\n", total_duration);


FIM_pilot_only = squeeze(mean(FIMs_pilot_only_ici, 1));
CRB_pilot_only = inv(FIM_pilot_only);

FIM_full_data = squeeze(mean(FIMs_full_data_ici, 1));
CRB_full_data = inv(FIM_full_data);

FIM_OTFS = squeeze(mean(FIMs_OTFS, 1));
CRB_OTFS = inv(FIM_OTFS);

colors = lines(P);

figure;
for p = 1:7
    CRB_delay_pilot_only_ici = CRB_pilot_only(p, p);
    CRB_delay_pilot_only_ici = sqrt(CRB_delay_pilot_only_ici);
    CRB_delay_pilot_only_ici = CRB_delay_pilot_only_ici * sqrt(N0/2);

    semilogy(SNR_dbs, CRB_delay_pilot_only_ici,"Color", colors(p, :), "LineStyle", '--', 'LineWidth', 2); hold on;
end
legend("1", "2", "3", "4 " ,"5", "6", "7", "8", "9");
grid on;
xlabel("SNR [dB]");
ylabel("Delay RMSE [s]");
title("OFDM Pilot Only");
theme(gcf, "light");

figure;
for p = 1:7
    CRB_delay_full_data_ici = CRB_full_data(p, p);
    CRB_delay_full_data_ici = sqrt(CRB_delay_full_data_ici);
    CRB_delay_full_data_ici = CRB_delay_full_data_ici * sqrt(N0/2);

    semilogy(SNR_dbs, CRB_delay_full_data_ici,"Color", colors(p, :), "LineStyle", '--', 'LineWidth', 2); hold on;
end
legend("1", "2", "3", "4 " ,"5", "6", "7", "8", "9");
grid on;
xlabel("SNR [dB]");
ylabel("Delay RMSE [s]");
title("OFDM Full Data");
theme(gcf, "light");

figure;
for p = 1:7
    CRB_delay_OTFS = CRB_full_data(p, p);
    CRB_delay_OTFS = sqrt(CRB_delay_OTFS);
    CRB_delay_OTFS = CRB_delay_OTFS * sqrt(N0/2);

    semilogy(SNR_dbs, CRB_delay_OTFS,"Color", colors(p, :), "LineStyle", '--', 'LineWidth', 2); hold on;
end
legend("1", "2", "3", "4 " ,"5", "6", "7", "8", "9");
grid on;
xlabel("SNR [dB]");
ylabel("Delay RMSE [s]");
title("OTFS");
theme(gcf, "light");
