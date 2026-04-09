clear;
% close all;
% clc;

type = "OFDM";

initialize_parameters;

% generate_indices;

% nMonteCarlo = 1;
MCRBs = zeros(4*P, 4*P, length(SNR_dbs), nMonteCarlo);

tic
for SNR_idx = 1:length(SNR_dbs)
    SNR_db = SNR_dbs(SNR_idx);

    % Calculate noise variance N0 for a unit-power signal and unit-power channel
    SNR_lin = db2pow(SNR_db);
    N0 = 1 / SNR_lin;
    for mc_idx = 1:nMonteCarlo     
        symbols = randi(M_order, N, M)-1;
        X = qammod(symbols, M_order, UnitAveragePower=true);
        X_TF = isfft(X);
    
        % Generate random phases and doppler shifts for THIS specific packet
        current_doppler = nu_max * cos(2*pi*rand(1, length(path_gains_amp)) - pi);
        current_path_gains = path_gains_amp .* exp(1j*rand(1, length(path_gains_amp))*2*pi);
    
        X_hat = X; % Assume no mismatch for now
    
        mu_gt = channel(X, current_path_gains, time_delays, current_doppler, deltaf, Ts, 0); % noiseless
        mu_gt = mu_gt(:);
    
        MCRBs(:, :, SNR_idx, mc_idx) = getMCRB_TF(X_hat, mu_gt, deltaf, Ts, N0, current_path_gains, current_doppler, time_delays, P);
    end
end
total_duration = toc;

fprintf("Completed in %f seconds\n", total_duration);

MCRBs = mean(MCRBs, 4);

colors = lines(P);

figure;
for p = 1:7
    MCRB_p = MCRBs(4*p, 4*p, :);
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

