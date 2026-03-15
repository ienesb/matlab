clear;
% close all;
% clc;

type = "OFDM";

initialize_parameters;

generate_indices;

bers = zeros(length(SNR_dbs), nMonteCarlo);

for SNR_idx = 1:length(SNR_dbs)
    SNR_db = SNR_dbs(SNR_idx);
    SNR_lin = SNR_lins(SNR_idx);
    % rng(0);
    tic
    parfor mc_idx = 1:nMonteCarlo
        alpha_gt = sqrt(SNR_lin * sigma2) * exp(1j*rand*2*pi);
        
        symbols = randi(M_order, N, M)-1;
        X = qammod(symbols, M_order, UnitAveragePower=true);

        Y = channel(X, alpha_gt, tau_gt, nu_gt, deltaf, Ts, sigma2);
        y = Y(:);

        %% Pilot Only Parameter Estimation
        X_hat = zeros(N, M);
        
        X_hat(pilot_indices) = X(pilot_indices);

        eta0 = [nu_gt-5; tau_gt-0.5e-7];   % initial guess

        lb = [nu_gt-100000; tau_gt-1e-6];
        ub = [nu_gt+100000; tau_gt+1e-6];
        if SNR_db >= 20
            lb = [nu_gt-10; tau_gt-1e-7];
            ub = [nu_gt+10; tau_gt+1e-7];
        end
        [eta_opt, ~] = fmincon(@(eta) cost_tau_nu_pilots(eta, deltaf, Ts, X_hat, y, pilot_indices), ...
                  eta0, [], [], [], [], lb, ub, [], options);

        nu_hat  = eta_opt(1);
        tau_hat = eta_opt(2);
        
        % C = channel_noiseless(X_hat, 1, tau_hat, nu_hat, deltaf, Ts);
        % 
        % Cp = C(pilot_indices);
        % Yp = Y(pilot_indices);
        % 
        % alpha_hat = (Cp' * Yp) / (Cp' * Cp);
        alpha_hat = alpha_gt;
        % H_hat = alpha_hat * C;
        
        %% X_hat Calculation
        % X_hat = (Y .* conj(H_hat)) ./ (abs(H_hat).^2 + sigma2);
        X_hat = getChannelMatrixPinv(N, M, tau_hat, nu_hat, alpha_hat, deltaf, Ts) * y; % 10 1 100
        X_hat = reshape(X_hat, [N, M]);
       
        transmitted_bits = qamdemod(X(data_indices), M_order, OutputType="bit", UnitAveragePower=true);
        received_bits = qamdemod(X_hat(data_indices), M_order, OutputType="bit", UnitAveragePower=true);
        bers(SNR_idx, mc_idx) = mean(transmitted_bits ~= received_bits);
        
    end
    toc
    a = 0;
end

bers = mean(bers, 2);

figure;
semilogy(SNR_dbs, bers, "LineStyle", '-', 'LineWidth', 2);
grid on;
xlabel("SNR [dB]");
ylabel("BER");
title("OFDM");
theme(gcf, "light");

