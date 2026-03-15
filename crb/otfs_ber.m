clear;
% close all;
% clc;

type = "OTFS";

initialize_parameters;

generate_indices;

pilot_indices = zeros(N, M);
pilot_indices(35, 25) = 1;
pilot_indices = logical(pilot_indices);

guard_indices = zeros(N, M);
guard_indices(25:45, 15:35) = 1;
guard_indices(35, 25) = 0;
guard_indices = logical(guard_indices);

Xp = sqrt(sum(guard_indices(:)) + 1);

data_indices = guard_indices + pilot_indices;
data_indices = ~data_indices;

bers = zeros(length(SNR_dbs), nMonteCarlo);

tic
for SNR_idx = 1:length(SNR_dbs)
    SNR_db = SNR_dbs(SNR_idx);
    SNR_lin = SNR_lins(SNR_idx);
    % rng(0);
    parfor mc_idx = 1:nMonteCarlo
        alpha_gt = sqrt(SNR_lin * sigma2) * exp(1j*rand*2*pi);
        
        symbols = randi(M_order, N, M)-1;
        X_DD = qammod(symbols, M_order, UnitAveragePower=true);
        X_DD(guard_indices) = 0;
        X_DD(pilot_indices) = Xp;

        X_DDp = X_DD;
        X_DDp(data_indices) = 0;

        X_DDd = X_DD;
        X_DDd(pilot_indices) = 0;

        X_TF = isfft(X_DD);
        X_TFp = isfft(X_DDp);
        X_TFd = isfft(X_DDd);

        Y_TF = channel(X_TF, alpha_gt, tau_gt, nu_gt, deltaf, T, sigma2);
        y = Y_TF(:);

        Y_DD = sfft(Y_TF);

        %% Pilot Only Parameter Estimation
        % [tau_hat, nu_hat] = my_otfs_ce(Y_TF, deltaf, T);

        % eta0 = [nu_hat; tau_hat];   % initial guess
        eta0 = [nu_gt-5; tau_gt-0.5e-7];   % initial guess
        lb = [nu_gt-100000; tau_gt-1e-6];
        ub = [nu_gt+100000; tau_gt+1e-6];

        [eta_opt, Jmin] = fmincon( ...
            @(eta) cost_tau_nu(eta, deltaf, T, X_TFp, y), ...
            eta0, [], [], [], [], lb, ub, [], options);

        nu_hat  = eta_opt(1);
        tau_hat = eta_opt(2);

        alpha_hat = alpha_gt; % Remove this line!!!!!!!!!!!!!!!!!!

        %% Channel Estimation
        % H_TF_hat = alpha_hat * exp(-1j*2*pi*deltaf*ns*tau_hat) * exp(1j*2*pi*T*ms*nu_hat);

        %% X_hat Calculation
        % X_TF_hat = (Y_TF .* conj(H_TF_hat)) ./ (abs(H_TF_hat).^2 + sigma2);

        X_TF_hat = getChannelMatrixPinv(N, M, tau_hat, nu_hat, alpha_hat, deltaf, Ts) * y; % 10 1 100
        X_TF_hat = reshape(X_TF_hat, [N, M]);
        X_DD_hat = sfft(X_TF_hat);

        sym_hat = qamdemod(X_DD_hat, M_order, UnitAveragePower=true);
        X_DD_hat   = qammod(sym_hat, M_order, UnitAveragePower=true);
        X_DD_hat(guard_indices) = X_DD(guard_indices);
        X_DD_hat(pilot_indices) = X_DD(pilot_indices);
        X_TF_hat = isfft(X_DD_hat);

        transmitted_bits = qamdemod(X_DD(data_indices), M_order, OutputType="bit", UnitAveragePower=true);
        received_bits = qamdemod(X_DD_hat(data_indices), M_order, OutputType="bit", UnitAveragePower=true);
        bers(SNR_idx, mc_idx) = mean(transmitted_bits ~= received_bits);
        
    end
end
toc

bers = mean(bers, 2);

figure;
semilogy(SNR_dbs, bers, "LineStyle", '-', 'LineWidth', 2);
grid on;
xlabel("SNR [dB]");
ylabel("BER");
title("OTFS");
theme(gcf, "light");

