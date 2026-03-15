%% CRB and MLE for pilot-only and full data OFDM (QPSK) (CRB Derived)
clear;
% close all;
% clc;

type = "OFDM";

initialize_parameters;

generate_indices;

CRBs_pilot_only = zeros(length(SNR_dbs), nMonteCarlo, 4, 4);
CRBs_pilot_only_ici = zeros(length(SNR_dbs), nMonteCarlo, 4, 4);

CRBs_full_data = zeros(length(SNR_dbs), nMonteCarlo, 4, 4);
CRBs_full_data_ici = zeros(length(SNR_dbs), nMonteCarlo, 4, 4);

tau_errors_pilot_only = zeros(length(SNR_dbs), nMonteCarlo);
tau_errors_full_data = zeros(length(SNR_dbs), nMonteCarlo);

nu_errors_pilot_only = zeros(length(SNR_dbs), nMonteCarlo);
nu_errors_full_data = zeros(length(SNR_dbs), nMonteCarlo);

tic
for SNR_idx = 1:length(SNR_dbs)
    SNR_db = SNR_dbs(SNR_idx);
    SNR_lin = SNR_lins(SNR_idx);
    % rng(0);
    parfor mc_idx = 1:nMonteCarlo
        alpha_gt = sqrt(SNR_lin * sigma2) * exp(1j*rand*2*pi);
        
        symbols = randi(M_order, N, M)-1;
        X = qammod(symbols, M_order, UnitAveragePower=true);
        
        ns = (0:(N-1)).';
        ms = 0:(M-1);
        
        % H = alpha_gt * exp(-1j*2*pi*deltaf*ns*tau_gt) * exp(1j*2*pi*T*ms*nu_gt);
        % 
        % Z = (randn(N, M) + 1j*randn(N, M)) * sqrt(sigma2/2);
        % % Z = 0;
        % Y = H .* X + Z;
        Y = channel(X, alpha_gt, tau_gt, nu_gt, deltaf, T, sigma2);
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
        [eta_opt, ~] = fmincon(@(eta) cost_tau_nu_pilots(eta, deltaf, T, X_hat, y, pilot_indices), ...
                  eta0, [], [], [], [], lb, ub, [], options);

        nu_hat  = eta_opt(1);
        tau_hat = eta_opt(2);

        tau_errors_pilot_only(SNR_idx, mc_idx) = abs(tau_gt - tau_hat);
        nu_errors_pilot_only(SNR_idx, mc_idx) = abs(nu_gt - nu_hat);

        CRBs_pilot_only(SNR_idx, mc_idx, :, :) = getCRB_OFDM(X, deltaf, T, sigma2, alpha_gt, nu_gt, tau_gt, pilot_indices);
        CRBs_pilot_only_ici(SNR_idx, mc_idx, :, :) = getCRB_OFDM_ICI(X, deltaf, T, sigma2, alpha_gt, nu_gt, tau_gt, pilot_indices);

        %% Full Data Parameter Estimation
        [eta_opt, ~] = fmincon( ...
            @(eta) cost_tau_nu(eta, deltaf, T, X, y), ...
            eta0, [], [], [], [], lb, ub, [], options);
        
        nu_opt = eta_opt(1);
        tau_opt = eta_opt(2);

        tau_errors_full_data(SNR_idx, mc_idx) = abs(tau_gt - tau_opt);
        nu_errors_full_data(SNR_idx, mc_idx) = abs(nu_gt - nu_opt);

        CRBs_full_data(SNR_idx, mc_idx, :, :) = getCRB_OFDM(X, deltaf, T, sigma2, alpha_gt, nu_gt, tau_gt, ones(N, M));
        CRBs_full_data_ici(SNR_idx, mc_idx, :, :) = getCRB_OFDM_ICI(X, deltaf, T, sigma2, alpha_gt, nu_gt, tau_gt, ones(N, M));

    end
end
toc

tau_errors_pilot_only = sqrt(mean(tau_errors_pilot_only.^2, 2));
range_errors_pilot_only = tau_errors_pilot_only .* c;

nu_errors_pilot_only = sqrt(mean(nu_errors_pilot_only.^2, 2));
velocity_errors_pilot_only = nu_errors_pilot_only .* (lambda / (2*cos(beta/2)));

tau_errors_full_data = sqrt(mean(tau_errors_full_data.^2, 2));
range_errors_full_data = tau_errors_full_data .* c;

nu_errors_full_data = sqrt(mean(nu_errors_full_data.^2, 2));
velocity_errors_full_data = nu_errors_full_data .* (lambda / (2*cos(beta/2)));

CRB_delay_pilot_onlyv2 = CRBs_pilot_only(:, :, 4, 4);
CRB_delay_pilot_onlyv2 = mean(CRB_delay_pilot_onlyv2, 2);
CRB_delay_pilot_onlyv2 = sqrt(CRB_delay_pilot_onlyv2);
CRB_range_pilot_onlyv2 = CRB_delay_pilot_onlyv2 .* c;

CRB_doppler_pilot_onlyv2 = CRBs_pilot_only(:, :, 3, 3);
CRB_doppler_pilot_onlyv2 = mean(CRB_doppler_pilot_onlyv2, 2);
CRB_doppler_pilot_onlyv2 = sqrt(CRB_doppler_pilot_onlyv2);
CRB_velocity_pilot_onlyv2 = CRB_doppler_pilot_onlyv2 .* (lambda / (2*cos(beta/2)));

CRB_delay_pilot_only_ici = CRBs_pilot_only_ici(:, :, 4, 4);
CRB_delay_pilot_only_ici = mean(CRB_delay_pilot_only_ici, 2);
CRB_delay_pilot_only_ici = sqrt(CRB_delay_pilot_only_ici);
CRB_range_pilot_only_ici = CRB_delay_pilot_only_ici .* c;

CRB_doppler_pilot_only_ici = CRBs_pilot_only_ici(:, :, 3, 3);
CRB_doppler_pilot_only_ici = mean(CRB_doppler_pilot_only_ici, 2);
CRB_doppler_pilot_only_ici = sqrt(CRB_doppler_pilot_only_ici);
CRB_velocity_pilot_only_ici = CRB_doppler_pilot_only_ici .* (lambda / (2*cos(beta/2)));


CRB_delay_full_datav2 = CRBs_full_data(:, :, 4, 4);
CRB_delay_full_datav2 = mean(CRB_delay_full_datav2, 2);
CRB_delay_full_datav2 = sqrt(CRB_delay_full_datav2);
CRB_range_full_datav2 = CRB_delay_full_datav2 .* c;

CRB_doppler_full_datav2 = CRBs_full_data(:, :, 3, 3);
CRB_doppler_full_datav2 = mean(CRB_doppler_full_datav2, 2);
CRB_doppler_full_datav2 = sqrt(CRB_doppler_full_datav2);
CRB_velocity_full_datav2 = CRB_doppler_full_datav2 .* (lambda / (2*cos(beta/2)));

CRB_delay_full_data_ici = CRBs_full_data_ici(:, :, 4, 4);
CRB_delay_full_data_ici = mean(CRB_delay_full_data_ici, 2);
CRB_delay_full_data_ici = sqrt(CRB_delay_full_data_ici);
CRB_range_full_data_ici = CRB_delay_full_data_ici .* c;

CRB_doppler_full_data_ici = CRBs_full_data_ici(:, :, 3, 3);
CRB_doppler_full_data_ici = mean(CRB_doppler_full_data_ici, 2);
CRB_doppler_full_data_ici = sqrt(CRB_doppler_full_data_ici);
CRB_velocity_full_data_ici = CRB_doppler_full_data_ici .* (lambda / (2*cos(beta/2)));

%% CRB Pilot Only
[CRB_delay_pilot_only, CRB_range_pilot_only, CRB_doppler_pilot_only, CRB_velocity_pilot_only] = getCRB_OFDM_derived(N, M, np, mp, c, lambda, deltaf, T, beta, SNR_lins);

%% CRB Full Data
[CRB_delay_full_data, CRB_range_full_data, CRB_doppler_full_data, CRB_velocity_full_data] = getCRB_OFDM_derived(N, M, 1, 1, c, lambda, deltaf, T, beta, SNR_lins);

%% Figure

colors = lines(7);

figure;
semilogy(SNR_dbs, tau_errors_pilot_only, "Color", colors(2, :), "LineStyle", '-', 'LineWidth', 2); hold on;
semilogy(SNR_dbs, CRB_delay_pilot_only, "Color", colors(2, :), "LineStyle", '--', 'LineWidth', 2);
semilogy(SNR_dbs, tau_errors_full_data, "Color", colors(3, :), "LineStyle", '-', 'LineWidth', 2);
semilogy(SNR_dbs, CRB_delay_full_data, "Color", colors(3, :), "LineStyle", '--', 'LineWidth', 2);
semilogy(SNR_dbs, CRB_delay_full_datav2, "Color", colors(4, :), "LineStyle", '--', 'LineWidth', 2);
semilogy(SNR_dbs, CRB_delay_full_data_ici, "Color", colors(5, :), "LineStyle", '--', 'LineWidth', 2);
semilogy(SNR_dbs, CRB_delay_pilot_onlyv2, "Color", colors(6, :), "LineStyle", '--', 'LineWidth', 2);
semilogy(SNR_dbs, CRB_delay_pilot_only_ici, "Color", colors(7, :), "LineStyle", '--', 'LineWidth', 2);
grid on;
legend("Pilot Only MLE", "Pilot Only CRB", "Full Data MLE", "Full Data CRB", "Full Data CRB v2", "Full Data CRB ICI", "Pilot Only CRB v2", "Pilot Only CRB ICI");
xlabel("SNR [dB]");
ylabel("Delay RMSE [s]");
title("OFDM");
theme(gcf, "light");

figure;
semilogy(SNR_dbs, nu_errors_pilot_only, "Color", colors(2, :), "LineStyle", '-', 'LineWidth', 2); hold on;
semilogy(SNR_dbs, CRB_doppler_pilot_only, "Color", colors(2, :), "LineStyle", '--', 'LineWidth', 2);
semilogy(SNR_dbs, nu_errors_full_data, "Color", colors(3, :), "LineStyle", '-', 'LineWidth', 2);
semilogy(SNR_dbs, CRB_doppler_full_data, "Color", colors(3, :), "LineStyle", '--', 'LineWidth', 2);
semilogy(SNR_dbs, CRB_doppler_full_datav2, "Color", colors(4, :), "LineStyle", '--', 'LineWidth', 2);
semilogy(SNR_dbs, CRB_doppler_full_data_ici, "Color", colors(5, :), "LineStyle", '--', 'LineWidth', 2);
semilogy(SNR_dbs, CRB_doppler_pilot_onlyv2, "Color", colors(6, :), "LineStyle", '--', 'LineWidth', 2);
semilogy(SNR_dbs, CRB_doppler_pilot_only_ici, "Color", colors(7, :), "LineStyle", '--', 'LineWidth', 2);
grid on;
legend("Pilot Only MLE", "Pilot Only CRB", "Full Data MLE", "Full Data CRB", "Full Data CRB v2", "Full Data CRB ICI", "Pilot Only CRB v2", "Pilot Only CRB ICI");
xlabel("SNR [dB]");
ylabel("Doppler RMSE [Hz]");
title("OFDM");
theme(gcf, "light");

figure;
semilogy(SNR_dbs, range_errors_pilot_only, "Color", colors(2, :), "LineStyle", '-', 'LineWidth', 2); hold on;
semilogy(SNR_dbs, CRB_range_pilot_only, "Color", colors(2, :), "LineStyle", '--', 'LineWidth', 2);
semilogy(SNR_dbs, range_errors_full_data, "Color", colors(3, :), "LineStyle", '-', 'LineWidth', 2);
semilogy(SNR_dbs, CRB_range_full_data, "Color", colors(3, :), "LineStyle", '--', 'LineWidth', 2);
semilogy(SNR_dbs, CRB_range_full_datav2, "Color", colors(4, :), "LineStyle", '--', 'LineWidth', 2);
semilogy(SNR_dbs, CRB_range_full_data_ici, "Color", colors(5, :), "LineStyle", '--', 'LineWidth', 2);
semilogy(SNR_dbs, CRB_range_pilot_onlyv2, "Color", colors(6, :), "LineStyle", '--', 'LineWidth', 2);
semilogy(SNR_dbs, CRB_range_pilot_only_ici, "Color", colors(7, :), "LineStyle", '--', 'LineWidth', 2);
grid on;
legend("Pilot Only MLE", "Pilot Only CRB", "Full Data MLE", "Full Data CRB", "Full Data CRB v2", "Full Data CRB ICI", "Pilot Only CRB v2", "Pilot Only CRB ICI");
xlabel("SNR [dB]");
ylabel("Range RMSE [m]");
title("OFDM");
theme(gcf, "light");

figure;
semilogy(SNR_dbs, velocity_errors_pilot_only, "Color", colors(2, :), "LineStyle", '-', 'LineWidth', 2); hold on;
semilogy(SNR_dbs, CRB_velocity_pilot_only, "Color", colors(2, :), "LineStyle", '--', 'LineWidth', 2);
semilogy(SNR_dbs, velocity_errors_full_data, "Color", colors(3, :), "LineStyle", '-', 'LineWidth', 2);
semilogy(SNR_dbs, CRB_velocity_full_data, "Color", colors(3, :), "LineStyle", '--', 'LineWidth', 2);
semilogy(SNR_dbs, CRB_velocity_full_datav2, "Color", colors(4, :), "LineStyle", '--', 'LineWidth', 2);
semilogy(SNR_dbs, CRB_velocity_full_data_ici, "Color", colors(5, :), "LineStyle", '--', 'LineWidth', 2);
semilogy(SNR_dbs, CRB_velocity_pilot_onlyv2, "Color", colors(6, :), "LineStyle", '--', 'LineWidth', 2);
semilogy(SNR_dbs, CRB_velocity_pilot_only_ici, "Color", colors(7, :), "LineStyle", '--', 'LineWidth', 2);
grid on;
legend("Pilot Only MLE", "Pilot Only CRB", "Full Data MLE", "Full Data CRB", "Full Data CRB v2", "Full Data CRB ICI", "Pilot Only CRB v2", "Pilot Only CRB ICI");
xlabel("SNR [dB]");
ylabel("Velocity RMSE [m/s]");
title("OFDM");
theme(gcf, "light");