%% MCRB and MLE for pilot-only, data aided, and genie OFDM (QPSK) (CRB derived)
clear;
% close all;
clc;

type = "OFDM";

initialize_parameters;

generate_indices;

MCRBs = zeros(4, 4, length(SNR_dbs), nMonteCarlo);
tau_errors = zeros(length(SNR_dbs), nMonteCarlo);
tau_errors_pilot_only = zeros(length(SNR_dbs), nMonteCarlo);
tau_errors_full_data = zeros(length(SNR_dbs), nMonteCarlo);

nu_errors = zeros(length(SNR_dbs), nMonteCarlo);
nu_errors_pilot_only = zeros(length(SNR_dbs), nMonteCarlo);
nu_errors_full_data = zeros(length(SNR_dbs), nMonteCarlo);

% pilot_periods = [10, 2, 2, 1; 5, 5, 1, 1];
np = 2;
mp = 5;

pilot_indices = zeros(N, M);
pilot_indices(1:np:end, 1:mp:end) = 1;
pilot_indices = logical(pilot_indices);

data_indices = ~pilot_indices;

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
    
        Y = channel(X, alpha_gt, tau_gt, nu_gt, deltaf, Ts, sigma2);
        y = Y(:);

        %% Channel Estimation
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

        tau_errors_pilot_only(SNR_idx, mc_idx) = abs(tau_gt - tau_hat);
        nu_errors_pilot_only(SNR_idx, mc_idx) = abs(nu_gt - nu_hat);
        
        % E = exp(-1j*2*pi*deltaf*ns*tau_hat) * exp(1j*2*pi*Ts*ms*nu_hat);
        % 
        % C = (E .* X_hat);        % use true pilots (X on pilots)
        % 
        % Cp = C(pilot_indices);
        % Yp = Y(pilot_indices);
        % 
        % alpha_hat = (Cp' * Yp) / (Cp' * Cp);
        alpha_hat = alpha_gt;
        
        % H_hat = alpha_hat * E;     % full grid channel (phase terms apply everywhere)


        %% Full Data Parameter Estimation
        [eta_opt, ~] = fmincon( ...
            @(eta) cost_tau_nu(eta, deltaf, Ts, X, y), ...
            eta0, [], [], [], [], lb, ub, [], options);
        
        nu_opt = eta_opt(1);
        tau_opt = eta_opt(2);

        tau_errors_full_data(SNR_idx, mc_idx) = abs(tau_gt - tau_opt);
        nu_errors_full_data(SNR_idx, mc_idx) = abs(nu_gt - nu_opt);

        %% Data-Aided Parameter Estimation
        % X_hat = (Y .* conj(H_hat)) ./ (abs(H_hat).^2 + sigma2);
        % X_soft = (Y .* conj(H_hat)) ./ (abs(H_hat).^2 + sigma2);
        X_soft = getChannelMatrixPinv(N, M, tau_hat, nu_hat, alpha_hat, deltaf, Ts) * y; % 10 1 100
        X_soft = reshape(X_soft, [N, M]);
        sym_hat = qamdemod(X_soft, M_order, UnitAveragePower=true);
        X_hat   = qammod(sym_hat, M_order, UnitAveragePower=true);
        
        X_hat(pilot_indices) = X(pilot_indices);

        [eta_opt, Jmin] = fmincon( ...
            @(eta) cost_tau_nu(eta, deltaf, Ts, X_hat, y), ...
            eta0, [], [], [], [], lb, ub, [], options);
        
        nu_opt = eta_opt(1);
        tau_opt = eta_opt(2);

        tau_errors(SNR_idx, mc_idx) = abs(tau_gt - tau_opt);
        nu_errors(SNR_idx, mc_idx) = abs(nu_gt - nu_opt);

        %% MCRB Calculation
        
        MCRB = getMCRB_OFDM(X_hat, mu_gt, deltaf, Ts, sigma2, alpha_gt, nu_gt, tau_gt);
        MCRBs(:, :, SNR_idx, mc_idx) = MCRB;
    end
end
toc

MCRBs = mean(MCRBs, 4);

MCRB_range = sqrt(MCRBs(4, 4, :) * c^2);
MCRB_range = MCRB_range(:);
MCRB_delay = MCRB_range ./ c;

MCRB_velocity = sqrt(MCRBs(3, 3, :));
MCRB_velocity = MCRB_velocity(:) * (lambda / (2*cos(beta/2)));
MCRB_doppler = MCRB_velocity ./ (lambda / (2*cos(beta/2)));

tau_errors = sqrt(mean(tau_errors.^2, 2));
range_errors = tau_errors .* c;

nu_errors = sqrt(mean(nu_errors.^2, 2));
velocity_errors = nu_errors .* (lambda / (2*cos(beta/2)));

tau_errors_pilot_only = sqrt(mean(tau_errors_pilot_only.^2, 2));
range_errors_pilot_only = tau_errors_pilot_only .* c;

nu_errors_pilot_only = sqrt(mean(nu_errors_pilot_only.^2, 2));
velocity_errors_pilot_only = nu_errors_pilot_only .* (lambda / (2*cos(beta/2)));

tau_errors_full_data = sqrt(mean(tau_errors_full_data.^2, 2));
range_errors_full_data = tau_errors_full_data .* c;

nu_errors_full_data = sqrt(mean(nu_errors_full_data.^2, 2));
velocity_errors_full_data = nu_errors_full_data .* (lambda / (2*cos(beta/2)));


%% CRB Pilot Only
K = ceil((N - 1) / np);
L = ceil((M - 1) / mp);

absP = (K + 1) * (L + 1);

CRB_range_pilot_only = 12 / (K * (K + 2) * absP * np^2) * c^2 ./ (8 * pi^2 * SNR_lins * deltaf^2);
CRB_range_pilot_only = sqrt(CRB_range_pilot_only);
CRB_delay_pilot_only = CRB_range_pilot_only ./ c;

CRB_velocity_pilot_only = 12 / (L * (L + 2) * absP * mp^2) * lambda^2 ./ (32 * pi^2 * SNR_lins * Ts^2 * cos(beta/2)^2);
CRB_velocity_pilot_only = sqrt(CRB_velocity_pilot_only);
CRB_doppler_pilot_only = CRB_velocity_pilot_only ./ (lambda / (2*cos(beta/2)));

%% CRB Full Data
np = 1;
mp = 1;
K = ceil((N - 1) / np);
L = ceil((M - 1) / mp);

absP = (K + 1) * (L + 1);

CRB_range_full_data = 12 / (K * (K + 2) * absP * np^2) * c^2 ./ (8 * pi^2 * SNR_lins * deltaf^2);
CRB_range_full_data = sqrt(CRB_range_full_data);
CRB_delay_full_data = CRB_range_full_data ./ c;

CRB_velocity_full_data = 12 / (L * (L + 2) * absP * mp^2) * lambda^2 ./ (32 * pi^2 * SNR_lins * Ts^2 * cos(beta/2)^2);
CRB_velocity_full_data = sqrt(CRB_velocity_full_data);
CRB_doppler_full_data = CRB_velocity_full_data ./ (lambda / (2*cos(beta/2)));

%% Figure

colors = lines(3);

figure;
semilogy(SNR_dbs, tau_errors, "Color", colors(1, :), "LineStyle", '-', 'LineWidth', 2); hold on;
semilogy(SNR_dbs, MCRB_delay, "Color", colors(1, :), "LineStyle", '--', 'LineWidth', 2);
semilogy(SNR_dbs, tau_errors_pilot_only, "Color", colors(2, :), "LineStyle", '-', 'LineWidth', 2);
semilogy(SNR_dbs, CRB_delay_pilot_only, "Color", colors(2, :), "LineStyle", '--', 'LineWidth', 2);
semilogy(SNR_dbs, tau_errors_full_data, "Color", colors(3, :), "LineStyle", '-', 'LineWidth', 2);
semilogy(SNR_dbs, CRB_delay_full_data, "Color", colors(3, :), "LineStyle", '--', 'LineWidth', 2);
grid on;
legend("Data Aided MLE", "Data Aided MCRB", "Pilot Only MLE", "Pilot Only CRB", "Full Data MLE", "Full Data CRB");
xlabel("SNR [dB]");
ylabel("Delay RMSE [s]");
title("OFDM");
theme(gcf, "light");

figure;
semilogy(SNR_dbs, nu_errors, "Color", colors(1, :), "LineStyle", '-', 'LineWidth', 2); hold on;
semilogy(SNR_dbs, MCRB_doppler, "Color", colors(1, :), "LineStyle", '--', 'LineWidth', 2);
semilogy(SNR_dbs, nu_errors_pilot_only, "Color", colors(2, :), "LineStyle", '-', 'LineWidth', 2);
semilogy(SNR_dbs, CRB_doppler_pilot_only, "Color", colors(2, :), "LineStyle", '--', 'LineWidth', 2);
semilogy(SNR_dbs, nu_errors_full_data, "Color", colors(3, :), "LineStyle", '-', 'LineWidth', 2);
semilogy(SNR_dbs, CRB_doppler_full_data, "Color", colors(3, :), "LineStyle", '--', 'LineWidth', 2);
grid on;
legend("Data Aided MLE", "Data Aided MCRB", "Pilot Only MLE", "Pilot Only CRB", "Full Data MLE", "Full Data CRB");
xlabel("SNR [dB]");
ylabel("Doppler RMSE [Hz]");
title("OFDM");
theme(gcf, "light");

figure;
semilogy(SNR_dbs, range_errors, "Color", colors(1, :), "LineStyle", '-', 'LineWidth', 2); hold on;
semilogy(SNR_dbs, MCRB_range, "Color", colors(1, :), "LineStyle", '--', 'LineWidth', 2);
semilogy(SNR_dbs, range_errors_pilot_only, "Color", colors(2, :), "LineStyle", '-', 'LineWidth', 2);
semilogy(SNR_dbs, CRB_range_pilot_only, "Color", colors(2, :), "LineStyle", '--', 'LineWidth', 2);
semilogy(SNR_dbs, range_errors_full_data, "Color", colors(3, :), "LineStyle", '-', 'LineWidth', 2);
semilogy(SNR_dbs, CRB_range_full_data, "Color", colors(3, :), "LineStyle", '--', 'LineWidth', 2);
grid on;
legend("Data Aided MLE", "Data Aided MCRB", "Pilot Only MLE", "Pilot Only CRB", "Full Data MLE", "Full Data CRB");
xlabel("SNR [dB]");
ylabel("Range RMSE [m]");
title("OFDM");
theme(gcf, "light");

figure;
semilogy(SNR_dbs, velocity_errors, "Color", colors(1, :), "LineStyle", '-', 'LineWidth', 2); hold on;
semilogy(SNR_dbs, MCRB_velocity, "Color", colors(1, :), "LineStyle", '--', 'LineWidth', 2);
semilogy(SNR_dbs, velocity_errors_pilot_only, "Color", colors(2, :), "LineStyle", '-', 'LineWidth', 2);
semilogy(SNR_dbs, CRB_velocity_pilot_only, "Color", colors(2, :), "LineStyle", '--', 'LineWidth', 2);
semilogy(SNR_dbs, velocity_errors_full_data, "Color", colors(3, :), "LineStyle", '-', 'LineWidth', 2);
semilogy(SNR_dbs, CRB_velocity_full_data, "Color", colors(3, :), "LineStyle", '--', 'LineWidth', 2);
grid on;
legend("Data Aided MLE", "Data Aided MCRB", "Pilot Only MLE", "Pilot Only CRB", "Full Data MLE", "Full Data CRB");
xlabel("SNR [dB]");
ylabel("Velocity RMSE [m/s]");
title("OFDM");
theme(gcf, "light");