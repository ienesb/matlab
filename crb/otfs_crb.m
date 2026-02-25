%% CRB and MLE for pilot-only and full data OTFS (QPSK) (CRB Computed)
clear;
% close all;
% clc;

options = optimoptions('fmincon', ...
                  'Algorithm','interior-point', ...
                  'Display','off', ...
                  'MaxFunctionEvaluations',1e5, ...
                  'FiniteDifferenceType','central', ...
                  'FiniteDifferenceStepSize',[1e-2; 1e-10], ...  % (Hz step, seconds step)
                  'StepTolerance',1e-14, ...
                  'OptimalityTolerance',1e-12);

nMonteCarlo = 100;

N = 70;
M = 50;

M_order = 4;

c = 3e8;
fc = 30*1e9; % 30 GHz
lambda = c / fc;

deltaf = 200*1e3; % 200 kHz
Tcp = 1*1e-6; % 1 us
T = 1/deltaf + Tcp;

p_tx = [-40, 0];
p_rx = [0, 40];

p = [rand*20+80, rand*20-100]; % m

v = rand*60-30; % m/s
delta = rand*10-5; % degree

d_tx = norm(p - p_tx);
d_rx = norm(p - p_rx);
D = norm(p_tx - p_rx);

d_bis = (d_tx + d_rx);
v_bis = v * cosd(delta);

beta = acos((d_tx^2 + d_rx^2 - D^2)/(2*d_tx*d_rx));

tau_gt = d_bis / c;
nu_gt = 2*v_bis*cos(beta/2)/lambda;

delta_tau = 1/(N*deltaf);
delta_nu = 1/(M*T);

SNR_dbs = -40:2:30;
SNR_lins = db2pow(SNR_dbs);

sigma2 = 1;

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


CRBs_pilot_only = zeros(length(SNR_dbs), nMonteCarlo, 4, 4);
CRBs_full_data = zeros(length(SNR_dbs), nMonteCarlo, 4, 4);

tau_errors_pilot_only = zeros(length(SNR_dbs), nMonteCarlo);
nu_errors_pilot_only = zeros(length(SNR_dbs), nMonteCarlo);

tau_errors_full_data = zeros(length(SNR_dbs), nMonteCarlo);
nu_errors_full_data = zeros(length(SNR_dbs), nMonteCarlo);


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

        X_DDd = X_DD
        X_DDd(pilot_indices) = 0;

        X_TF = isfft(X_DD);
        X_TFp = isfft(X_DDp);
        X_TFd = isfft(X_DDd);

        ns = (0:(N-1)).';
        ms = 0:(M-1);
        
        H_TF = alpha_gt * exp(-1j*2*pi*deltaf*ns*tau_gt) * exp(1j*2*pi*T*ms*nu_gt);
    
        Z = (randn(N, M) + 1j*randn(N, M)) * sqrt(sigma2/2);
        % Z = 0;
        Y_TF = H_TF .* X_TF + Z;
        y = Y_TF(:);

        Y_DD = sfft(Y_TF);

        %% Pilot Only Parameter Estimation
        % [tau_hat, nu_hat] = my_otfs_ce(Y_TF, deltaf, T);

        % eta0 = [nu_hat; tau_hat];   % initial guess
        eta0 = [nu_gt-5; tau_gt-0.5e-7];   % initial guess
        lb = [nu_gt-100000; tau_gt-1e-6];
        ub = [nu_gt+100000; tau_gt+1e-6];

        [eta_opt, ~] = fmincon( ...
            @(eta) cost_tau_nu(eta, deltaf, T, X_TFp, y), ...
            eta0, [], [], [], [], lb, ub, [], options);

        nu_hat  = eta_opt(1);
        tau_hat = eta_opt(2);

        tau_errors_pilot_only(SNR_idx, mc_idx) = abs(tau_gt - tau_hat);
        nu_errors_pilot_only(SNR_idx, mc_idx) = abs(nu_gt - nu_hat);

        CRB_pilot_only = getCRB_OTFS(X_DD, deltaf, T, sigma2, alpha_gt, nu_gt, tau_gt, pilot_indices);
        CRBs_pilot_only(SNR_idx, mc_idx, :, :) = CRB_pilot_only;
        
        %% Full Data Parameter Estimation

        [eta_opt, ~] = fmincon( ...
            @(eta) cost_tau_nu(eta, deltaf, T, X_TF, y), ...
            eta0, [], [], [], [], lb, ub, [], options);

        nu_hat  = eta_opt(1);
        tau_hat = eta_opt(2);

        tau_errors_full_data(SNR_idx, mc_idx) = abs(tau_gt - tau_hat);
        nu_errors_full_data(SNR_idx, mc_idx) = abs(nu_gt - nu_hat);

        CRB_full_data = getCRB_OTFS(X_DD, deltaf, T, sigma2, alpha_gt, nu_gt, tau_gt, ones(N, M));
        CRBs_full_data(SNR_idx, mc_idx, :, :) = CRB_full_data;
        
    end
end
toc

tau_errors_pilot_only = sqrt(mean(tau_errors_pilot_only.^2, 2));
range_errors_pilot_only = tau_errors_pilot_only .* c;

tau_errors_full_data = sqrt(mean(tau_errors_full_data.^2, 2));
range_errors_full_data = tau_errors_full_data .* c;

nu_errors_pilot_only = sqrt(mean(nu_errors_pilot_only.^2, 2));
velocity_errors_pilot_only = nu_errors_pilot_only .* (lambda / (2*cos(beta/2)));

nu_errors_full_data = sqrt(mean(nu_errors_full_data.^2, 2));
velocity_errors_full_data = nu_errors_full_data .* (lambda / (2*cos(beta/2)));

CRB_delay_pilot_only = CRBs_pilot_only(:, :, 4, 4);
CRB_delay_pilot_only = mean(CRB_delay_pilot_only, 2);
CRB_delay_pilot_only = sqrt(CRB_delay_pilot_only);
CRB_range_pilot_only = CRB_delay_pilot_only .* c;

CRB_delay_full_data = CRBs_full_data(:, :, 4, 4);
CRB_delay_full_data = mean(CRB_delay_full_data, 2);
CRB_delay_full_data = sqrt(CRB_delay_full_data);
CRB_range_full_data = CRB_delay_full_data .* c;

CRB_doppler_pilot_only = CRBs_pilot_only(:, :, 3, 3);
CRB_doppler_pilot_only = mean(CRB_doppler_pilot_only, 2);
CRB_doppler_pilot_only = sqrt(CRB_doppler_pilot_only);
CRB_velocity_pilot_only = CRB_doppler_pilot_only .* (lambda / (2*cos(beta/2)));

CRB_doppler_full_data = CRBs_full_data(:, :, 3, 3);
CRB_doppler_full_data = mean(CRB_doppler_full_data, 2);
CRB_doppler_full_data = sqrt(CRB_doppler_full_data);
CRB_velocity_full_data = CRB_doppler_full_data .* (lambda / (2*cos(beta/2)));

colors = lines(4);

figure;
semilogy(SNR_dbs, tau_errors_pilot_only, "Color", colors(2, :), "LineStyle", '-', 'LineWidth', 2); hold on;
semilogy(SNR_dbs, CRB_delay_pilot_only, "Color", colors(2, :), "LineStyle", '--', 'LineWidth', 2);
semilogy(SNR_dbs, tau_errors_full_data, "Color", colors(3, :), "LineStyle", '-', 'LineWidth', 2);
semilogy(SNR_dbs, CRB_delay_full_data, "Color", colors(3, :), "LineStyle", '--', 'LineWidth', 2);
grid on;
legend("Pilot Only MLE", "Pilot Only CRB", "Full Data MLE", "Full Data CRB");
xlabel("SNR [dB]");
ylabel("Delay RMSE [s]");
title("OTFS");
theme(gcf, "light");

figure;
semilogy(SNR_dbs, nu_errors_pilot_only, "Color", colors(2, :), "LineStyle", '-', 'LineWidth', 2); hold on;
semilogy(SNR_dbs, CRB_doppler_pilot_only, "Color", colors(2, :), "LineStyle", '--', 'LineWidth', 2);
semilogy(SNR_dbs, nu_errors_full_data, "Color", colors(3, :), "LineStyle", '-', 'LineWidth', 2);
semilogy(SNR_dbs, CRB_doppler_full_data, "Color", colors(3, :), "LineStyle", '--', 'LineWidth', 2);
grid on;
legend("Pilot Only MLE", "Pilot Only CRB", "Full Data MLE", "Full Data CRB");
xlabel("SNR [dB]");
ylabel("Doppler RMSE [Hz]");
title("OTFS");
theme(gcf, "light");

figure;
semilogy(SNR_dbs, range_errors_pilot_only, "Color", colors(2, :), "LineStyle", '-', 'LineWidth', 2); hold on;
semilogy(SNR_dbs, CRB_range_pilot_only, "Color", colors(2, :), "LineStyle", '--', 'LineWidth', 2);
semilogy(SNR_dbs, range_errors_full_data, "Color", colors(3, :), "LineStyle", '-', 'LineWidth', 2);
semilogy(SNR_dbs, CRB_range_full_data, "Color", colors(3, :), "LineStyle", '--', 'LineWidth', 2);
grid on;
legend("Pilot Only MLE", "Pilot Only CRB", "Full Data MLE", "Full Data CRB");
xlabel("SNR [dB]");
ylabel("Range RMSE [m]");
title("OTFS");
theme(gcf, "light");

figure;
semilogy(SNR_dbs, velocity_errors_pilot_only, "Color", colors(2, :), "LineStyle", '-', 'LineWidth', 2); hold on;
semilogy(SNR_dbs, CRB_velocity_pilot_only, "Color", colors(2, :), "LineStyle", '--', 'LineWidth', 2);
semilogy(SNR_dbs, velocity_errors_full_data, "Color", colors(3, :), "LineStyle", '-', 'LineWidth', 2);
semilogy(SNR_dbs, CRB_velocity_full_data, "Color", colors(3, :), "LineStyle", '--', 'LineWidth', 2);
grid on;
legend("Pilot Only MLE", "Pilot Only CRB", "Full Data MLE", "Full Data CRB");
xlabel("SNR [dB]");
ylabel("Velocity RMSE [m/s]");
title("OTFS");
theme(gcf, "light");
