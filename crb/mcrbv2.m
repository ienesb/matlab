%% MCRB and MLE for pilot-only, data aided, and genie OFDM + BER (QPSK)
clear;
% close all;
clc;

options = optimoptions('fmincon', ...
              'Algorithm','interior-point', ...
              'Display','off', ...
              'MaxFunctionEvaluations',1e5, ...
              'FiniteDifferenceType','central', ...
              'FiniteDifferenceStepSize',[1e-2; 1e-10], ...  % (Hz step, seconds step)
              'StepTolerance',1e-14, ...
              'OptimalityTolerance',1e-12);

nMonteCarlo = 10;

N = 70;
M = 50;

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

SNR_dbs = -30:30;
SNR_lins = db2pow(SNR_dbs);

MCRBs = zeros(4, 4, length(SNR_dbs), nMonteCarlo);
tau_errors = zeros(length(SNR_dbs), nMonteCarlo);
tau_errors_pilot_only = zeros(length(SNR_dbs), nMonteCarlo);
tau_errors_genie = zeros(length(SNR_dbs), nMonteCarlo);

nu_errors = zeros(length(SNR_dbs), nMonteCarlo);
nu_errors_pilot_only = zeros(length(SNR_dbs), nMonteCarlo);
nu_errors_genie = zeros(length(SNR_dbs), nMonteCarlo);

BER_pilotCE   = zeros(length(SNR_dbs), 1);   % BER using pilot-only channel (first detection)
BER_dataAided = zeros(length(SNR_dbs), 1);   % BER after data-aided param est + channel rebuild (second detection)

% pilot_periods = [10, 2, 2, 1; 5, 5, 1, 1];
pilot_period = [4; 7];

np = pilot_period(1);
mp = pilot_period(2);


pilot_indices = zeros(N, M);
pilot_indices(1:np:end, 1:mp:end) = 1;
pilot_indices = logical(pilot_indices);

data_indices = ~pilot_indices;
nBits_perMC  = nnz(data_indices) * 2; % QPSK => 2 bits/symbol

tic
for SNR_idx = 1:length(SNR_dbs)
    SNR_db = SNR_dbs(SNR_idx);
    SNR_lin = SNR_lins(SNR_idx);
    % rng(0);
    bitErrs_pilotCE   = zeros(nMonteCarlo,1);
    bitErrs_dataAided = zeros(nMonteCarlo,1);
    parfor mc_idx = 1:nMonteCarlo
        sigma2 = 1;
        alpha_gt = sqrt(SNR_lin * sigma2) * exp(1j*rand*2*pi);
        
        symbols = randi(4, N, M)-1;
        X = qammod(symbols, 4, UnitAveragePower=true);
        
        ns = (0:(N-1)).';
        ms = 0:(M-1);
        
        H = alpha_gt * exp(-1j*2*pi*deltaf*ns*tau_gt) * exp(1j*2*pi*T*ms*nu_gt);
        
        mu_gt = H .* X;
        mu_gt = mu_gt(:);
    
        Z = (randn(N, M) + 1j*randn(N, M)) * sqrt(sigma2/2);
        Y = H .* X + Z;
        y = Y(:);

        %% Channel Estimation            
        
        % Ground-truth bits on data positions
        bits_tx = qamdemod(X(data_indices), 4, UnitAveragePower=true, OutputType="bit");

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
        
        E = exp(-1j*2*pi*deltaf*ns*tau_hat) * exp(1j*2*pi*T*ms*nu_hat);
        
        C = (E .* X_hat);        % use true pilots (X on pilots)
                    
        Cp = C(pilot_indices);
        Yp = Y(pilot_indices);
        
        alpha_hat = (Cp' * Yp) / (Cp' * Cp);
        
        H_hat = alpha_hat * E;     % full grid channel (phase terms apply everywhere)

        %% Genie
        [eta_opt, ~] = fmincon( ...
            @(eta) cost_tau_nu(eta, deltaf, T, X, y), ...
            eta0, [], [], [], [], lb, ub, [], options);
        
        nu_opt = eta_opt(1);
        tau_opt = eta_opt(2);

        tau_errors_genie(SNR_idx, mc_idx) = abs(tau_gt - tau_opt);
        nu_errors_genie(SNR_idx, mc_idx) = abs(nu_gt - nu_opt);

        %% Data-Aided Parameter Estimation
        % X_hat = (Y .* conj(H_hat)) ./ (abs(H_hat).^2 + sigma2);
        X_soft = (Y .* conj(H_hat)) ./ (abs(H_hat).^2 + sigma2);
        sym_hat = qamdemod(X_soft, 4, UnitAveragePower=true);
        X_hat   = qammod(sym_hat, 4, UnitAveragePower=true);
        
        X_hat(pilot_indices) = X(pilot_indices);


        bits_hat1 = qamdemod(X_soft(data_indices), 4, UnitAveragePower=true, OutputType="bit");
        bitErrs_pilotCE(mc_idx) = sum(bits_hat1 ~= bits_tx);

        [eta_opt, Jmin] = fmincon( ...
            @(eta) cost_tau_nu(eta, deltaf, T, X_hat, y), ...
            eta0, [], [], [], [], lb, ub, [], options);
        
        nu_opt = eta_opt(1);
        tau_opt = eta_opt(2);

        tau_errors(SNR_idx, mc_idx) = abs(tau_gt - tau_opt);
        nu_errors(SNR_idx, mc_idx) = abs(nu_gt - nu_opt);

        % ---------- Rebuild channel using refined (tau_opt, nu_opt) ----------
        E2 = exp(-1j*2*pi*deltaf*ns*tau_opt) * exp(1j*2*pi*T*ms*nu_opt);
        
        % Estimate alpha again using TRUE pilots (important)
        X_pil = zeros(N,M);
        X_pil(pilot_indices) = X(pilot_indices);
        
        C2  = E2 .* X_pil;
        Cp2 = C2(pilot_indices);
        Yp  = Y(pilot_indices);              % already defined earlier, but safe to reuse
        
        alpha_hat2 = (Cp2' * Yp) / (Cp2' * Cp2);
        H_hat2     = alpha_hat2 * E2;
        
        % ---------- BER (data-aided refined channel, second detection) ----------
        X_soft2 = (Y .* conj(H_hat2)) ./ (abs(H_hat2).^2 + sigma2);
        
        bits_hat2 = qamdemod(X_soft2(data_indices), 4, UnitAveragePower=true, OutputType="bit");
        bitErrs_dataAided(mc_idx) = sum(bits_hat2 ~= bits_tx);

        %% MCRB Calculation
        
        MCRB = getMCRB(X_hat, mu_gt, deltaf, T, sigma2, alpha_gt, nu_gt, tau_gt);
        MCRBs(:, :, SNR_idx, mc_idx) = MCRB;
    end
    BER_pilotCE(SNR_idx)   = sum(bitErrs_pilotCE)   / (nMonteCarlo * nBits_perMC);
    BER_dataAided(SNR_idx) = sum(bitErrs_dataAided) / (nMonteCarlo * nBits_perMC);
end
toc

MCRBs = mean(MCRBs, 4);

MCRB_range = sqrt(MCRBs(4, 4, :) * c^2);
MCRB_range = MCRB_range(:);

MCRB_velocity = sqrt(MCRBs(3, 3, :));
MCRB_velocity = MCRB_velocity(:) * (lambda / (2*cos(beta/2)));

tau_errors = sqrt(mean(tau_errors.^2, 2));
range_errors = tau_errors .* c;
nu_errors = sqrt(mean(nu_errors.^2, 2));
velocity_errors = nu_errors .* (lambda / (2*cos(beta/2)));

tau_errors_pilot_only = sqrt(mean(tau_errors_pilot_only.^2, 2));
range_errors_pilot_only = tau_errors_pilot_only .* c;
nu_errors_pilot_only = sqrt(mean(nu_errors_pilot_only.^2, 2));
velocity_errors_pilot_only = nu_errors_pilot_only .* (lambda / (2*cos(beta/2)));

tau_errors_genie = sqrt(mean(tau_errors_genie.^2, 2));
range_errors_genie = tau_errors_genie .* c;
nu_errors_genie = sqrt(mean(nu_errors_genie.^2, 2));
velocity_errors_genie = nu_errors_genie .* (lambda / (2*cos(beta/2)));

%% CRB Pilot Only
K = ceil((N - 1) / np);
L = ceil((M - 1) / mp);

absP = (K + 1) * (L + 1);

CRB_range_pilot_only = 12 / (K * (K + 2) * absP * np^2) * c^2 ./ (8 * pi^2 * SNR_lins * deltaf^2);
CRB_range_pilot_only = sqrt(CRB_range_pilot_only);

CRB_velocity_pilot_only = 12 / (L * (L + 2) * absP * mp^2) * lambda^2 ./ (32 * pi^2 * SNR_lins * T^2 * cos(beta/2)^2);
CRB_velocity_pilot_only = sqrt(CRB_velocity_pilot_only);

%% CRB Full Data
np = 1;
mp = 1;
K = ceil((N - 1) / np);
L = ceil((M - 1) / mp);

absP = (K + 1) * (L + 1);

CRB_range_genie = 12 / (K * (K + 2) * absP * np^2) * c^2 ./ (8 * pi^2 * SNR_lins * deltaf^2);
CRB_range_genie = sqrt(CRB_range_genie);

CRB_velocity_genie = 12 / (L * (L + 2) * absP * mp^2) * lambda^2 ./ (32 * pi^2 * SNR_lins * T^2 * cos(beta/2)^2);
CRB_velocity_genie = sqrt(CRB_velocity_genie);

%% Figure

colors = lines(3);

figure;
semilogy(SNR_dbs, range_errors, "Color", colors(1, :), 'LineWidth', 2); hold on;
semilogy(SNR_dbs, MCRB_range, "Color", colors(1, :), "LineStyle", '--', 'LineWidth', 2);
semilogy(SNR_dbs, range_errors_pilot_only, "Color", colors(2, :),  'LineWidth', 2);
semilogy(SNR_dbs, CRB_range_pilot_only, "Color", colors(2, :), "LineStyle", '--', 'LineWidth', 2);
semilogy(SNR_dbs, range_errors_genie, "Color", colors(3, :), 'LineWidth', 2);
semilogy(SNR_dbs, CRB_range_genie, "Color", colors(3, :), "LineStyle", '--', 'LineWidth', 2);
legend("MLE (Data-aided)", "MCRB (Data-aided)", "MLE (Pilot only)", "CRB (Pilot only)", "MLE (Full Data)", "CRB (Full Data)");
grid on;
xlabel("SNR [dB]");
ylabel("Range RMSE [m]");


figure;
semilogy(SNR_dbs, velocity_errors, "Color", colors(1, :), 'LineWidth', 2); hold on;
semilogy(SNR_dbs, MCRB_velocity, "Color", colors(1, :), "LineStyle", '--', 'LineWidth', 2);
semilogy(SNR_dbs, velocity_errors_pilot_only, "Color", colors(2, :),  'LineWidth', 2);
semilogy(SNR_dbs, CRB_velocity_pilot_only, "Color", colors(2, :), "LineStyle", '--', 'LineWidth', 2);
semilogy(SNR_dbs, velocity_errors_genie, "Color", colors(3, :), 'LineWidth', 2);
semilogy(SNR_dbs, CRB_velocity_genie, "Color", colors(3, :), "LineStyle", '--', 'LineWidth', 2);
legend("MLE (Data-aided)", "MCRB (Data-aided)", "MLE (Pilot only)", "CRB (Pilot only)", "MLE (Full Data)", "CRB (Full Data)");
grid on;
xlabel("SNR [dB]");
ylabel("Velocity RMSE [m/s]");

% BER
figure;
semilogy(SNR_dbs, BER_pilotCE,   'LineWidth', 2); hold on;
semilogy(SNR_dbs, BER_dataAided, 'LineWidth', 2);
grid on;
xlabel("SNR [dB]");
ylabel("BER");
legend("BER (Pilot-only CE, 1st detection)", "BER (Data-aided refined CE, 2nd detection)");
ylim([1e-6, 1]);