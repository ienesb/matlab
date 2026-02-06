%% mcrb.m -> MCRB and MLE for pilot-only, data aided, and genie OFDM (M-QAM ranging) (MLE parts not complete!!!)
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

nMonteCarlo = 100;

N = 70;
M = 50;

M_orders = [4, 16, 64];

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

SNR_dbs = -30:60;
SNR_lins = db2pow(SNR_dbs);

MCRBs = zeros(4, 4, length(SNR_dbs), length(M_orders), nMonteCarlo);
tau_errors = zeros(length(SNR_dbs), length(M_orders), nMonteCarlo);

% pilot_periods = [10, 2, 2, 1; 5, 5, 1, 1];
pilot_period = [2; 5];

np = pilot_period(1);
mp = pilot_period(2);

tic
for M_order_idx = 1:length(M_orders)
    M_order = M_orders(M_order_idx);
    for SNR_idx = 1:length(SNR_dbs)
        SNR_db = SNR_dbs(SNR_idx);
        SNR_lin = SNR_lins(SNR_idx);
        % rng(0);
        parfor mc_idx = 1:nMonteCarlo
            sigma2 = 1;
            alpha_gt = sqrt(SNR_lin * sigma2) * exp(1j*rand*2*pi);
            
            symbols = randi(M_order, N, M)-1;
            X = qammod(symbols, M_order, UnitAveragePower=true);
            
            ns = (0:(N-1)).';
            ms = 0:(M-1);
            
            H = alpha_gt * exp(-1j*2*pi*deltaf*ns*tau_gt) * exp(1j*2*pi*T*ms*nu_gt);
            
            mu_gt = H .* X;
            mu_gt = mu_gt(:);
        
            Z = (randn(N, M) + 1j*randn(N, M)) * sqrt(sigma2/2);
            Y = H .* X + Z;
            y = Y(:);
    
            %% Channel Estimation            
            pilot_indices = zeros(N, M);
            pilot_indices(1:np:end, 1:mp:end) = 1;
            pilot_indices = logical(pilot_indices);
    
            X_hat = zeros(N, M);
            
            X_hat(pilot_indices) = X(pilot_indices);
    
            eta0 = [nu_gt-5; tau_gt-0.5e-7];   % initial guess
    
            lb = [nu_gt-10; tau_gt-1e-7];
            ub = [nu_gt+10; tau_gt+1e-7];
    
            [eta_opt, ~] = fmincon(@(eta) cost_tau_nu_pilots(eta, deltaf, T, X_hat, y, pilot_indices), ...
                      eta0, [], [], [], [], lb, ub, [], options);
    
            nu_hat  = eta_opt(1);
            tau_hat = eta_opt(2);
            
            E = exp(-1j*2*pi*deltaf*ns*tau_hat) * exp(1j*2*pi*T*ms*nu_hat);
            
            C = (E .* X_hat);        % use true pilots (X on pilots)
                        
            Cp = C(pilot_indices);
            Yp = Y(pilot_indices);
            
            alpha_hat = (Cp' * Yp) / (Cp' * Cp);
            
            H_hat = alpha_hat * E;     % full grid channel (phase terms apply everywhere)
    
            %% Data-Aided Parameter Estimation
            % X_hat = (Y .* conj(H_hat)) ./ (abs(H_hat).^2 + sigma2);
            X_soft = (Y .* conj(H_hat)) ./ (abs(H_hat).^2 + sigma2);
            sym_hat = qamdemod(X_soft, M_order, UnitAveragePower=true);
            X_hat   = qammod(sym_hat, M_order, UnitAveragePower=true);
            
            X_hat(pilot_indices) = X(pilot_indices);
    
            [eta_opt, Jmin] = fmincon( ...
                @(eta) cost_tau_nu(eta, deltaf, T, X_hat, y), ...
                eta0, [], [], [], [], lb, ub, [], options);
            
            nu_opt = eta_opt(1);
            tau_opt = eta_opt(2);
    
            tau_errors(SNR_idx, M_order_idx, mc_idx) = abs(tau_gt - tau_opt);
    
            %% MCRB Calculation
            
            MCRB = getMCRB(X_hat, mu_gt, deltaf, T, sigma2, alpha_gt, nu_gt, tau_gt);
            MCRBs(:, :, SNR_idx, M_order_idx, mc_idx) = MCRB;
        end
    end
end
toc

MCRBs = mean(MCRBs, 5);
MCRB_ranges = sqrt(MCRBs(4, 4, :, :) * c^2);
MCRB_ranges = squeeze(MCRB_ranges);

tau_errors = sqrt(mean(tau_errors.^2, 3));
range_errors = tau_errors .* c;

%% CRB Pilot Only
K = ceil((N - 1) / np);
L = ceil((M - 1) / mp);

absP = (K + 1) * (L + 1);

CRB_range_pilot_only = 12 / (K * (K + 2) * absP * np^2) * c^2 ./ (8 * pi^2 * SNR_lins * deltaf^2);
CRB_range_pilot_only = sqrt(CRB_range_pilot_only);


%% CRB Full Data
np = 1;
mp = 1;
K = ceil((N - 1) / np);
L = ceil((M - 1) / mp);

absP = (K + 1) * (L + 1);

CRB_range_genie = 12 / (K * (K + 2) * absP * np^2) * c^2 ./ (8 * pi^2 * SNR_lins * deltaf^2);
CRB_range_genie = sqrt(CRB_range_genie);


%% Figure
figure;
h_mcrb = semilogy(SNR_dbs, MCRB_ranges, 'LineWidth', 2); hold on;
h_crlb_pilot_only = semilogy(SNR_dbs, CRB_range_pilot_only, '--', 'LineWidth', 2);
h_crlb_genie = semilogy(SNR_dbs, CRB_range_genie, '--', 'LineWidth', 2);
legend("MCRB 4-QAM", "MCRB 16-QAM", "MCRB 64-QAM", "CRB Pilot only (QPSK)", "CRB Full Data (QPSK)");
grid on;
xlabel("SNR [dB]");
ylabel("Range RMSE [m]");
