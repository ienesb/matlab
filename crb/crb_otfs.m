%% CRB for OTFS
clear;
% close all;
clc;

nMonteCarlo = 10;

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

tau_gt = tau_gt * 10;
nu_gt = nu_gt * 10;

SNR_dbs = -30:2:30;
% SNR_dbs = 0;
SNR_lins = db2pow(SNR_dbs);

% pilot_periods = [10, 2, 2, 1; 5, 5, 1, 1];
pilot_period = [2; 5];

np = pilot_period(1);
mp = pilot_period(2);

pilot_indices = zeros(N, M);
pilot_indices(1:np:end, 1:mp:end) = 1;
pilot_indices = logical(pilot_indices);

data_indices = ~pilot_indices;

CRBs_mono = zeros(length(SNR_dbs), nMonteCarlo, 4, 4);
CRBs_mono_ICI = zeros(length(SNR_dbs), nMonteCarlo, 4, 4);
CRBs_pilot_only = zeros(length(SNR_dbs), nMonteCarlo, 4, 4);

tic
for SNR_idx = 1:length(SNR_dbs)
    SNR_db = SNR_dbs(SNR_idx);
    SNR_lin = SNR_lins(SNR_idx);
    % rng(0);
    parfor mc_idx = 1:nMonteCarlo
        sigma2 = 1;
        alpha_gt = sqrt(SNR_lin * sigma2) * exp(1j*rand*2*pi);
        
        symbols = randi(M_order, N, M)-1;
        X_DD = qammod(symbols, M_order, UnitAveragePower=true);
        
        S = ifft(X_DD, M, 2) * sqrt(M);
        
        S = S(:);

        CRB_pilot_only = getCRB_OTFS(X_DD, deltaf, T, sigma2, alpha_gt, nu_gt, tau_gt, pilot_indices);
        CRBs_pilot_only(SNR_idx, mc_idx, :, :) = CRB_pilot_only;

        CRB_mono = getCRB_OTFS(X_DD, deltaf, T, sigma2, alpha_gt, nu_gt, tau_gt, true(N, M));
        CRBs_mono(SNR_idx, mc_idx, :, :) = CRB_mono;

        CRB_mono_ICI = getCRB_OTFS_ICI(S, N, M, deltaf, T, sigma2, alpha_gt, nu_gt, tau_gt);
        CRBs_mono_ICI(SNR_idx, mc_idx, :, :) = CRB_mono_ICI;
    end
end
toc

CRB_delay_mono = CRBs_mono(:, :, 4, 4);
CRB_delay_mono = mean(CRB_delay_mono, 2);
CRB_delay_mono = sqrt(CRB_delay_mono);

CRB_delay_mono_ICI = CRBs_mono_ICI(:, :, 4, 4);
CRB_delay_mono_ICI = mean(CRB_delay_mono_ICI, 2);
CRB_delay_mono_ICI = sqrt(CRB_delay_mono_ICI);

CRB_doppler_mono = CRBs_mono(:, :, 3, 3);
CRB_doppler_mono = mean(CRB_doppler_mono, 2);
CRB_doppler_mono = sqrt(CRB_doppler_mono);

CRB_doppler_mono_ICI = CRBs_mono_ICI(:, :, 3, 3);
CRB_doppler_mono_ICI = mean(CRB_doppler_mono_ICI, 2);
CRB_doppler_mono_ICI = sqrt(CRB_doppler_mono_ICI);


colors = lines(4);

figure;
semilogy(SNR_dbs, CRB_delay_mono, "Color", colors(3, :), "LineStyle", '--', 'LineWidth', 2); hold on;
semilogy(SNR_dbs, CRB_delay_mono_ICI, "Color", colors(4, :), "LineStyle", '--', 'LineWidth', 2);
grid on;
legend("Genie", "Genie (ISI/ICI)");
xlabel("SNR [dB]");
ylabel("Delay RMSE [s]");
title("OTFS");

figure;
semilogy(SNR_dbs, CRB_doppler_mono, "Color", colors(3, :), "LineStyle", '--', 'LineWidth', 2); hold on;
semilogy(SNR_dbs, CRB_doppler_mono_ICI, "Color", colors(4, :), "LineStyle", '--', 'LineWidth', 2);
grid on;
legend("Genie", "Genie (ISI/ICI)");
xlabel("SNR [dB]");
ylabel("Doppler RMSE [Hz]");
title("OTFS");
