%% OFDM vs OTFS in terms of BER assumed full CSI (can be removed)
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
v = v * 100;
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

BER_OTFS = zeros(length(SNR_dbs), 1);
BER_OFDM = zeros(length(SNR_dbs), 1);

pilot_indices = zeros(N, M);
% pilot_indices(1:np:end, 1:mp:end) = 1;
pilot_indices = logical(pilot_indices);

data_indices = ~pilot_indices;
nBits_perMC  = nnz(data_indices) * log2(M_order); % QPSK => 2 bits/symbol

tic
for SNR_idx = 1:length(SNR_dbs)
    SNR_db = SNR_dbs(SNR_idx);
    SNR_lin = SNR_lins(SNR_idx);
    % rng(0);
    bitErrs_OTFS = zeros(nMonteCarlo,1);
    bitErrs_OFDM = zeros(nMonteCarlo,1);
    parfor mc_idx = 1:nMonteCarlo
        sigma2 = 1;
        alpha_gt = sqrt(SNR_lin * sigma2) * exp(1j*rand*2*pi);
        
        symbols = randi(M_order, N, M)-1;
        X_OFDM = qammod(symbols, M_order, UnitAveragePower=true);
        X_DD_OTFS = qammod(symbols, M_order, UnitAveragePower=true);
        X_TF_OTFS = isfft(X_DD_OTFS);

        bits_tx = qamdemod(X_DD_OTFS(data_indices), M_order, UnitAveragePower=true, OutputType="bit");
        
        ns = (0:(N-1)).';
        ms = 0:(M-1);
        
        H_TF = alpha_gt * exp(-1j*2*pi*deltaf*ns*tau_gt) * exp(1j*2*pi*T*ms*nu_gt);
    
        Z = (randn(N, M) + 1j*randn(N, M)) * sqrt(sigma2/2);
        Y_TF_OTFS = H_TF .* X_TF_OTFS + Z;
        Y_OFDM = H_TF .* X_OFDM + Z;

        X_soft_TF_OTFS = (Y_TF_OTFS .* conj(H_TF)) ./ (abs(H_TF).^2 + sigma2);
        X_soft_DD_OTFS = sfft(X_soft_TF_OTFS);

        X_soft_OFDM = (Y_OFDM .* conj(H_TF)) ./ (abs(H_TF).^2 + sigma2);

        bits_hat_OTFS = qamdemod(X_soft_DD_OTFS(data_indices), M_order, UnitAveragePower=true, OutputType="bit");
        bitErrs_OTFS(mc_idx) = sum(bits_hat_OTFS ~= bits_tx);

        bits_hat_OFDM = qamdemod(X_soft_OFDM(data_indices), M_order, UnitAveragePower=true, OutputType="bit");
        bitErrs_OFDM(mc_idx) = sum(bits_hat_OFDM ~= bits_tx);
    end
    BER_OTFS(SNR_idx) = sum(bitErrs_OTFS) / (nMonteCarlo * nBits_perMC);
    BER_OFDM(SNR_idx) = sum(bitErrs_OFDM) / (nMonteCarlo * nBits_perMC);
end
toc


% BER
figure;
semilogy(SNR_dbs, BER_OTFS, 'LineWidth', 2); hold on;
semilogy(SNR_dbs, BER_OFDM, 'LineWidth', 2); hold on;
grid on;
xlabel("SNR [dB]");
ylabel("BER");
legend("BER (OTFS)", "BER (OFDM)");
ylim([1e-6, 1]);