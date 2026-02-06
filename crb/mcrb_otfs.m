%% MCRB not complete!!!!!!!!!!
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

nMonteCarlo = 1;

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

delta_tau = 1/(N*deltaf);
delta_nu = 1/(M*T);

FN = dftmtx(N) / sqrt(N);
FM = dftmtx(M) / sqrt(M);
F = kron(FM, FN);

SNR_dbs = -30:5:30;
% SNR_dbs = 0;
SNR_lins = db2pow(SNR_dbs);

pilot_indices = zeros(N, M);
pilot_indices(35, 25) = 1;
pilot_indices = logical(pilot_indices);

guard_indices = zeros(N, M);
guard_indices(25:45, 15:35) = 1;
guard_indices(35, 25) = 0;
guard_indices = logical(guard_indices);

data_indices = guard_indices + pilot_indices;
data_indices = ~data_indices;

nBits_perMC  = nnz(data_indices) * log2(M_order); % QPSK => 2 bits/symbol

tau_errors_pilot_only = zeros(length(SNR_dbs), nMonteCarlo);

CRBs_mono = zeros(length(SNR_dbs), nMonteCarlo, 4, 4);
CRBs_pilot_only = zeros(length(SNR_dbs), nMonteCarlo, 4, 4);
MCRBs = zeros(length(SNR_dbs), nMonteCarlo, 4, 4);

tic
for SNR_idx = 1:length(SNR_dbs)
    SNR_db = SNR_dbs(SNR_idx);
    SNR_lin = SNR_lins(SNR_idx);
    % rng(0);
    bitErrs_OFDM = zeros(nMonteCarlo,1);
    for mc_idx = 1:nMonteCarlo
        sigma2 = 1;
        alpha_gt = sqrt(SNR_lin * sigma2) * exp(1j*rand*2*pi);
        
        symbols = randi(M_order, N, M)-1;
        X_DD = qammod(symbols, M_order, UnitAveragePower=true);
        X_DD(guard_indices) = 0;
        X_DD(pilot_indices) = 1;

        x_DD = X_DD(:);
        
        X_TF = isfft(X_DD);

        bits_tx = qamdemod(X_DD(data_indices), M_order, UnitAveragePower=true, OutputType="bit");
        
        ns = (0:(N-1)).';
        ms = 0:(M-1);
        
        H_TF = alpha_gt * exp(-1j*2*pi*deltaf*ns*tau_gt) * exp(1j*2*pi*T*ms*nu_gt);
        H_DD = sfft(H_TF);
    
        Z = (randn(N, M) + 1j*randn(N, M)) * sqrt(sigma2/2);
        % Z = 0;
        Y_TF = H_TF .* X_TF + Z;

        Y_DD1 = sfft(Y_TF);
        y_DD1 = Y_DD1(:);
        % Y_DD2 = circular_conv2d(H_DD, X_DD);

        % diff = sum(abs(Y_DD1(:) - Y_DD2(:)));

        kernel = ones(21, 21);
        kernel = -kernel;
        kernel(11, 11) = 30;

        B = real(ifft2( fft2(abs(Y_DD1)) .* fft2(kernel, N, M) ));

        % figure;
        % mesh(abs(X_DD))
        % figure;
        % mesh(abs(Y_DD1))
        % figure;
        % mesh(B)

        [maxVal, linearIdx] = max(B(:));
        [rowIdx, colIdx] = ind2sub([N, M], linearIdx);

        rowIdx = rowIdx - 10;
        colIdx = colIdx - 10;

        tau_hat = (rowIdx-35)*delta_tau;
        nu_hat = (colIdx-25)*delta_nu;

        tau_errors_pilot_only(SNR_idx, mc_idx) = abs(tau_gt - tau_hat);

        C_TF = exp(-1j*2*pi*deltaf*ns*tau_hat) * exp(1j*2*pi*T*ms*nu_hat);
        C_DD = sfft(C_TF);

        H_DD_hat = alpha_gt * C_DD;
        
        % L = FN * H_DD_hat * FM;
        L = FN * H_DD * FM;
        H = F' * diag(L(:)) * F;
        
        x_DD_hat = (H' * H + eye(N*M)) \ H' * y_DD1;
        X_DD_hat = reshape(x_DD_hat, [N, M]);
        X_DD_hat(guard_indices) = 0;
        X_DD_hat(pilot_indices) = 1;

        MCRB = getMCRB_OTFS(X_DD, X_DD_hat, deltaf, T, sigma2, alpha_gt, nu_gt, tau_gt);
        MCRBs(SNR_idx, mc_idx, :, :) = MCRB;

        CRB_pilot_only = getCRB(X_DD, deltaf, T, sigma2, alpha_gt, nu_gt, tau_gt, pilot_indices);
        CRBs_pilot_only(SNR_idx, mc_idx, :, :) = CRB_pilot_only;

        CRB_mono = getCRB(X_DD, deltaf, T, sigma2, alpha_gt, nu_gt, tau_gt, true(N, M));
        CRBs_mono(SNR_idx, mc_idx, :, :) = CRB_mono;
    end
end
toc

tau_errors_pilot_only = sqrt(mean(tau_errors_pilot_only.^2, 2));
range_errors_pilot_only = tau_errors_pilot_only .* c;

MCRB_range = MCRBs(:, :, 4, 4);
MCRB_range = mean(MCRB_range, 2);
MCRB_range = sqrt(MCRB_range * c^2);

CRB_range_pilot_only = CRBs_pilot_only(:, :, 4, 4);
CRB_range_pilot_only = mean(CRB_range_pilot_only, 2);
CRB_range_pilot_only = sqrt(CRB_range_pilot_only * c^2);

CRB_range_mono = CRBs_mono(:, :, 4, 4);
CRB_range_mono = mean(CRB_range_mono, 2);
CRB_range_mono = sqrt(CRB_range_mono * c^2);


colors = lines(3);

figure;
semilogy(SNR_dbs, MCRB_range, "Color", colors(1, :), "LineStyle", '--', 'LineWidth', 2); hold on;
semilogy(SNR_dbs, range_errors_pilot_only, "Color", colors(2, :), "LineStyle", '-', 'LineWidth', 2);
semilogy(SNR_dbs, CRB_range_pilot_only, "Color", colors(2, :), "LineStyle", '--', 'LineWidth', 2);
semilogy(SNR_dbs, CRB_range_mono, "Color", colors(3, :), "LineStyle", '--', 'LineWidth', 2);
grid on;
legend("MCRB", "Pilot Only Estimation MSE", "Pilot Only CRLB", "Genie CRLB");
xlabel("SNR [dB]");
ylabel("Range RMSE [m]");