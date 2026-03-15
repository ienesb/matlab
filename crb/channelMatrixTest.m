clear;
close all;
clc;

type = "OFDM";

initialize_parameters;

generate_indices;

SNR_idx = 21;
SNR_db = SNR_dbs(SNR_idx);
SNR_lin = SNR_lins(SNR_idx);

% rng(0);

alpha_gt = sqrt(SNR_lin * sigma2) * exp(1j*rand*2*pi);

symbols = randi(M_order, N, M)-1;
X = qammod(symbols, M_order, UnitAveragePower=true);
x = X(:);

tic
Y1 = channel_noiseless(X, alpha_gt, tau_gt, nu_gt, deltaf, Ts);
toc
Y1 = Y1(:);

tic
H = getChannelMatrix(N, M, tau_gt, nu_gt, alpha_gt, deltaf, Ts);    
Y2 = H * x;
toc

tic
pinvH1 = pinv(H);
toc

tic
pinvH2 = getChannelMatrixPinv(N, M, tau_gt, nu_gt, alpha_gt, deltaf, Ts);
toc



diff = sum(abs(Y2 - Y1))
diffRatio = diff / sum(abs(Y1))

figure; hold on;
plot(real(Y1));
plot(real(Y2));
legend("Y1", "Y2");