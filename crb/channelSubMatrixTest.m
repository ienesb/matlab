clear;
close all;
clc;

type = "OFDM";

initialize_parametersv2;

generate_indices;

SNR_idx = 1;

SNR_db = SNR_dbs(SNR_idx);
path_gains_db = SNR_db + relative_powers;
path_gains = db2pow(path_gains_db); % add random phase shift !!!!!!

% rng(0);

symbols = randi(M_order, N, M)-1;
X = qammod(symbols, M_order, UnitAveragePower=true);
x = X(:);

tic
Y1 = channel(X, path_gains, time_delays, doppler_shifts, deltaf, Ts, 0);
toc

tic
Y2 = channelv2(X, path_gains, time_delays, doppler_shifts, deltaf, Ts, 0);
toc

tic
X_hat = channel_inverse(Y1, path_gains, time_delays, doppler_shifts, deltaf, Ts);
toc

Y1 = Y1(:);
Y2 = Y2(:);


diff = sum(abs(Y2 - Y1))
diffRatio = diff / sum(abs(Y1))

figure; hold on;
plot(real(Y1));
plot(real(Y2));
legend("Y1", "Y2");

figure; hold on;
plot(imag(Y1));
plot(imag(Y2));
legend("Y1", "Y2");

figure; hold on;
plot(real(X(:)));
plot(real(X_hat(:)));
legend("X", "X hat");

figure; hold on;
plot(imag(X(:)));
plot(imag(X_hat(:)));
legend("X", "X hat");