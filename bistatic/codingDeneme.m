clear;
close all;
clc;

rcs2_dB = -5;
rcs2 = db2pow(rcs2_dB);

nIter = 1;
pilot_ratio = 0.05;
monteCarlo = 1;
is_genie = 0;

params = init_simulation_params(rcs2_dB, pilot_ratio, nIter, monteCarlo);

v2struct(params);

v2struct(params);
[X, pilot_mask] = generate_ofdm_symbols(params);

H = generate_channel(params);
Y = generate_received_signal(X, H, params);

%% Stage 1 - Initial Channel Estimation
% Step 1: Channel Estimation at Pilot Locations
% H_hat = initial_channel_estimation(X, Y, pilot_mask);
H_hat = zeros(N, M);
% SNRh = sum(abs(alphas).^2) / sigma2;
H_hat = estimate_channel(X, Y, H_hat, pilot_mask, "RF"); % Initial Channel Estimation
H_hat = pilotBasedChannelNormalization(H_hat, pilot_mask);

% Step 2: Target Detection/Estimation from Channel Estimate 
HDD = generate_dd_map(H_hat, params);

% plotDDMap(HDD, params, 0, 1)
% plotRangeProfile(HDD, params);

[tau_hats, nu_hats] = detect_targets(HDD, params);

K_hat = length(tau_hats);

for k = 1:K_hat
    [tau_hats(k), nu_hats(k)] = refine_parameters(tau_hats(k), nu_hats(k), H_hat, params);
    % [targetIdx, false_alarm] = getDetectedTarget(tau_hats(k), nu_hats(k), params);
end

% Step 3: Channel Reconstruction from Target Estimates

A = zeros(N*M, K_hat);
for k = 1:K_hat
    b_tau = getb(tau_hats(k), params);
    c_nu = getc(nu_hats(k), params);
    a_k = kron((c_nu), b_tau);
    A(:, k) = a_k;
end

alpha_hat = (abs(pinv(A)*H_hat(:))) / pilot_ratio;

for k = 1:K_hat
    b_tau = getb(tau_hats(k), params);
    c_nu = getc(nu_hats(k), params);
    temp = alpha_hat(k) * (b_tau * c_nu.');
    temp(pilot_mask) = 0;

    H_hat = H_hat + temp; % pilot_mask entrileri ayni kaliyor.
end

is_ref_detected = zeros(nIter, 1);
%% Stage 2 - Data Demodulation
X_hat = data_demodulation(Y, H_hat, X, pilot_mask, params);
[rmse, ber] = getSer(X, X_hat, pilot_ratio);
% X_hat = hard_demod(X_hat, ~pilot_mask);