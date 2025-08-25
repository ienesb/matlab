clear;
clc;
close all;

addpath('modules');  % Optional folder for functions

%% Simulation Parameters
params = init_simulation_params();
v2struct(params);

%% Transmit Signal
[X, Xp, Xd, pilot_mask, data_mask] = generate_ofdm_symbols(params);

%% Channel Model
H = generate_channel(params);
Y = generate_received_signal(X, H, params);

%% Stage 1 - Initial Channel Estimation
% Step 1: Channel Estimation at Pilot Locations
H_hat = initial_channel_estimation(X, Y, pilot_mask);

% Step 2: Target Detection/Estimation from Channel Estimate Assuming single
% target detected !!!!!!!!!!!!!!!
HDD = generate_dd_map(H_hat, params);

[tau_hat, nu_hat] = detect_targets(HDD, params);

% tau_idx_hat = tau2idx(tau_hat, Tsym, N_fft);
% nu_idx_hat = nu2idx(nu_hat, delta_f, M_fft);

[targetIdx, false_alarm] = getDetectedTarget(tau_hat, nu_hat, params)

return

% Step 3: Channel Reconstruction from Target Estimates

b_tau = getb(tau_hat, params);
c_nu = getc(nu_hat, params);
A = kron((c_nu), b_tau);

alpha_hat = (abs(pinv(A)*H_hat(:))) / pilot_ratio;
% alpha_hat = 5.3963e-7;

H_hat = zeros(size(H_hat));

H_hat = H_hat + alpha_hat * (b_tau * c_nu.');


%% Stage 2 - Data Demodulation
X_hat = data_demodulation(Y, H_hat, Xp, data_mask, params);
ser = getSer(X, X_hat, data_mask)


%% Stage 3 - Iterative Refinement (1 iteration shown)
H_hat = refine_channel_estimate(X_hat, Y, pilot_mask);
X_hat = data_demodulation(Y, H_hat, Xp, data_mask, params);
ser = getSer(X, X_hat, data_mask)
return
%% Stage 4 - Final Detection (Delay-Doppler Map)
[DD_map, peaks] = delay_doppler_detection(H, params);

%% Display
visualize_results(DD_map);

%% Evaluation Metric
Pd = evaluate_detection(peaks, params.targets);
disp(["Probability of Detection: ", num2str(Pd)]);
