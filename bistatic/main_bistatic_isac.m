% Bistatic OFDM ISAC - Modular Implementation
% Based on: "Bridging the Gap via Data-Aided Sensing"

% === Main Script ===
clear;
close all;

addpath('modules');  % Optional folder for functions

c = 3e8;

%% Simulation Parameters
params = init_simulation_params();
lambda = params.lambda;
N = params.N;
M = params.M;

%% Transmit Signal
[X, Xp, Xd, pilot_mask, data_mask] = generate_ofdm_symbols(params);

%% Channel Model
[H, Y] = generate_channel_and_received_signal(X, params);

%% Stage 1 - Initial Channel Estimation (Pilot-only)
H_hat = initial_channel_estimation(X, Y, pilot_mask);
HDD = generate_dd_map(H_hat, params);

d0 = norm(params.pT - params.pR);
tau0 = d0 / c;
alpha0 = lambda / (4*pi*d0);
b0 = exp(-1j * 2*pi * (0:N-1)' * params.delta_f * tau0);
c0 = ones(M, 1);
A = kron(conj(c0), b0);


for k = 1:params.K
    tau_k = params.delays(k);
    nu_k = params.dopplers(k);
    b_tau_k = exp(-1j * 2*pi * (0:N-1)' * params.delta_f * tau_k);
    c_nu_k = exp(1j * 2*pi * (0:M-1)' * params.Tsym * nu_k);
    a_k = kron(conj(c_nu_k), b_tau_k);
    A = [A a_k];
end

alpha_hat = A'*H_hat(:);

%% Stage 2 - Data Demodulation
X_hat = data_demodulation(Y, H_hat, Xp, data_mask, params);

%% Stage 3 - Iterative Refinement (1 iteration shown)
H_hat = refine_channel_estimate(X_hat, Y, pilot_mask);
X_hat = data_demodulation(Y, H_hat, Xp, data_mask, params);

%% Stage 4 - Final Detection (Delay-Doppler Map)
[DD_map, peaks] = delay_doppler_detection(H, params);

%% Display
visualize_results(DD_map);

%% Evaluation Metric
Pd = evaluate_detection(peaks, params.targets);
disp(["Probability of Detection: ", num2str(Pd)]);
