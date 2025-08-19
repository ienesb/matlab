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
% [H, Y] = generate_channel_and_received_signal(X, params);

H = generate_channel(params);
Y = generate_received_signal(X, H, params);

% Fn = dftmtx(N)./sqrt(N);
% Fm = dftmtx(M)./sqrt(M);
% 
% DD_map_complex = Fn' * H * Fm;
% DD_map = abs(DD_map_complex).^2;

%% Stage 1 - Initial Channel Estimation
% Step 1: Channel Estimation at Pilot Locations
H_hat = initial_channel_estimation(X, Y, pilot_mask);

% Step 2: Target Detection/Estimation from Channel Estimate
HDD = generate_dd_map(H_hat, params);
% implement CFAR here

% Step 3: Channel Reconstruction from Target Estimates
d0 = norm(params.pT - params.pR);
tau0 = d0 / c;
alpha0 = lambda / (4*pi*d0);
b0 = exp(-1j * 2*pi * (0:N-1)' * params.delta_f * tau0);
c0 = ones(M, 1);
A = kron(conj(c0), b0);

for k = 1:params.K
    tau_k = params.delays(k);
    nu_k = params.dopplers(k);
    b_tau_k = getb(tau_k, params);
    c_nu_k = getc(nu_k, params);
    a_k = kron((c_nu_k), b_tau_k);
    A = [A a_k];
end

alpha_hat = abs(pinv(A)*H(:));

H_hat2 = zeros(size(H_hat));

for k = 1:params.K
    pk = params.targets(k, :);
    vk = params.velocities(k, :);

    d1 = norm(params.pT - pk);
    d2 = norm(pk - params.pR);
    tau = (d1 + d2) / c;

    v_rel = dot(vk, (params.pR - pk)) / norm(params.pR - pk);
    nu = 2 * v_rel / lambda;

    alpha = alpha_hat(k+1);

    b_tau = getb(tau, params);
    c_nu = getc(nu, params);

    H_hat2 = H_hat2 + alpha * (b_tau * c_nu.');
end

% Add LOS path (k=0)
d0 = norm(params.pT - params.pR);
tau0 = d0 / c;
nu0 = 0;
alpha0 = alpha_hat(1);
b0 = getb(tau0, params);
c0 = getc(nu0, params);

H_hat2 = H_hat2 + alpha0 * (b0 * c0.');

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
