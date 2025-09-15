clear;
close all;
clc;

load("results/sim.mat");


%% Stage 1 - Initial Channel Estimation
% Step 1: Channel Estimation at Pilot Locations
H_hat = initial_channel_estimation(X, Y, pilot_mask);

% Step 2: Target Detection/Estimation from Channel Estimate 
HDD = generate_dd_map(H_hat, params);

tau = 0.2332*1e-6;
nu = -1.6457;

tic
max_inner = 0;
phis = zeros(N, M, 9417);
for idx1 = 62:280 % 219
    tau = delay_array(idx1);
    b = getb(tau, params);
    for idx2 = 108:150 % 43
        nu = doppler_array(idx2);
        c = getc(nu, params);
        phi = b * c.';
        temp = abs(trace(H_hat' * phi));
        phis(:, :, (idx1-62)*43+idx2-108+1) = phi;
        if temp > max_inner
            max_inner = temp;
            % idx1
            % idx2
        end
    end
end
toc

I = repmat(eye(M), 1, 1, 9417);
tic
temp = I .* pagemtimes(H_hat', phis);
temp = sum(temp, [1, 2]);
temp = abs(temp(:));
[~, max_idx] = max(temp);
toc

return
tic
[tau_hats, nu_hats] = detect_targets(HDD, params);
toc
K_hat = length(tau_hats);

for k = 1:K_hat
    [tau_hats(k), nu_hats(k)] = refine_parameters(tau_hats(k), nu_hats(k), H_hat, params);
    [targetIdx, false_alarm] = getDetectedTarget(tau_hats(k), nu_hats(k), params);
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

