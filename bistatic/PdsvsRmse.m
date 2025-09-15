clear;
close all;
clc;

rcs2_dB = 1.5;
nIter = 1;
pilot_ratio = 0;
monteCarlo = 50;
is_genie = 0;
is_data_only = 0;

params = init_simulation_params(rcs2_dB, pilot_ratio, nIter, monteCarlo, is_genie, is_data_only);
v2struct(params);

sigma2s = 0.75:0.05:1.25 + 0.5;

Pds = zeros(size(sigma2s));
rmse_dBs = zeros(size(sigma2s));

tic
for idx = 1:length(sigma2s)
    sigma2x = sigma2s(idx);
    Pd = 0;
    rmse = 0;
    parfor m = 1:monteCarlo % parfor
        [X, ~] = generate_ofdm_symbols(params);
        % pilot_mask = zeros(size(X));
        pilot_mask = logical(zeros(size(X))); %#ok<*LOGL>

        Z = sqrt(sigma2x/2) * (randn(N, M) + 1j*randn(N, M));
        X_hat = X + Z;
        [~, r] = getSer(X, X_hat, pilot_mask);
        rmse = rmse + r;
        
        H = generate_channel(params);
        Y = generate_received_signal(X, H, params);
        
        
        %% Stage 3 - Iterative Refinement
        % H_hat = refine_channel_estimate(X_hat, Y, pilot_mask, alphas, sigma2); % bunun yerine estimate_channel kullan!!!
        H_hat = zeros(size(X));
        SNRh = sum(abs(alphas).^2) / sigma2;
        H_hat = estimate_channel(X_hat, Y, H_hat, ~pilot_mask, "LMMSE", SNRh);
        % H_hat = pilotBasedChannelNormalization(H_hat, pilot_mask);
        
        %% Stage 4 - Final Detection (Delay-Doppler Map)
        HDD = generate_dd_map(H_hat, params);
        
        [tau_hats, nu_hats] = detect_targets(HDD, params);
        
        K_hat = length(tau_hats);
        is_ref_detected = 0;
        for k = 1:K_hat
            [tau_hats(k), nu_hats(k)] = refine_parameters(tau_hats(k), nu_hats(k), H_hat, params);
            [targetIdx, false_alarm] = getDetectedTarget(tau_hats(k), nu_hats(k), params);
            if targetIdx == ref_target_idx && ~false_alarm
                is_ref_detected = 1;
            end
        end
        Pd = Pd + is_ref_detected;
    end
    
    Pds(idx) = Pd / monteCarlo;
    rmse_dBs(idx) = pow2db(rmse / monteCarlo);
    
end
toc

figure;
plot(rmse_dBs, Pds);