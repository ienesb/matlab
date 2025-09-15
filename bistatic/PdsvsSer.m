clear;
close all;
clc;

rcs2_dB = 1.5;
nIter = 1;
pilot_ratio = 0;
monteCarlo = 500;
is_genie = 0;
is_data_only = 0;

params = init_simulation_params(rcs2_dB, pilot_ratio, nIter, monteCarlo, is_genie, is_data_only);
v2struct(params);

sers = 0.4:0.01:0.6;

Pds = zeros(size(sers));
tic
for idx = 1:length(sers)
    ser = sers(idx);
    Pd = 0;
    parfor m = 1:monteCarlo % parfor
        
        [X, ~] = generate_ofdm_symbols(params);
        % pilot_mask = zeros(size(X));
        
        H = generate_channel(params);
        Y = generate_received_signal(X, H, params);
        
        X_hat = X;
                
        nErrors = round(N*M*ser);
        linIdx = randperm(N*M, nErrors).';
        
        error = randi([1 3], nErrors, 1);
        
        data = qamdemod(X_hat(linIdx), 4);
        
        data = data + error;
        data = mod(data, 4);
        
        X_hat(linIdx) = qammod(data, 4);
        % X_hat(linIdx) = 0;
        
        %% Stage 3 - Iterative Refinement
        % H_hat = refine_channel_estimate(X_hat, Y, pilot_mask, alphas, sigma2); % bunun yerine estimate_channel kullan!!!
        H_hat = zeros(size(X));
        mask = logical(ones(size(X_hat))); %#ok<LOGL>
        SNRh = sum(abs(alphas).^2) / sigma2;
        H_hat = estimate_channel(X_hat, Y, H_hat, mask, "LMMSE", SNRh);
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
end
toc
figure;
plot(sers, Pds);