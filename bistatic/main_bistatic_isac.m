%% Transmit Signal
[X, Xp, Xd, pilot_mask, data_mask] = generate_ofdm_symbols(params);

%% Channel Model
H = generate_channel(params);
Y = generate_received_signal(X, H, params);

%% Stage 1 - Initial Channel Estimation
% Step 1: Channel Estimation at Pilot Locations
H_hat = initial_channel_estimation(X, Y, pilot_mask);

% Step 2: Target Detection/Estimation from Channel Estimate 
if ~is_data_only
    HDD = generate_dd_map(H_hat, params);
    
    % plotDDMap(HDD, params, 0, 1)
    % plotRangeProfile(HDD, params);
    
    [tau_hats, nu_hats] = detect_targets(HDD, params);
    
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
    
    is_ref_detected_array = zeros(1, nIter);
    for iter = 1:nIter
        %% Stage 2 - Data Demodulation
        X_hat = data_demodulation(Y, H_hat, Xp, data_mask, params);
        % ser = getSer(X, X_hat, data_mask);
        
        %% Stage 3 - Iterative Refinement
        H_hat = refine_channel_estimate(X_hat, Y, pilot_mask, alphas, sigma2);

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
        is_ref_detected_array(iter) = is_ref_detected;
    end
end



