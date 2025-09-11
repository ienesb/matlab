function is_ref_detected = genie_sim(params)
    v2struct(params);
    [X, pilot_mask] = generate_ofdm_symbols(params);
    
    H = generate_channel(params);
    Y = generate_received_signal(X, H, params);
    
    % H_hat = initial_channel_estimation(X, Y, pilot_mask);
    H_hat = zeros(N, M);
    H_hat = estimate_channel(X, Y, H_hat, pilot_mask, "RF"); % Initial Channel Estimation
    
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

end


