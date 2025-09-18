function is_ref_detected = pilot_only_sim(params, is_genie)
    if nargin < 2
        is_genie = 0;
    end
    v2struct(params);
    [X, pilot_mask] = generate_ofdm_symbols(params);
    
    H = generate_channel(params);
    Y = generate_received_signal(X, H, params);
    
    H_hat = zeros(N, M);
    if is_genie
        mask = logical(ones(N, M));
    else
        mask = pilot_mask;
    end
    H_hat = estimate_channel(X, Y, H_hat, mask, "RF"); % Initial Channel Estimation
    
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


