function [target_Idx, false_alarm] = getDetectedTarget(tau_hat, nu_hat, params)
    v2struct(params);

    tau_idx_hat = tau2idx(tau_hat, Tsym, N_fft);
    nu_idx_hat = nu2idx(nu_hat, delta_f, M_fft);

    diff_tau_idx = tau_idx - tau_idx_hat;
    diff_nu_idx = nu_idx - nu_idx_hat;

    [~, target_Idx] = min(diff_tau_idx.^2 + diff_nu_idx.^2);

    if abs(taus(target_Idx) - tau_hat) < tau_res && abs(nus(target_Idx) - nu_hat) < nu_res
        false_alarm = 0;
    else
        false_alarm = 1;
    end
end