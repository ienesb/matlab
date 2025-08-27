function [tau_hat, nu_hat] = refine_parameters(tau_hat, nu_hat, H, params)
    v2struct(params);

    N_tau = 2000;
    N_nu = 2000;

    tau_points = linspace(-tau_res/2, tau_res/2, N_tau) + tau_hat;
    nu_points = linspace(-nu_res/2, nu_res/2, N_nu) + nu_hat;

    b_vectors = getb(tau_points, params);
    c_vectors = getc(nu_points, params);

    mf = abs(b_vectors' * H * conj(c_vectors));
    
    % figure;
    % mesh(mf);
    
    [~, idx] = max(mf(:));
    [tau_idx, nu_idx] = ind2sub(size(mf), idx);
    
    tau_hat = tau_points(tau_idx);
    nu_hat = nu_points(nu_idx);
end