function [Y, H_true] = generate_channel_and_received_signal(X, params)
    % Generate the bistatic channel and simulate the received signal Y
    
    N = params.N;
    M = params.M;
    K = params.K;
    fc = params.fc;
    c = 3e8;
    lambda = c / fc;
    
    H_true = zeros(N, M);
    Z = sqrt(params.sigma2/2) * (randn(N, M) + 1j*randn(N, M));
    
    for k = 1:K
        pk = params.targets(k, :);
        vk = params.velocities(k, :);
    
        d1 = norm(params.pT - pk);
        d2 = norm(pk - params.pR);
        tau = (d1 + d2) / c

        v_rel = dot(vk, (params.pR - pk)) / norm(params.pR - pk);
        nu = 2 * v_rel / lambda
    
        rcs = 10^(params.rcs_dB(k)/10);
        alpha = lambda * sqrt(rcs) / ((4*pi)^(3/2) * d1 * d2);
    
        b_tau = exp(-1j * 2*pi * (0:N-1)' * params.delta_f * tau);
        c_nu = exp(1j * 2*pi * (0:M-1) * params.Tsym * nu);
    
        H_true = H_true + alpha * (b_tau * c_nu);
    end
    
    % Add LOS path (k=0)
    d0 = norm(params.pT - params.pR);
    tau0 = d0 / c;
    alpha0 = lambda / (4*pi*d0)
    b0 = exp(-1j * 2*pi * (0:N-1)' * params.delta_f * tau0);
    c0 = ones(1, M);
    
    H_true = H_true + alpha0 * (b0 * c0);
    
    Y = X .* H_true + Z;

end
