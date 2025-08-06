function H = generate_channel(params)
    N = params.N;
    M = params.M;
    K = params.K;
    fc = params.fc;
    c = 3e8;
    lambda = c / fc;
    
    H = zeros(N, M);
    
    for k = 1:K
        pk = params.targets(k, :);
        vk = params.velocities(k, :);

        d1 = norm(params.pT - pk);
        d2 = norm(pk - params.pR);
        tau = (d1 + d2) / c;

        v_rel = dot(vk, (params.pR - pk)) / norm(params.pR - pk);
        nu = 2 * v_rel / lambda;

        rcs = 10^(params.rcs_dB(k)/10);
        alpha = lambda * sqrt(rcs) / ((4*pi)^(3/2) * d1 * d2);

        b_tau = getb(tau, params);
        c_nu = getc(nu, params);

        H = H + alpha * (b_tau * c_nu.');
    end
    
    % Add LOS path (k=0)
    d0 = norm(params.pT - params.pR);
    tau0 = d0 / c;
    nu0 = 0;
    alpha0 = lambda / (4*pi*d0);
    b0 = getb(tau0, params);
    c0 = getc(nu0, params);
    
    H = H + alpha0 * (b0 * c0.');
end
