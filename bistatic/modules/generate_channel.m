function H = generate_channel(params)
    v2struct(params);
    H = zeros(N, M);

    for k = 1:K
        pk = targets(k, :);
        vk = velocities(k, :);

        d1 = norm(pT - pk);
        d2 = norm(pk - pR);
        tau = (d1 + d2) / c;
        tau = taus(k);

        v_rel = dot(vk, (pR - pk)) / norm(pR - pk);
        nu = 2 * v_rel / lambda;
        nu = nus(k);

        rcs = 10^(rcs_dB(k)/10);
        alpha = lambda * sqrt(rcs) / ((4*pi)^(3/2) * d1 * d2);
        alpha = alphas(k);

        b_tau = getb(tau, params);
        c_nu = getc(nu, params);

        H = H + alpha * (b_tau * c_nu.');
    end
    
    % Add LOS path (k=0)
    d0 = norm(pT - pR);
    tau0 = d0 / c;
    nu0 = 0;
    alpha0 = lambda / (4*pi*d0);
    b0 = getb(tau0, params);
    c0 = getc(nu0, params);
    
    % H = H + alpha0 * (b0 * c0.');
end
