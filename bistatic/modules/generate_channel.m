function H = generate_channel(params) % N, M, taus, nus, alphas, K, 
    v2struct(params);
    H = zeros(N, M);

    for k = 1:K
        tau = taus(k);
        nu = nus(k);
        alpha = alphas(k);

        b_tau = getb(tau, params);
        c_nu = getc(nu, params);

        H = H + alpha * (b_tau * c_nu.');
    end

    % Add LOS path (k=0)
    % d0 = norm(pT - pR);
    % tau0 = d0 / c;
    % nu0 = 0;
    % alpha0 = lambda / (4*pi*d0);
    % b0 = getb(tau0, params);
    % c0 = getc(nu0, params);
    % H = H + alpha0 * (b0 * c0.');
end
