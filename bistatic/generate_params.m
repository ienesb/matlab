function [taus, nus, alphas] = generate_params(positions, velocities, pT, pR, lambda, rcs)
    taus = zeros(K, 1);
    nus = zeros(K, 1);
    alphas = zeros(K, 1);
    for k = 1:K
        pk = positions(k, :);
        vk = velocities(k, :);

        tau = pos2tau(pk, pT, pR);
        taus(k) = tau;

        % los_vector = (pR - pk) / norm(pR - pk);
        % v_rel = dot(vk, los_vector);
        % nu = 2 * v_rel / lambda;
        nu = velocity2nu(vk.', pk.', pT.', pR.', lambda);
        nus(k) = nu;

        alpha = lambda * sqrt(rcs(k)) / ((4*pi)^1.5 * norm(pk - pT) * norm(pk - pR));
        alphas(k) = alpha;
    end
end