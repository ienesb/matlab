function perm = match_paths(tau_true, nu_true, tau_hat, nu_hat)
% Optimal one-to-one assignment between true and estimated paths,
% minimizing the sum of normalized (tau, nu) squared distances.
% Returns perm such that tau_hat(perm(i)) pairs with tau_true(i).

    tau_true = tau_true(:);
    nu_true  = nu_true(:);
    tau_hat  = tau_hat(:);
    nu_hat   = nu_hat(:);

    P = length(tau_true);

    tau_scale = max(max(abs(tau_true)), 1e-12);
    nu_scale  = max(max(abs(nu_true)),  1);

    cost = ((tau_true - tau_hat.') / tau_scale).^2 + ...
           ((nu_true  - nu_hat.' ) / nu_scale ).^2;

    pairs = matchpairs(cost, 1e12);   % Kx2: [true_idx, est_idx]

    perm = zeros(P, 1);
    for k = 1:size(pairs, 1)
        perm(pairs(k, 1)) = pairs(k, 2);
    end
end
