function kl = kl_distance(eta, deltaf, T, X_hat, mu_gt)
    nu  = eta(1);
    tau = eta(2);

    [N, M] = size(X_hat);
    ns = (0:(N-1)).';
    ms = 0:(M-1);


    c_tilde = (exp(-1j*2*pi*deltaf*ns*tau) * exp(1j*2*pi*T*ms*nu)) .* X_hat;
    c_tilde = c_tilde(:);

    alpha = pinv(c_tilde) * mu_gt;
    
    mu_tilde = alpha * c_tilde;

    kl = norm(mu_gt - mu_tilde)^2;
end