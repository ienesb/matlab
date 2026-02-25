function kl = kl_distancev2(eta, deltaf, T, X_hat, mu_gt)
    nu  = eta(1);
    tau = eta(2);
    alphaR = eta(3);
    alphaI = eta(4);
    alpha = alphaR + 1j*alphaI;

    [N, M] = size(X_hat);
    ns = (0:(N-1)).';
    ms = 0:(M-1);

    c_tilde = (exp(-1j*2*pi*deltaf*ns*tau) * exp(1j*2*pi*T*ms*nu)) .* X_hat;
    c_tilde = c_tilde(:);
    
    mu_tilde = alpha * c_tilde;

    kl = norm(mu_gt - mu_tilde)^2;
end