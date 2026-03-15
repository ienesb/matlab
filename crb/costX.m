function J = costX(eta, tau, nu, alpha, deltaf, Ts, Y) % eta: X
    X = eta;
    Y_hat = channel_noiseless(X, alpha, tau, nu, deltaf, Ts);
    J = sum(abs(Y(:) - Y_hat(:)).^2);
end