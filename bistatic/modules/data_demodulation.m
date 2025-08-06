function X_hat = data_demodulation(Y, H_est, Xp, data_mask, params)
    % LMMSE data demodulation for data symbols
    
    N = params.N;
    M = params.M;
    X_hat = zeros(N, M);
    
    H_data = H_est(data_mask);
    Y_data = Y(data_mask);
    
    X_hat(data_mask) = (Y_data .* conj(H_data)) ./ (abs(H_data).^2 + params.sigma2);
    X_hat = X_hat + Xp;
end
