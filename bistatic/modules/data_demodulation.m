function X_hat = data_demodulation(Y, H_est, Xp, data_mask, params)
    % LMMSE data demodulation for data symbols
    v2struct(params);

    X_hat = zeros(N, M);
    
    H_data = H_est(data_mask);
    Y_data = Y(data_mask);
    
    X_hat(data_mask) = (Y_data .* conj(H_data)) ./ (abs(H_data).^2 + sigma2);
    
    X_hat = X_hat + Xp;

    % data = qamdemod(X_hat(data_mask), 4);
    % X_hat(data_mask) = qammod(data, 4, 'UnitAveragePower', true);
end
