function X_hat = data_demodulation(Y, H, X, pilot_mask, params) % sigma2
    sigma2 = params.sigma2;

    data_mask = ~pilot_mask;
    X_hat = X;
    
    H_data = H(data_mask);
    Y_data = Y(data_mask);
    
    X_hat(data_mask) = (Y_data .* conj(H_data)) ./ (abs(H_data).^2 + sigma2);
    
    % data = qamdemod(X_hat(data_mask), 4);
    % X_hat(data_mask) = qammod(data, 4, 'UnitAveragePower', true);
end
