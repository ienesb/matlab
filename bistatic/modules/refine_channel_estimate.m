function H_refined = refine_channel_estimate(X_hat, Y, pilot_mask)
    % Refine the channel estimate using data symbols via matched filtering
    data_mask = ~pilot_mask;

    H_refined = initial_channel_estimation(X_hat, Y, pilot_mask);

    SNR = 0.244;
    H_refined(data_mask) = ( conj(X_hat(data_mask)) .* Y(data_mask) ) ./ ( abs(X_hat(data_mask)).^2 + 1/SNR ); % LMMSE estimator
end
