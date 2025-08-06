function H_refined = refine_channel_estimate(X_hat, Y, pilot_mask)
    % Refine the channel estimate using data symbols via matched filtering
    data_mask = ~pilot_mask;

    H_refined = zeros(size(Y));
    H_refined(pilot_mask) = Y(pilot_mask) ./ X_hat(pilot_mask);
    H_refined(data_mask) = Y(data_mask) .* conj(X_hat(data_mask));
end
