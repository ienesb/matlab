function H_est = initial_channel_estimation(X, Y, pilot_mask)
    % Initial pilot-only channel estimate using reciprocal filtering
    H_est = zeros(size(Y));
    H_est(pilot_mask) = Y(pilot_mask) ./ X(pilot_mask);
end
