function H_est = initial_channel_estimation(X, Y, pilot_mask)
    % Initial pilot-only channel estimate using reciprocal filtering
    H_est = zeros(size(Y));

    % H_est(pilot_mask) = Y(pilot_mask) ./ X(pilot_mask); % Reciprocal Filtering (RF)
    
    % H_est(pilot_mask) = Y(pilot_mask) .* conj(X(pilot_mask)); % Matched Filtering (MF)
    
    SNR = 0.244;
    H_est(pilot_mask) = ( conj(X(pilot_mask)) .* Y(pilot_mask) ) ./ ( abs(X(pilot_mask)).^2 + 1/SNR ); % LMMSE estimator
end
