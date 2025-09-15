function H_hat = estimate_channel(X, Y, H_hat, mask, type, SNR)
    
    switch type
        case "RF"
            H_hat(mask) = Y(mask) ./ X(mask); % Reciprocal Filtering (RF)
        case "MF"
            H_hat(mask) = Y(mask) .* conj(X(mask)); % Matched Filtering (MF)
        case "LMMSE"
            % SNR = 0.3334;
            H_hat(mask) = ( conj(X(mask)) .* Y(mask) ) ./ ( abs(X(mask)).^2 + 1/SNR ); % LMMSE estimator
    end
end
