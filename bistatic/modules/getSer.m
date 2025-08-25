function ser = getSer(X, X_hat, data_mask)

    RMSE_dB = pow2db((mean(abs(X_hat(:) - X(:)).^2)))
    errors = X ~= X_hat;
    ser = sum(errors, "all") / sum(data_mask, "all");
end