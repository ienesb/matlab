function ser = getSer(X, X_hat, data_mask)
    errors = X ~= X_hat;
    ser = sum(errors, "all") / sum(data_mask, "all");
end