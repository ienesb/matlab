function [ser, rmse] = getSer(X, X_hat, pilot_mask, bPlot)
    if nargin < 4
        bPlot = 0;
    end
    data_mask = ~pilot_mask;
    rmse = mean(abs(X_hat(:) - X(:)).^2);
    errors = X ~= X_hat;
    ser = sum(errors(:)) / sum(data_mask(:));
    if bPlot
        scatterplot(X_hat(:));
    end
end