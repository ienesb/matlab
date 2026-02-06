function [RMSE, ber] = getSer(X, Xhat, pilot_ratio)
    if nargin < 3
        pilot_ratio = 0;
    end
    [N, M] = size(X);
    receivedBits = qamdemod(Xhat, 4, "OutputType", "bit");
    trueBits = qamdemod(X, 4, "OutputType", "bit");

    RMSE = sum(abs(Xhat(:) - X(:)).^2)/(N*M*(1-pilot_ratio));
    ber = sum(receivedBits ~= trueBits, "all")/(N*M*2*(1-pilot_ratio));
end