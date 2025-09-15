function X = hard_demod(X, data_mask)
    data = qamdemod(X(data_mask), 4);
    X(data_mask) = qammod(data, 4, 'UnitAveragePower', true);
end