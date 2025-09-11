function [X, pilot_mask] = generate_ofdm_symbols(params) % N, M, M_order, pilot_ratio
    v2struct(params);

    pilot_mask = false(N, M);
    
    num_pilots = round(pilot_ratio * N * M);
    pilot_indices = randperm(N * M, num_pilots);
    pilot_mask(pilot_indices) = true;
    
    Xp = zeros(N, M);
    Xp(pilot_mask) = 1;
    
    Xd = zeros(N, M);
    data_mask = ~pilot_mask;
    
    data_symbols = qammod(randi([0 M_order-1], sum(data_mask(:)), 1), M_order, 'UnitAveragePower', true);
    Xd(data_mask) = data_symbols;
    
    X = Xp + Xd;
end
