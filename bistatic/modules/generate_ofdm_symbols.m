function [X, Xp, Xd, pilot_mask, data_mask] = generate_ofdm_symbols(params)
    % Generate OFDM symbols with pilots and data, and separate outputs
    
    NM = params.N * params.M;
    pilot_mask = false(params.N, params.M);
    
    % Randomly assign pilots
    num_pilots = round(params.pilot_ratio * NM);
    pilot_indices = randperm(NM, num_pilots);
    pilot_mask(pilot_indices) = true;
    
    % Assign pilot symbols
    Xp = zeros(params.N, params.M);
    Xp(pilot_mask) = 1;
    
    % Assign data symbols
    Xd = zeros(params.N, params.M);
    data_mask = ~pilot_mask;
    
    switch params.modulation
        case 'QPSK'
            M_order = 4;
        case '16QAM'
            M_order = 16;
        otherwise
            M_order = 4;
    end
    
    data_symbols = qammod(randi([0 M_order-1], sum(data_mask(:)), 1), M_order, 'UnitAveragePower', true);
    Xd(data_mask) = data_symbols;
    
    % Combine to form full transmit matrix
    X = Xp + Xd;
end
