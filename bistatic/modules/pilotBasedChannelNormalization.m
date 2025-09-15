function H = pilotBasedChannelNormalization(H, pilot_mask)
    data_mask = ~pilot_mask;

    Hp = H(pilot_mask);
    Hd = H(data_mask);
    avgPilotPower = sum(abs(Hp(:)).^2) / sum(pilot_mask(:));
    avgDataPower = sum(abs(Hd(:)).^2) / sum(data_mask(:));
    
    if avgDataPower ~= 0
        H(data_mask) = H(data_mask) * sqrt( avgPilotPower / avgDataPower );
    end
end