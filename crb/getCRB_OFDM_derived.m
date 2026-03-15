function [CRB_delay, CRB_range, CRB_doppler, CRB_velocity] = getCRB_OFDM_derived(N, M, np, mp, c, lambda, deltaf, T, beta, SNR_lins)
    K = ceil((N - 1) / np);
    L = ceil((M - 1) / mp);
    
    absP = (K + 1) * (L + 1);
    
    CRB_range = 12 / (K * (K + 2) * absP * np^2) * c^2 ./ (8 * pi^2 * SNR_lins * deltaf^2);
    CRB_range = sqrt(CRB_range);
    CRB_delay = CRB_range ./ c;
    
    CRB_velocity = 12 / (L * (L + 2) * absP * mp^2) * lambda^2 ./ (32 * pi^2 * SNR_lins * T^2 * cos(beta/2)^2);
    CRB_velocity = sqrt(CRB_velocity);
    CRB_doppler = CRB_velocity ./ (lambda / (2*cos(beta/2)));
end