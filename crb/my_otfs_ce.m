function [tau_hat, nu_hat] = my_otfs_ce(Y_TF, deltaf, T)
    [N, M] = size(Y_TF);
    
    delta_tau = 1/(N*deltaf);
    delta_nu = 1/(M*T);
    % 1. take ifft/fft with N_fft, M_fft    
    N_fft = 2^12;
    M_fft = 2^12;

    temp = sqrt(N_fft) * ifft(Y_TF, N_fft, 1);
    Y_DD = 1/sqrt(M_fft) * fft(temp, M_fft, 2);

    % 2. get the maximum point
    [~, linIdx] = max(abs(Y_DD(:)));
    [rowIdx, colIdx] = ind2sub(size(Y_DD), linIdx);
    
    % 3. idx / N/fft * N, idx / M/fft * M,
    rowIdx = (rowIdx - 1)/N_fft*N + 1;
    colIdx = (colIdx - 1)/M_fft*M + 1;

    % 4. remove threshold
    rowIdx = rowIdx - 35;
    colIdx = colIdx - 25;
    
    % 5. multiply with delta_tau, delta_nu.
    tau_hat = rowIdx * delta_tau;
    nu_hat = colIdx * delta_nu;
end