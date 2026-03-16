function [est_tau, est_nu, est_alpha] = est_parameters(Y_DD_guard, data_mask_ce, l_p_ce, k_p_ce, deltaf, T, x_p)
    N1 = 41;
    M1 = 128;
    N2 = N1 * 10;
    M2 = M1 * 100;    

    Y = Y_DD_guard(~data_mask_ce);
    Y = reshape(Y, [N1, M1]);
    
    y = ifft2(Y, N1, M1) * sqrt(N1*M1);
    
    winN = hamming(N1);
    winM = hamming(M1);
    window2D = winN * winM'; 
    
    y_windowed = y .* window2D;
    
    Y_interp = fft2(y_windowed, N2, M2) / sqrt(N2*M2);
      
    th = 0.03;
    
    Y_interp_th = Y_interp;
    Y_interp_th(abs(Y_interp) < th) = 0;
    
    
    % TF = islocalmax(Y_interp_th, 1) & islocalmax(Y_interp_th, 2);
    TF = imregionalmax(abs(Y_interp_th));
    
    peakValues = abs(Y_interp_th(TF));
    linearIdx = find(TF);
    
    [~, sortOrder] = sort(peakValues, 'descend');
    
    sortedLinearIdx = linearIdx(sortOrder);
    
    % est_alpha = Y_interp_th(sortedLinearIdx) * sqrt(1000);
    
    [peak_rows, peak_cols] = ind2sub(size(Y_interp_th), sortedLinearIdx);
        
    est_tau = (peak_rows - l_p_ce*10 - 1) / (N2*deltaf);
    est_nu = (peak_cols - k_p_ce*100 - 1) / (M2*T);

    % Estimate complex path gain for each detected path.
    % Read from the oversampled windowed spectrum at the detected peak position.
    % At the true (fractional-Doppler) peak, the window phase is zero, so
    % the complex value carries the correct alpha phase with no fractional-
    % Doppler-induced rotation.
    W_N = sum(hamming(N1));
    W_M = sum(hamming(M1));
    norm_factor = sqrt(N1 * M1 * N2 * M2) / (W_N * W_M);
    est_alpha = Y_interp(sortedLinearIdx) .* norm_factor / x_p;
end