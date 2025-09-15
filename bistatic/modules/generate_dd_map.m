function HDD = generate_dd_map(H, params) % N_fft, M_fft, is_fftshifted
    v2struct(params);

    % window2D = hamming(size(H, 1)) * (hamming(size(H, 2))).';
    % H = H .* window2D;

    % Step 1: IFFT over columns to get the delay profile
    H_temp = ifft(H, N_fft) * sqrt(N_fft); % N_fft-by-M
    
    % Step 2: FFT over rows to get the Doppler profile
    H_temp = fft(H_temp.', M_fft) * 1/sqrt(M_fft);

    HDD = H_temp.';
    HDD = abs(HDD).^2;
    
    if is_fftshifted
        HDD = fftshift(HDD, 2);
    end
end