function HDD = generate_dd_map(H, params)
    v2struct(params);

    % Fn = dftmtx(N_fft)./sqrt(N_fft);
    % Fm = dftmtx(M_fft)./sqrt(M_fft);
    
    % % Delay-Doppler transformation: FFT over frequency, IFFT over time
    % H_padded = zeros(N_fft, M_fft);
    % H_padded(1:size(H, 1), 1:size(H, 2)) = H;
    % DD_map_complex = Fn' * H_padded * Fm;
    % HDD = abs(DD_map_complex).^2;

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