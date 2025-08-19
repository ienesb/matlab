function HDD = generate_dd_map(H, params)
    % N = params.N;
    % M = params.M;
    N = 4096;
    M = 256;

    Fn = dftmtx(N)./sqrt(N);
    Fm = dftmtx(M)./sqrt(M);
    
    % Delay-Doppler transformation: FFT over frequency, IFFT over time
    H_padded = zeros(N, M);
    H_padded(1:size(H, 1), 1:size(H, 2)) = H;
    DD_map_complex = Fn' * H_padded * Fm;
    HDD = abs(DD_map_complex).^2;
end