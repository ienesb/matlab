function HDD = generate_dd_map(H, params)
    N = params.N;
    M = params.M;

    Fn = dftmtx(N)./sqrt(N);
    Fm = dftmtx(M)./sqrt(M);
    
    % Delay-Doppler transformation: FFT over frequency, IFFT over time
    DD_map_complex = Fn' * H * Fm;
    HDD = abs(DD_map_complex).^2;
end