function CRB = getCRB_OTFS_ICI (S, N, M, deltaf, T, sigma2, alpha, nu, tau) % monostatic [alpha_R, alpha_I, nu, tau]
    % F = dftmtx(N*M) / sqrt(N*M);

    k = (0:(N*M-1)).';
    b = exp(-1j*2*pi*deltaf*tau/M*k);
    c = exp(1j*2*pi*T*nu/N*k);
    
    B = diag(b);
    C = diag(c);
    
    % Y = alpha * C * F' * B * F * S;
    
    temp = fft(S, N*M, 1) / sqrt(N*M);
    temp = B * temp;
    
    temp2 = (-1j * 2 * pi * deltaf / M * k) .* temp;
    
    temp2 = ifft(temp2, N*M, 1) * sqrt(N*M);
    temp = ifft(temp, N*M, 1) * sqrt(N*M);


    g0 = C * temp;
    g1 = 1j * g0;
    g2 = alpha * (1j * 2 * pi * T / N * k) .* g0;
    g3 = alpha * C * temp2;
    

    G = [g0 g1 g2 g3];
    FIM = 2/sigma2 * real(G' * G);
    CRB = inv(FIM);
end