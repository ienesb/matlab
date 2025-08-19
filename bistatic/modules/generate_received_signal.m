function Y = generate_received_signal(X, H, params)
    N = params.N;
    M = params.M;
    
    Z = sqrt(params.sigma2/2) * (randn(N, M) + 1j*randn(N, M));
    Y = X .* H + Z;
    % Y = X .* H;
end