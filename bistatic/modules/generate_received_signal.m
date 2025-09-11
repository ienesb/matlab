function Y = generate_received_signal(X, H, params, is_noisy) % sigma2 N, M
    if nargin < 4 
        is_noisy = 1;
    end
    v2struct(params);
 
    Y = X .* H;
    if is_noisy
        Z = sqrt(sigma2/2) * (randn(N, M) + 1j*randn(N, M));
        Y = Y + Z;
    end
end