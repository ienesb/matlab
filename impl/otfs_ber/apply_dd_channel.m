function Y = apply_dd_channel(X, G, sigma2)
% APPLY_DD_CHANNEL  Apply the exact 2D delay-Doppler channel G to DD-domain
% symbols X and add AWGN.  Implements
%   Y[k,l] = sum_{k',l'} G[k',l'] X[[k-k']_N, [l-l']_M] + v[k,l]
% efficiently as one Doppler-domain circular convolution per occupied delay
% tap (FFT along the Doppler dimension).
%
% X, G : N x M.   sigma2 : noise variance (per complex sample).

[N, M] = size(X);
Y = zeros(N, M);
cols = find(any(G ~= 0, 1));          % occupied delay taps l' (1-based)

for c = cols
    lp = c - 1;                       % delay tap l'
    Xs = circshift(X, [0, lp]);       % Xs(:,l) = X(:, l-l')
    Y  = Y + ifft(fft(G(:,c)) .* fft(Xs, [], 1), [], 1);
end

% AWGN, v ~ CN(0, sigma2)
Y = Y + sqrt(sigma2/2)*(randn(N,M) + 1i*randn(N,M));
end
