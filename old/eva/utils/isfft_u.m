function X = isfft_u(x)
% Unitary ISFFT mapping a delay-Doppler grid x[l,k] (N x M) to a
% time-frequency grid X[n,m] (N x M), matching cp_otfs_io.tex eq. (1):
%   X[n,m] = (1/sqrt(NM)) sum_{l,k} x[l,k] exp(j 2pi (mk/M - nl/N))
%
% Verified inverse of sfft_u: SFFT_u(ISFFT_u(x)) = x.

    [N, M] = size(x);
    X = sqrt(M/N) * fft(ifft(x, [], 2), [], 1);
end
