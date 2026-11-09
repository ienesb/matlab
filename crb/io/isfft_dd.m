function X = isfft_dd(x)
% ISFFT_DD  Unitary inverse symplectic FFT, delay-Doppler x[l,k] -> time-frequency
% X[n,m], both N x M (delay/subcarrier rows, Doppler/symbol columns).
% latex/io/cp_otfs_io.tex eq. (1):
%   X[n,m] = (1/sqrt(NM)) sum_{l,k} x[l,k] exp( j2pi (mk/M - nl/N) )
% Exact inverse of sfft_dd.  (Identical to the legacy isfft_u.)
    [N, M] = size(x);
    X = sqrt(M/N) * fft(ifft(x, [], 2), [], 1);
end
