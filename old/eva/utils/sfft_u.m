function y = sfft_u(Y)
% Unitary SFFT mapping a time-frequency grid Y[n,m] (N x M) back to a
% delay-Doppler grid y[l',k'] (N x M), matching cp_otfs_io.tex:
%   y[l',k'] = (1/sqrt(NM)) sum_{n,m} Y[n,m] exp(-j 2pi (mk'/M - nl'/N))
%
% Verified inverse of isfft_u.

    [N, M] = size(Y);
    y = sqrt(N/M) * ifft(fft(Y, [], 2), [], 1);
end
