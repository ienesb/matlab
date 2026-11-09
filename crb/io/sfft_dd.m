function y = sfft_dd(Y)
% SFFT_DD  Unitary symplectic FFT, time-frequency Y[n,m] -> delay-Doppler y[l',k'].
% latex/io/cp_otfs_io.tex:
%   y[l',k'] = (1/sqrt(NM)) sum_{n,m} Y[n,m] exp( -j2pi (mk'/M - nl'/N) )
% Exact inverse of isfft_dd.  (Identical to the legacy sfft_u.)
    [N, M] = size(Y);
    y = sqrt(N/M) * ifft(fft(Y, [], 2), [], 1);
end
