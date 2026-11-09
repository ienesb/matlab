function [y, s, r] = rcp_otfs_time_sim(x, g, cr)
% RCP_OTFS_TIME_SIM  Brute-force sampled RCP-OTFS chain (oracle for io_rcp_otfs),
% latex/io/rcp_otfs_io.tex eqs. RCP_frame, RCP_time, RCP_circ, RCP_rx:
%   s[mN+q] = (1/sqrt(M)) sum_k x[q,k] e^{j2pi mk/M}        (unitary-normalized samples)
%   r[c]    = sum_p alpha_p e^{j2pi nu_p c Ts} s[(c - l_p) mod NM]   (frame-circular)
%   y[l',k']= (1/sqrt(M)) sum_m r[mN+l'] e^{-j2pi m k'/M}
% Exact in Doppler and in ISI for sample-spaced delays tau_p = l_p Ts.
    [N, M] = size(x);
    S = sqrt(M) * ifft(x, [], 2);                 % N x M, S(q+1, m+1) = s[mN+q]
    s = S(:);                                     % column-major: c = mN + q
    c = (0:N*M-1).';
    r = zeros(N*M, 1);
    for p = 1:numel(cr.alpha)
        lp = round(cr.tau(p)/g.Ts);
        assert(abs(cr.tau(p)/g.Ts - lp) < 1e-6, 'rcp_otfs_time_sim:offgrid', 'integer delays only');
        r = r + cr.alpha(p) * exp(1j*2*pi*cr.nu(p)*c*g.Ts) .* circshift(s, lp);
    end
    R = reshape(r, N, M);                         % R(l'+1, m+1) = r[mN+l']
    y = fft(R, [], 2) / sqrt(M);
end
