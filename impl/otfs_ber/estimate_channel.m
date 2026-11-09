function Ghat = estimate_channel(Y, p, xp, thr)
% ESTIMATE_CHANNEL  Threshold-based embedded-pilot channel estimation for
% the fractional-Doppler FULL-guard case (Embedded Sec. IV-B, eq. 9).
%
% For each received symbol in the pilot delay band l in [l_p, l_p+l_tau]
% and every Doppler index k:
%   if |y[k,l]| >= thr :  Ghat[[k-k_p]_N, l-l_p] = y[k,l]/x_p
%   else               :  Ghat[...] = 0
% The guard band guarantees these received symbols are free of data
% interference, so Ghat is a direct (sparse) estimate of the effective
% DD response G used by the detector.
%
% Y : N x M received DD grid.  xp : pilot amplitude.  thr : detection
% threshold on |y| (= 3*sigma).  Returns Ghat (N x M).

N = p.N; M = p.M;
Ghat = zeros(N, M);

for lp = 0:p.l_tau                       % effective delay tap l'
    l   = p.l_p + lp;                    % grid delay column
    col = Y(:, l+1);                     % received pilot column (Nx1)
    det = abs(col) >= thr;               % threshold detection
    if any(det)
        % Map grid Doppler k -> effective Doppler tap k' = (k-k_p) mod N
        kpr = mod((0:N-1).' - p.k_p, N); % k' for each k
        idx = find(det);
        Ghat(kpr(idx)+1, lp+1) = col(idx)/xp;
    end
end
end
