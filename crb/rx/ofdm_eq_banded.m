function Xsoft = ofdm_eq_banded(Y, g, cr, N0, W)
% OFDM_EQ_BANDED  Genie banded-ICI MMSE equalizer with the KNOWN CP-OFDM channel
% (latex/io/cp_ofdm_io.tex): per symbol m, y_m = H_m x_m + w with
%   H_m[n,n'] = sum_p alpha_p e^{j2pi m TS nu_p} e^{-j2pi n' df tau_p} I_p(n'-n)
% kept for |n'-n| <= W (default 8).  Port of old/eva/utils/ofdm_banded_mmse.m
% with ici_kernel and TS.
    if nargin < 5 || isempty(W), W = 8; end
    [N, M] = size(Y);  P = numel(cr.alpha);
    Deltas = (-W:W).';  nd = numel(Deltas);
    ns = (0:N-1).';
    Xsoft = zeros(N, M);
    base = zeros(N, nd, P);
    for p = 1:P
        Ip = ici_kernel(Deltas, cr.nu(p), g.T);                       % I_p(Delta), Delta = n' - n
        for d = 1:nd
            base(:, d, p) = cr.alpha(p) * Ip(d) * exp(-1j*2*pi*(ns + Deltas(d))*g.deltaf*cr.tau(p));
        end
    end
    for m = 0:M-1
        dopp = exp(1j*2*pi*m*g.TS*cr.nu(:)).';                          % 1 x P
        ii = zeros(N*nd, 1);  jj = ii;  vv = ii;  cnt = 0;
        for d = 1:nd
            np = ns + Deltas(d);
            ok = (np >= 0) & (np <= N-1);
            nn = ns(ok);  npv = np(ok);
            val = zeros(numel(nn), 1);
            for p = 1:P, val = val + base(ok, d, p) * dopp(p); end
            r = cnt + (1:numel(nn));
            ii(r) = nn + 1;  jj(r) = npv + 1;  vv(r) = val;  cnt = cnt + numel(nn);
        end
        Hm = sparse(ii(1:cnt), jj(1:cnt), vv(1:cnt), N, N);
        Xsoft(:, m+1) = (Hm'*Hm + N0*speye(N)) \ (Hm' * Y(:, m+1));
    end
end
