function Xhat = ofdm_banded_mmse(Y, deltaf, Ts, alpha, tau, nu, N0, W)
% OFDM_BANDED_MMSE  Banded MMSE equalisation of OFDM with the KNOWN ICI channel.
% Per OFDM symbol m the received subcarrier vector is y_m = H_m x_m + w, with the
% banded ICI matrix
%   H_m[n, n'] = sum_p alpha_p e^{j2pi m Ts nu_p} e^{-j2pi n' deltaf tau_p} I_p(n'-n),
% I_p from get_ici_coeffs, kept for |n'-n| <= W (the Dirichlet ICI decays fast).
% Returns the equalised (soft) symbol grid Xhat (caller demaps / hard-decides).
    [N, M] = size(Y);  P = numel(alpha);  T = 1/deltaf;
    Deltas = (-W:W).';  nd = numel(Deltas);
    ns = (0:N-1).';
    Xhat = zeros(N, M);

    % Per-path, per-Delta base column b_{p,Delta}[n] = alpha_p I_p(Delta) e^{-j2pi(n+Delta)deltaf tau_p}
    % (the only m-dependence is the scalar Doppler factor e^{j2pi m Ts nu_p}).
    base = zeros(N, nd, P);
    for p = 1:P
        for d = 1:nd
            D = Deltas(d);
            Ip = get_ici_coeffs(D, nu(p), T);
            base(:,d,p) = alpha(p) * Ip * exp(-1j*2*pi*(ns+D)*deltaf*tau(p));
        end
    end

    for m = 0:M-1
        dopp = exp(1j*2*pi*m*Ts*nu);                 % 1 x P
        ii = zeros(N*nd,1); jj = zeros(N*nd,1); vv = zeros(N*nd,1);  cnt = 0;
        for d = 1:nd
            D = Deltas(d);
            np = ns + D;
            ok = (np >= 0) & (np <= N-1);
            nn = ns(ok);  npv = np(ok);
            val = zeros(numel(nn),1);
            for p = 1:P
                val = val + base(ok,d,p) * dopp(p);
            end
            r = cnt + (1:numel(nn));
            ii(r) = nn + 1;  jj(r) = npv + 1;  vv(r) = val;  cnt = cnt + numel(nn);
        end
        Hm = sparse(ii(1:cnt), jj(1:cnt), vv(1:cnt), N, N);
        Xhat(:,m+1) = (Hm'*Hm + N0*speye(N)) \ (Hm' * Y(:,m+1));
    end
end
