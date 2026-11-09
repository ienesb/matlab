function J = fim_cp_otfs_direct(x, g, cr, use_ici, sigma2)
% FIM_CP_OTFS_DIRECT  Independent delay-Doppler-domain CP-OTFS FIM built from the
% Dirichlet kernels A_p, B_p (no ICI; cp_otfs_fim.tex eqs. 2, 4, 9-11) and the
% two-index kernels K_p, K_p^{(1)}, K_p^{I'} (ICI; eqs. 43-45, 50-52).  Test
% oracle only: O(N^3) per target, used to validate the ISFFT-unitarity route.
    [N, M] = size(x);
    P = numel(cr.alpha);
    G = zeros(N*M, 4*P);
    ns = (0:N-1).';  ms = (0:M-1).';
    [Lp, L] = ndgrid(0:N-1, 0:N-1);   % [A]_{l',l} = A(l'-l)
    [K, Kp] = ndgrid(0:M-1, 0:M-1);   % [B]_{k,k'} = B(k-k')
    for p = 1:P
        a = cr.alpha(p);  tau = cr.tau(p);  nu = cr.nu(p);
        Bm  = arrayfun(@(d) sum(       exp(1j*2*pi*ms*(d/M + g.TS*nu))), K - Kp);
        Bpm = arrayfun(@(d) sum(ms  .* exp(1j*2*pi*ms*(d/M + g.TS*nu))), K - Kp);
        if ~use_ici
            Am  = arrayfun(@(d) sum(       exp(1j*2*pi*ns*(d/N - g.deltaf*tau))), Lp - L);
            Apm = arrayfun(@(d) sum(ns  .* exp(1j*2*pi*ns*(d/N - g.deltaf*tau))), Lp - L);
            Phi = Am  * x * Bm  / (N*M);
            Psi = Apm * x * Bm  / (N*M);
            Om  = Am  * x * Bpm / (N*M);
            gnu = 1j*2*pi*g.TS*a * Om;
        else
            [NP, NN] = ndgrid(0:N-1, 0:N-1);
            [I, Ip] = ici_kernel(NP - NN, nu, g.T);           % I(n'+1, n+1) = I(n'-n)
            V  = exp(1j*2*pi*(0:N-1).' * (0:N-1) / N);        % V(n, l') = e^{j2pi n l'/N}
            U  = exp(-1j*2*pi*g.deltaf*ns*tau) .* exp(-1j*2*pi*ns*(0:N-1)/N);  % U(n', l)
            % K(l',l) = sum_{n,n'} I(n'-n) e^{-j2pi n'[df tau + l/N]} e^{j2pi n l'/N}
            %         = sum_n V(n,l') sum_{n'} I(n'-n) U(n',l) = (V.' * (I.' * U))(l', l)
            Km  = V.' * (I.'  * U);
            K1  = V.' * (I.'  * (ns .* U));
            KIp = V.' * (Ip.' * U);
            Phi = Km  * x * Bm / (N*M);
            Psi = K1  * x * Bm / (N*M);
            Lam = (g.TS * (Km * x * Bpm) + KIp * x * Bm) / (N*M);
            gnu = 1j*2*pi*a * Lam;
        end
        col = 4*(p-1);
        G(:, col+1) = Phi(:);
        G(:, col+2) = 1j*Phi(:);
        G(:, col+3) = -1j*2*pi*g.deltaf*a * Psi(:);
        G(:, col+4) = gnu(:);
    end
    J = (2/sigma2) * real(G'*G);
end
