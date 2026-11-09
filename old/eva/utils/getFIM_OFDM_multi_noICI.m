function J = getFIM_OFDM_multi_noICI(X, deltaf, Ts, sigma2, alpha, tau, nu)
% Multi-target OFDM Fisher information matrix, no-ICI channel
% (Section 3.2.1 of ofdm_fim.tex).  Parameter ordering:
%   eta = [a1R, a1I, tau1, nu1,  a2R, a2I, tau2, nu2,  ...]
%
% Inputs / output mirror getFIM_OFDM_multi_ICI; J is 4P x 4P.

    [N, M] = size(X);
    P      = numel(alpha);

    ns = (0:N-1).';
    ms = (0:M-1);

    G = zeros(N*M, 4*P);
    for p = 1:P
        ap = alpha(p);  tp = tau(p);  np = nu(p);

        phase = exp(-1j*2*pi*deltaf*ns*tp) * exp(1j*2*pi*Ts*ms*np);
        g_aR  = X .* phase;
        g_tau = (-1j*2*pi*deltaf) * (ns .* g_aR) * ap;
        g_nu  = ( 1j*2*pi*Ts    ) * (g_aR .* ms) * ap;

        col = 4*(p-1);
        G(:, col+1) =       g_aR(:);
        G(:, col+2) = 1j  * g_aR(:);
        G(:, col+3) =       g_tau(:);
        G(:, col+4) =       g_nu(:);
    end

    J = (2/sigma2) * real(G' * G);
end
