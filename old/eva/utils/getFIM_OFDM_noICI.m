function J = getFIM_OFDM_noICI(X, deltaf, Ts, sigma2, alpha, tau, nu)
% Single-target OFDM Fisher information matrix in the no-ICI (orthogonal
% subcarrier) regime, Section 3.1.1 of ofdm_fim.tex.  Parameter ordering:
%   eta = [alphaR, alphaI, tau, nu].
%
% For unit-energy QPSK this is data-deterministic; for general X it is
% the realization-conditional FIM.

    [N, M] = size(X);

    ns = (0:N-1).';
    ms = (0:M-1);

    phase = exp(-1j*2*pi*deltaf*ns*tau) * exp(1j*2*pi*Ts*ms*nu);

    g_aR  = X .* phase;
    g_aI  = 1j * g_aR;
    g_tau = (-1j*2*pi*deltaf) * (ns .* g_aR) * alpha;
    g_nu  = ( 1j*2*pi*Ts    ) * (g_aR .* ms) * alpha;

    G = [g_aR(:), g_aI(:), g_tau(:), g_nu(:)];
    J = (2/sigma2) * real(G' * G);
end
