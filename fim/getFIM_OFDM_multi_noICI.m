function J = getFIM_OFDM_multi_noICI(X, deltaf, Ts, sigma2, alphas, taus, nus)
% Multi-target FIM for OFDM, no-ICI channel (Section 3.2.1 of ofdm_fim.tex).
% Parameter stacking: eta = [eta_1; ...; eta_P], eta_p = [alphaR_p, alphaI_p, tau_p, nu_p].
%
% Inputs:
%   X       : N x M transmitted symbol grid
%   deltaf  : subcarrier spacing [Hz]
%   Ts      : OFDM symbol duration (with CP) [s]
%   sigma2  : noise variance
%   alphas  : P x 1 complex path gains
%   taus    : P x 1 path delays  [s]
%   nus     : P x 1 path Dopplers [Hz]
%
% Output:
%   J       : 4P x 4P Fisher information matrix.

    [N, M] = size(X);
    P = numel(alphas);

    ns = (0:N-1).';
    ms = (0:M-1);

    G = zeros(N*M, 4*P);

    for p = 1:P
        phase = exp(-1j*2*pi*deltaf*ns*taus(p)) * exp(1j*2*pi*Ts*ms*nus(p));   % N x M
        Phi_p = X .* phase;                                                    % N x M

        g_aR  = Phi_p;
        g_aI  = 1j * Phi_p;
        g_tau = (-1j*2*pi*deltaf) * alphas(p) * (ns .* Phi_p);
        g_nu  = ( 1j*2*pi*Ts)     * alphas(p) * (Phi_p .* ms);

        G(:, 4*(p-1)+1) = g_aR(:);
        G(:, 4*(p-1)+2) = g_aI(:);
        G(:, 4*(p-1)+3) = g_tau(:);
        G(:, 4*(p-1)+4) = g_nu(:);
    end

    J = (2/sigma2) * real(G' * G);
end
