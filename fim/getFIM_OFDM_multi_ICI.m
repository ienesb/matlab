function J = getFIM_OFDM_multi_ICI(X, deltaf, Ts, sigma2, alphas, taus, nus)
% Multi-target FIM for OFDM, ICI channel (Section 3.2.2 of ofdm_fim.tex).
% Parameter stacking: eta = [eta_1; ...; eta_P], eta_p = [alphaR_p, alphaI_p, tau_p, nu_p].
%
% Uses the dimensionless ICI kernel of crb/EVA/get_ici_coeffs.m which satisfies
% I(0)=1 at nu=0, recovering the no-ICI Kronecker-delta limit.
%
% Inputs / output have the same shape as getFIM_OFDM_multi_noICI.

    [N, M] = size(X);
    P = numel(alphas);
    T = 1 / deltaf;

    ns = (0:N-1).';
    ms = (0:M-1);

    G = zeros(N*M, 4*P);

    for p = 1:P
        nu_p  = nus(p);
        tau_p = taus(p);
        a_p   = alphas(p);

        A = X .* exp(-1j*2*pi*deltaf*ns*tau_p) .* exp(1j*2*pi*Ts*ms*nu_p);   % N x M

        Phi   = zeros(N, M);
        Psi   = zeros(N, M);
        Omega = zeros(N, M);

        for n = 0:N-1
            I_p  = get_ici_coeffs           (ns - n, nu_p, T);   % N x 1
            dI_p = get_ici_coeffs_derivative(ns - n, nu_p, T);   % N x 1

            Phi  (n+1, :) = sum(A .* I_p, 1);
            Psi  (n+1, :) = sum(A .* (ns*deltaf) .* I_p, 1);
            Omega(n+1, :) = sum(A .* (ms*Ts .* I_p + dI_p), 1);
        end

        g_aR  = Phi;
        g_aI  = 1j * Phi;
        g_tau = -1j*2*pi*a_p * Psi;
        g_nu  =  1j*2*pi*a_p * Omega;

        G(:, 4*(p-1)+1) = g_aR(:);
        G(:, 4*(p-1)+2) = g_aI(:);
        G(:, 4*(p-1)+3) = g_tau(:);
        G(:, 4*(p-1)+4) = g_nu(:);
    end

    J = (2/sigma2) * real(G' * G);
end
