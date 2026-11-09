function J = getFIM_OTFS_multi(x, deltaf, Ts, sigma2, alphas, taus, nus)
% Multi-target FIM for OTFS (Section 3.2 of otfs_fim.tex).
% Parameter stacking: eta = [eta_1; ...; eta_P], eta_p = [alphaR_p, alphaI_p, tau_p, nu_p].
%
% Inputs:
%   x       : N x M delay-Doppler symbol grid
%   deltaf  : subcarrier spacing [Hz]
%   Ts      : OFDM symbol duration (with CP) [s]
%   sigma2  : noise variance
%   alphas  : P x 1 complex path gains
%   taus    : P x 1 path delays [s]
%   nus     : P x 1 path Dopplers [Hz]
%
% Output:
%   J       : 4P x 4P Fisher information matrix.

    [N, M] = size(x);
    P = numel(alphas);

    ns = (0:N-1).';
    ms = (0:M-1).';
    ds_A = (-(N-1):(N-1)).';
    ds_B = (-(M-1):(M-1)).';

    G = zeros(N*M, 4*P);

    for p = 1:P
        tau_p = taus(p);
        nu_p  = nus(p);
        a_p   = alphas(p);

        % Per-path delay/Doppler kernels A_p(d), A'_p(d), B_p(d), B'_p(d)
        E_A  = exp(1j*2*pi*ns*(ds_A.'/N - deltaf*tau_p));    % N x (2N-1)
        A_d  = sum(E_A, 1).';
        Ap_d = sum(ns .* E_A, 1).';

        E_B  = exp(1j*2*pi*ms*(ds_B.'/M + Ts*nu_p));         % M x (2M-1)
        B_d  = sum(E_B, 1).';
        Bp_d = sum(ms .* E_B, 1).';

        A_mat  = toeplitz(A_d (N:end), A_d (N:-1:1));
        Ap_mat = toeplitz(Ap_d(N:end), Ap_d(N:-1:1));
        B_mat  = toeplitz(B_d (M:end), B_d (M:-1:1));
        Bp_mat = toeplitz(Bp_d(M:end), Bp_d(M:-1:1));

        Phi_p   = (A_mat  * x * B_mat ) / (N*M);
        Psi_p   = (Ap_mat * x * B_mat ) / (N*M);
        Omega_p = (A_mat  * x * Bp_mat) / (N*M);

        g_aR  = Phi_p;
        g_aI  = 1j * Phi_p;
        g_tau = -1j*2*pi*deltaf * a_p * Psi_p;
        g_nu  =  1j*2*pi*Ts     * a_p * Omega_p;

        G(:, 4*(p-1)+1) = g_aR(:);
        G(:, 4*(p-1)+2) = g_aI(:);
        G(:, 4*(p-1)+3) = g_tau(:);
        G(:, 4*(p-1)+4) = g_nu(:);
    end

    J = (2/sigma2) * real(G' * G);
end
