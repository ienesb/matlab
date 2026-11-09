function J = getFIM_OFDM_multi_ICI(X, deltaf, Ts, sigma2, alpha, tau, nu)
% Multi-target OFDM Fisher information matrix with the full ICI kernel
% (Section 3.2.2 of ofdm_fim.tex).  Parameter ordering, stacking targets
% in order p = 1..P:
%   eta = [a1R, a1I, tau1, nu1,  a2R, a2I, tau2, nu2,  ...]
%
% Inputs:
%   X       : N x M transmitted symbol grid
%   alpha   : P x 1 complex path gains
%   tau     : P x 1 path delays   [s]
%   nu      : P x 1 path Dopplers [Hz]
%
% Output: J is 4P x 4P.

    [N, M] = size(X);
    P      = numel(alpha);
    T      = 1/deltaf;

    ns = (0:N-1).';
    ms = (0:M-1);

    [NP, NN] = ndgrid(0:N-1, 0:N-1);
    diff_n   = NP - NN;

    G = zeros(N*M, 4*P);

    for p = 1:P
        ap = alpha(p);  tp = tau(p);  np = nu(p);

        % Per-path ICI kernels and their nu-derivatives.
        I_mat  = get_ici_coeffs           (diff_n, np, T);   % N x N
        Ip_mat = get_ici_coeffs_derivative(diff_n, np, T);   % N x N

        % Pre-multiplied data: A_p[n',m] = X[n',m] * e^{-j2pi n' Df tau_p}
        %                                          * e^{+j2pi m Ts nu_p}
        Ap = X .* exp(-1j*2*pi*deltaf*ns*tp) .* exp(1j*2*pi*Ts*ms*np);

        % Doppler field uses I' := (1/j2pi) dI/dnu (ofdm_fim.pdf eq. 8-9);
        % get_ici_coeffs_derivative returns dI/dnu = j2pi*I', hence the 1/(j2pi).
        Phi_p   = I_mat.'  * Ap;
        Psi_p   = I_mat.'  * ((ns*deltaf) .* Ap);
        Omega_p = (ms*Ts) .* Phi_p + (Ip_mat.' * Ap)/(1j*2*pi);

        col = 4*(p-1);
        G(:, col+1) =          Phi_p(:);
        G(:, col+2) = 1j     * Phi_p(:);
        G(:, col+3) = -1j*2*pi*ap * Psi_p(:);
        G(:, col+4) =  1j*2*pi*ap * Omega_p(:);
    end

    J = (2/sigma2) * real(G' * G);
end
