function J = getFIM_OFDM_ICI(X, deltaf, Ts, sigma2, alpha, tau, nu)
% Single-target OFDM Fisher information matrix with the full
% Doppler-induced ICI kernel (Section 3.1.2 of ofdm_fim.tex).
% Parameter vector ordering: eta = [alphaR, alphaI, tau, nu].
%
% Inputs:
%   X       : N x M transmitted symbol grid (subcarrier x OFDM symbol)
%   deltaf  : subcarrier spacing [Hz]
%   Ts      : OFDM symbol duration including CP [s]
%   sigma2  : noise variance
%   alpha   : complex path gain (scalar)
%   tau     : path delay [s]
%   nu      : path Doppler [Hz]
%
% Output:
%   J       : 4 x 4 Fisher information matrix.
%
% Uses the dimensionless ICI kernel from crb/EVA/get_ici_coeffs:
%   I(n)  = (e^{j2pi nu T} - 1) / (j2pi (n + nu T))
%   I'(n) = d I(n)/d nu
% with I(0) = 1 when nu = 0 so the no-ICI limit is the Kronecker delta.

    [N, M] = size(X);
    T  = 1/deltaf;

    ns = (0:N-1).';         % n'  (subcarrier index, column)
    ms = (0:M-1);           % m   (OFDM-symbol index, row)

    % Build the N x N kernel and derivative matrices once.  Convention:
    %   I_mat(np+1, n+1) = I(np - n)   (np varies along dim 1, n along dim 2)
    [NP, NN] = ndgrid(0:N-1, 0:N-1);
    diff_n   = NP - NN;
    I_mat  = get_ici_coeffs(diff_n, nu, T);             % N x N
    Ip_mat = get_ici_coeffs_derivative(diff_n, nu, T);  % N x N

    % A[n', m] = X[n',m] * e^{-j2pi n' Df tau} * e^{+j2pi m Ts nu}
    A = X .* exp(-1j*2*pi*deltaf*ns*tau) .* exp(1j*2*pi*Ts*ms*nu);

    % Phi[n, m]   = sum_{n'} A[n',m] I(n'-n)        = (I_mat.' * A)
    % Psi[n, m]   = sum_{n'} A[n',m] n'Df I(n'-n)   = (I_mat.' * (ns*df.*A))
    % Omega[n, m] = sum_{n'} A[n',m] (mTs I + I')   = mTs*Phi + (Ip_mat.' * A)/(j2pi)
    % Per ofdm_fim.pdf (eq. 8-9) the Doppler field uses I' := (1/j2pi) dI/dnu,
    % while get_ici_coeffs_derivative returns the full dI/dnu = j2pi*I'; hence
    % the 1/(j2pi) scaling on the I' contribution below.
    Phi   = I_mat.'  * A;
    Psi   = I_mat.'  * ((ns*deltaf) .* A);
    Omega = (ms*Ts) .* Phi + (Ip_mat.' * A)/(1j*2*pi);

    g_aR  = Phi(:);
    g_aI  = 1j * g_aR;
    g_tau = -1j*2*pi*alpha * Psi(:);
    g_nu  =  1j*2*pi*alpha * Omega(:);

    G = [g_aR, g_aI, g_tau, g_nu];
    J = (2/sigma2) * real(G' * G);
end
