function J = getFIM_OTFS(x, deltaf, Ts, sigma2, alpha, tau, nu)
% Single-target FIM for OTFS (Section 3.1 of otfs_fim.tex).
% Parameter order: eta = [alphaR, alphaI, tau, nu].
%
% Inputs:
%   x       : N x M delay-Doppler symbol grid (delay rows, Doppler cols)
%   deltaf  : subcarrier spacing [Hz]
%   Ts      : OFDM symbol duration (with CP) [s]
%   sigma2  : noise variance
%   alpha   : complex path gain (scalar)
%   tau     : path delay [s]
%   nu      : path Doppler [Hz]
%
% Output:
%   J       : 4 x 4 Fisher information matrix.

    [N, M] = size(x);

    % ---- Delay and Doppler kernels A(d), A'(d), B(d), B'(d) ----
    % A(d)  = sum_{n=0..N-1} exp(j 2pi n (d/N - df*tau))
    % A'(d) = sum_{n=0..N-1} n * exp(j 2pi n (d/N - df*tau))
    ns = (0:N-1).';                          % N x 1
    ds_A = (-(N-1):(N-1)).';                 % (2N-1) x 1
    Earg_A = exp(1j*2*pi*ns*(ds_A.'/N - deltaf*tau));   % N x (2N-1)
    A_d  = sum(Earg_A, 1).';                 % (2N-1) x 1
    Ap_d = sum(ns .* Earg_A, 1).';

    % B(d)  = sum_{m=0..M-1} exp(j 2pi m (d/M + Ts*nu))
    % B'(d) = sum_{m=0..M-1} m * exp(j 2pi m (d/M + Ts*nu))
    ms = (0:M-1).';
    ds_B = (-(M-1):(M-1)).';
    Earg_B = exp(1j*2*pi*ms*(ds_B.'/M + Ts*nu));         % M x (2M-1)
    B_d  = sum(Earg_B, 1).';
    Bp_d = sum(ms .* Earg_B, 1).';

    % ---- Toeplitz matrices A_mat[l',l] = A(l'-l), B_mat[k,k'] = B(k-k') ----
    % toeplitz(col, row): first column carries d = 0, 1, ..., N-1; first row d = 0, -1, ..., -(N-1).
    A_mat  = toeplitz(A_d (N:end), A_d (N:-1:1));      % N x N
    Ap_mat = toeplitz(Ap_d(N:end), Ap_d(N:-1:1));
    B_mat  = toeplitz(B_d (M:end), B_d (M:-1:1));      % M x M (k row, k' col)
    Bp_mat = toeplitz(Bp_d(M:end), Bp_d(M:-1:1));

    % ---- Φ, Ψ, Ω as N x M matrices ----
    Phi   = (A_mat  * x * B_mat ) / (N*M);
    Psi   = (Ap_mat * x * B_mat ) / (N*M);
    Omega = (A_mat  * x * Bp_mat) / (N*M);

    % ---- Stack gradients and form FIM ----
    g_aR  = Phi(:);
    g_aI  = 1j * g_aR;
    g_tau = -1j*2*pi*deltaf * alpha * Psi(:);
    g_nu  =  1j*2*pi*Ts     * alpha * Omega(:);

    G = [g_aR, g_aI, g_tau, g_nu];
    J = (2/sigma2) * real(G' * G);
end
