function J = getFIM_OFDM_ICI(X, deltaf, Ts, sigma2, alpha, tau, nu)
% Single-target FIM for OFDM, ICI channel (Section 3.1.2 of ofdm_fim.tex).
% Parameter order: eta = [alphaR, alphaI, tau, nu].
%
% Uses the existing dimensionless ICI kernel:
%   I(n)  = (e^{j2pi nu T} - 1) / (j2pi (n + nu T))
%   I'(n) = d I(n) / d nu
% which satisfies I(0)=1 when nu=0, recovering the no-ICI Kronecker-delta limit.
%
% Inputs / Output: same shape as getFIM_OFDM_noICI.

    [N, M] = size(X);
    T = 1 / deltaf;

    ns = (0:N-1).';      % subcarrier index n' (column vector)
    ms = (0:M-1);        % OFDM symbol index m (row vector)

    % Frequency phase X * e^{-j2pi n' Df tau} -- depends on n' (rows), not on n.
    freq_phase = X .* exp(-1j*2*pi*deltaf*ns*tau);                  % N x M
    % Time phase e^{+j2pi m Ts nu} -- depends on m only.
    time_phase = exp(1j*2*pi*Ts*ms*nu);                             % 1 x M

    % Pre-multiply: A[n', m] = X[n',m] * e^{-j2pi n' Df tau} * e^{+j2pi m Ts nu}
    A = freq_phase .* time_phase;                                    % N x M

    Phi   = zeros(N, M);   % sum_{n'} A[n',m] I(n'-n)
    Psi   = zeros(N, M);   % sum_{n'} A[n',m] n' Df I(n'-n)
    Omega = zeros(N, M);   % sum_{n'} A[n',m] ( m Ts I + I' )

    for n = 0:N-1
        I  = get_ici_coeffs(ns - n, nu, T);              % N x 1
        dI = get_ici_coeffs_derivative(ns - n, nu, T);   % N x 1

        Phi(n+1, :)   = sum(A .* I, 1);
        Psi(n+1, :)   = sum(A .* (ns*deltaf) .* I, 1);
        Omega(n+1, :) = sum(A .* (ms*Ts .* I + dI), 1);
    end

    g_aR  = Phi(:);
    g_aI  = 1j * g_aR;
    g_tau = -1j*2*pi*alpha * Psi(:);
    g_nu  =  1j*2*pi*alpha * Omega(:);

    G = [g_aR, g_aI, g_tau, g_nu];                      % NM x 4

    J = (2/sigma2) * real(G' * G);
end
