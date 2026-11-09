function J = getFIM_OFDM_noICI_expected(N, M, deltaf, Ts, sigma2, alpha)
% Single-target EXPECTED OFDM Fisher information matrix, no-ICI channel.
% Closed form for i.i.d. unit-energy QPSK symbols (|X[n,m]|^2 = 1), i.e. the
% data-deterministic FIM of Section 3.1.1 / eqs. (98)-(107) of ofdm_fim.tex.
% Parameter ordering: eta = [alphaR, alphaI, tau, nu].
%
% The no-ICI FIM does NOT depend on the true (tau, nu): the tau-phase factors
% cancel pairwise and the no-ICI kernel is nu-independent.  It is therefore a
% function of the grid (N, M, deltaf, Ts), the noise variance, and the gain.
%
% Inputs:
%   N, M    : grid size (subcarriers x OFDM symbols)
%   deltaf  : subcarrier spacing [Hz]
%   Ts      : OFDM symbol duration including CP [s]
%   sigma2  : noise variance
%   alpha   : complex path gain (scalar)
%
% Output:
%   J       : 4 x 4 Fisher information matrix.

    aR = real(alpha);  aI = imag(alpha);  a2 = abs(alpha)^2;

    % Elementary power sums  P_a = sum_{n=0}^{N-1} n^a ,  Q_b = sum_{m=0}^{M-1} m^b.
    P0 = N;   P1 = N*(N-1)/2;   P2 = N*(N-1)*(2*N-1)/6;
    Q0 = M;   Q1 = M*(M-1)/2;   Q2 = M*(M-1)*(2*M-1)/6;   %#ok<NASGU>

    J = zeros(4,4);

    % alpha-alpha block (eqs. 98-99)
    J(1,1) = (2/sigma2) * P0 * Q0;        % = 2 N M / sigma2
    J(2,2) = J(1,1);
    J(1,2) = 0;  J(2,1) = 0;

    % tau-tau, nu-nu, tau-nu (eqs. 100-102)
    J(3,3) =  (8*pi^2*a2/sigma2) * deltaf^2 * (P2 * Q0);
    J(4,4) =  (8*pi^2*a2/sigma2) * Ts^2     * (P0 * Q2);
    J(3,4) = -(8*pi^2*a2/sigma2) * deltaf*Ts* (P1 * Q1);
    J(4,3) =  J(3,4);

    % alpha-tau, alpha-nu cross terms (eqs. 103-106)
    J(1,3) =  (4*pi*aI/sigma2) * deltaf * (P1 * Q0);
    J(1,4) = -(4*pi*aI/sigma2) * Ts     * (P0 * Q1);
    J(2,3) = -(4*pi*aR/sigma2) * deltaf * (P1 * Q0);
    J(2,4) =  (4*pi*aR/sigma2) * Ts     * (P0 * Q1);

    % Symmetric completion.
    J(3,1) = J(1,3);  J(4,1) = J(1,4);
    J(3,2) = J(2,3);  J(4,2) = J(2,4);
end
