function [A, B] = getAB_OFDM_noICI_expected(N, M, deltaf, Ts, sigma2, alpha, tau, nu, p)
% Single-target EXPECTED MCRB component matrices A and B for OFDM, no-ICI
% channel, i.i.d. unit-energy QPSK with BSC bit-error rate p (ofdm_mcrb.tex
% Sec. 1.1.1 / 2.1.1, "Expected ..." paragraphs).  Parameter order
% eta = [alphaR, alphaI, tau, nu].  tau, nu do not enter (no-ICI).
%
%   E{A} = -E{Jtilde} + E{A_Delta},     E{B} = E{Jtilde} + E{b} E{b}^T.
% Here E{Jtilde} = J (the QPSK FIM) and, using E{Gamma^(a,b)} = -2p df^a Ts^b P_a Q_b,
%   E{A_Delta}_kl = 2p [J]_kl  (kinematic & gain-kinematic),  0 (gain-gain),
% so  E{A}_kl = -(1-2p)[J]_kl  except the gain-gain block where E{A}_kl = -[J]_kl.
% The bias mean is E{b} = [-4p aR NM/s2, -4p aI NM/s2, 0, 0].
%
% NOTE (documented gap): the exact E{B} also contains the bias covariance
% Cov(b) = E{b b^T} - E{b}E{b}^T, which ofdm_mcrb.tex does NOT provide.  This
% function returns the doc's mean-bias form E{Jtilde} + E{b}E{b}^T; in
% particular the (tau,nu) bias-variance terms (whose mean is zero) are absent.
%
% Inputs: grid (N,M), deltaf, Ts, sigma2, complex gain alpha, tau, nu, BSC p.
% Output: 4x4 A, B.

    J = getFIM_OFDM_noICI_expected(N, M, deltaf, Ts, sigma2, alpha);

    % E{A}: scale everything by (1-2p) then restore the gain-gain block.
    A = -(1-2*p) * J;
    A(1:2,1:2) = -J(1:2,1:2);

    % E{b} and mean-bias E{B}.
    aR = real(alpha);  aI = imag(alpha);
    b = [ -4*p*aR*N*M/sigma2;
          -4*p*aI*N*M/sigma2;
           0;
           0 ];
    % Exact E{B} = E{Jtilde} + E{b}E{b}^T + Cov(b).
    Cb = getMCRB_covb_OFDM(N, M, deltaf, Ts, sigma2, alpha, tau, nu, p, false);
    B = J + b*b.' + Cb;
end
