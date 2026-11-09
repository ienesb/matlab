function [A, B] = getAB_CP_OTFS_ICI_expected(N, M, deltaf, Ts, sigma2, alpha, tau, nu, p)
% Single-target EXPECTED MCRB matrices A, B for CP-OTFS, with OFDM-side ICI,
% i.i.d. unit-energy QPSK with delay-Doppler BSC rate p (cp_otfs_mcrb.tex
% Sec. 1.1.2 / 2.1.2).  eta = [alphaR, alphaI, tau, nu].
%
%   E{A} = -E{Jtilde} + E{A_Delta}                  (exact; = OFDM ICI expected A)
%   E{B} = E{Jtilde} + E{b}E{b}^T + Cov(b)_CP-OTFS   (exact)
%
% E{A}, E{Jtilde}, E{b} reuse the OFDM ICI expected functions (unitarity).  The
% bias covariance is the CP-OTFS-specific one of cp_otfs_mcrb.tex Sec. 4
% (getMCRB_covb_CP_OTFS), not the OFDM Cov(b): the DD-domain BSC + ISFFT mixing
% makes it a single NM-dim bilinear form (with ICI even the gain block is
% nonzero, since K^H K != N^2 I).

    [A, Bofdm] = getAB_OFDM_ICI_expected(N, M, deltaf, Ts, sigma2, alpha, tau, nu, p);
    Cb_ofdm   = getMCRB_covb_OFDM   (N, M, deltaf, Ts, sigma2, alpha, tau, nu, p, true);
    Cb_cpotfs = getMCRB_covb_CP_OTFS(N, M, deltaf, Ts, sigma2, alpha, tau, nu, p, true);
    B = Bofdm - Cb_ofdm + Cb_cpotfs;     % = E{Jtilde} + E{b}E{b}^T + Cov(b)_CP-OTFS
end
