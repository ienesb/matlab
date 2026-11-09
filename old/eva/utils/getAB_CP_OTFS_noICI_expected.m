function [A, B] = getAB_CP_OTFS_noICI_expected(N, M, deltaf, Ts, sigma2, alpha, tau, nu, p)
% Single-target EXPECTED MCRB matrices A, B for CP-OTFS, no OFDM-side ICI,
% i.i.d. unit-energy QPSK with delay-Doppler BSC rate p (cp_otfs_mcrb.tex
% Sec. 1.1.1 / 2.1.1).  eta = [alphaR, alphaI, tau, nu].
%
%   E{A} = -E{Jtilde} + E{A_Delta}                       (exact; = OFDM expected A)
%   E{B} = E{Jtilde} + E{b}E{b}^T + Cov(b)_CP-OTFS        (exact)
%
% By (I)SFFT unitarity + i.i.d. QPSK, E{A}, E{Jtilde} and the mean bias E{b} all
% coincide with the OFDM no-ICI MCRB (cp_otfs_mcrb.tex Sec. 1-2), so they are
% reused from the OFDM expected functions.  The bias covariance Cov(b) is the
% CP-OTFS-specific one of cp_otfs_mcrb.tex Sec. 4 (getMCRB_covb_CP_OTFS): the BSC
% acts in the delay-Doppler domain and the ISFFT Dirichlet kernels mix the whole
% grid, so Cov(b) is a single NM-dim bilinear form (trace/Frobenius parts match
% OFDM by unitary invariance, the diagonal QPSK 4th-cumulant parts differ).

    [A, Bofdm] = getAB_OFDM_noICI_expected(N, M, deltaf, Ts, sigma2, alpha, tau, nu, p);
    Cb_ofdm   = getMCRB_covb_OFDM   (N, M, deltaf, Ts, sigma2, alpha, tau, nu, p, false);
    Cb_cpotfs = getMCRB_covb_CP_OTFS(N, M, deltaf, Ts, sigma2, alpha, tau, nu, p, false);
    B = Bofdm - Cb_ofdm + Cb_cpotfs;     % = E{Jtilde} + E{b}E{b}^T + Cov(b)_CP-OTFS
end
