function [A, B] = getAB_CP_OTFS_multi_ICI_expected(N, M, deltaf, Ts, sigma2, alpha, tau, nu, p)
% Multi-target EXPECTED MCRB matrices A, B for CP-OTFS, with OFDM-side ICI,
% i.i.d. unit-energy QPSK with delay-Doppler BSC rate p (cp_otfs_mcrb.tex
% Sec. 1.2.2 / 2.2.2).  Stacking eta=[eta_1;...;eta_P].
%
%   E{A} = -E{Jtilde} + E{A_Delta}                  (exact; = OFDM ICI multi A)
%   E{B} = E{Jtilde} + E{b}E{b}^T + Cov(b)_CP-OTFS   (exact)
%
% E{A}, E{Jtilde}, E{b} reuse the OFDM expected functions (unitarity); the bias
% covariance is the CP-OTFS-specific one of cp_otfs_mcrb.tex Sec. 4
% (getMCRB_covb_CP_OTFS, cross-target block form), not OFDM's.

    [A, Bofdm] = getAB_OFDM_multi_ICI_expected(N, M, deltaf, Ts, sigma2, alpha, tau, nu, p);
    Cb_ofdm   = getMCRB_covb_OFDM   (N, M, deltaf, Ts, sigma2, alpha, tau, nu, p, true);
    Cb_cpotfs = getMCRB_covb_CP_OTFS(N, M, deltaf, Ts, sigma2, alpha, tau, nu, p, true);
    B = Bofdm - Cb_ofdm + Cb_cpotfs;
end
