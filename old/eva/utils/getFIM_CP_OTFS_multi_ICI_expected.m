function J = getFIM_CP_OTFS_multi_ICI_expected(N, M, deltaf, Ts, sigma2, alpha, tau, nu)
% Multi-target EXPECTED CP-OTFS Fisher information matrix, with OFDM-side ICI,
% i.i.d. unit-energy QPSK (cp_otfs_fim.tex Sec. 3.2.2, "Expected FIM").
% Parameter stacking eta = [eta_1;...;eta_P], eta_p = [aR_p,aI_p,tau_p,nu_p].
%
% By (I)SFFT unitarity the cross-target kernel-summed moments
% S_a^{p,q}, T_a^{p,q}, Tcheck_a^{p,q}, U^{p,q} coincide with those of the OFDM
% ICI multi FIM (QPSK), so the expected CP-OTFS FIM equals the expected OFDM
% ICI multi FIM block-by-block.  Delegates to getFIM_OFDM_multi_ICI_expected.

    J = getFIM_OFDM_multi_ICI_expected(N, M, deltaf, Ts, sigma2, alpha, tau, nu);
end
