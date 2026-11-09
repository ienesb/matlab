function J = getFIM_CP_OTFS_multi_noICI_expected(N, M, deltaf, Ts, sigma2, alpha, tau, nu)
% Multi-target EXPECTED CP-OTFS Fisher information matrix, no OFDM-side ICI,
% i.i.d. unit-energy QPSK (cp_otfs_fim.tex Sec. 3.2.1, "Expected FIM").
% Parameter stacking eta = [eta_1;...;eta_P], eta_p = [aR_p,aI_p,tau_p,nu_p].
%
% By (I)SFFT unitarity the cross-target Dirichlet ambiguity moments equal the
% OFDM no-ICI cross-target moments Lambda_{p,q}^{(a,b)} (QPSK), so the expected
% CP-OTFS FIM equals the expected OFDM no-ICI multi FIM block-by-block.
% Delegates to getFIM_OFDM_multi_noICI_expected.

    J = getFIM_OFDM_multi_noICI_expected(N, M, deltaf, Ts, sigma2, alpha, tau, nu);
end
