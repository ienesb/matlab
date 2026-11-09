function J = getFIM_CP_OTFS_ICI_expected(N, M, deltaf, Ts, sigma2, alpha, tau, nu) %#ok<INUSL>
% Single-target EXPECTED CP-OTFS Fisher information matrix, with OFDM-side ICI,
% i.i.d. unit-energy QPSK (cp_otfs_fim.tex Sec. 3.1.2, "Expected FIM").
% Parameter order eta = [alphaR, alphaI, tau, nu]; depends on nu, not tau.
%
% By (I)SFFT unitarity, the two-index ICI delay kernel K(l',l) collapses under
% sum_{l',l} to N^2 times the OFDM kernel sums S_a, T_a, U, so the expected
% CP-OTFS-with-ICI FIM equals the expected OFDM ICI FIM entry-by-entry.
% Delegates to getFIM_OFDM_ICI_expected.  (Per realization the CP-OTFS FIM is
% random; this is its mean.)

    J = getFIM_OFDM_ICI_expected(N, M, deltaf, Ts, sigma2, alpha, nu);
end
