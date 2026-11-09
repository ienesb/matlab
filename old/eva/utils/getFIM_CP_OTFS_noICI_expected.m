function J = getFIM_CP_OTFS_noICI_expected(N, M, deltaf, Ts, sigma2, alpha, tau, nu) %#ok<INUSD>
% Single-target EXPECTED CP-OTFS Fisher information matrix, no OFDM-side ICI,
% i.i.d. unit-energy QPSK (cp_otfs_fim.tex Sec. 3.1.1, "Expected FIM").
% Parameter order eta = [alphaR, alphaI, tau, nu]; tau, nu do not enter.
%
% By (I)SFFT unitarity and i.i.d. unit-energy data, the expected CP-OTFS FIM
% equals the expected OFDM no-ICI FIM entry-by-entry (proven via the Bartlett
% identity in the document): the lag-summed Dirichlet moments A_{a,b}, B_{a,b}
% collapse to the elementary power sums P_a, Q_b.  Hence this delegates to
% getFIM_OFDM_noICI_expected.  NOTE: per realization the CP-OTFS FIM is random
% (it fluctuates even without ICI), but its MEAN is this deterministic matrix.

    J = getFIM_OFDM_noICI_expected(N, M, deltaf, Ts, sigma2, alpha);
end
