function J = getFIM_RCP_OTFS_multi_ICI(x, deltaf, Ts, sigma2, alpha, tau, nu)
% Multi-target RCP-OTFS Fisher information matrix WITH OFDM-side ICI
% (rcp_otfs_fim.tex).  Reduced-CP OTFS: a single cyclic prefix is prepended
% to the whole frame instead of one per OFDM symbol, which adds the
% intra-frame Doppler phase and the boundary-wrap phase to the delay-Doppler
% mean and replaces the per-symbol period T_S by the useful duration T.
%
% Thin wrapper over getFIM_RCP_OTFS_multi with the ICI kernel enabled.
% Parameter ordering: eta = [a1R,a1I,tau1,nu1, ...].  J is 4P x 4P.

    J = getFIM_RCP_OTFS_multi(x, deltaf, Ts, sigma2, alpha, tau, nu, true);
end
