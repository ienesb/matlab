function J = getFIM_RCP_OFDM_multi_ICI(X, deltaf, Ts, sigma2, alpha, tau, nu)
% Multi-target RCP-OFDM Fisher information matrix WITH the full Doppler
% ICI kernels (rcp_ofdm_fim.tex).  Reduced-CP OFDM: a single cyclic prefix
% is prepended to the whole frame instead of one per OFDM symbol, so delay
% additionally causes frame-cyclic ISI from the previous symbol.
%
% Thin wrapper over getFIM_RCP_OFDM_multi with the ICI kernel enabled.
% Parameter ordering: eta = [a1R,a1I,tau1,nu1, ...].  J is 4P x 4P.

    J = getFIM_RCP_OFDM_multi(X, deltaf, Ts, sigma2, alpha, tau, nu, true);
end
