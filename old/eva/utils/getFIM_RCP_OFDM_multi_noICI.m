function J = getFIM_RCP_OFDM_multi_noICI(X, deltaf, Ts, sigma2, alpha, tau, nu)
% Multi-target RCP-OFDM Fisher information matrix WITHOUT Doppler-induced
% ICI (rcp_ofdm_fim.tex): kernels frozen at nu = 0 and their
% nu-derivatives dropped; Doppler enters only through the per-window phase
% e^{j2pi nu m T}.  The delay-induced current/previous symbol split is
% retained -- it is structural to the reduced CP, present even at nu = 0.
%
% Thin wrapper over getFIM_RCP_OFDM_multi with the ICI kernel disabled.
% Parameter ordering: eta = [a1R,a1I,tau1,nu1, ...].  J is 4P x 4P.

    J = getFIM_RCP_OFDM_multi(X, deltaf, Ts, sigma2, alpha, tau, nu, false);
end
