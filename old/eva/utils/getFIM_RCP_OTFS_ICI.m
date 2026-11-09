function J = getFIM_RCP_OTFS_ICI(x, deltaf, Ts, sigma2, alpha, tau, nu)
% Single-target RCP-OTFS Fisher information matrix WITH OFDM-side ICI
% (rcp_otfs_fim.tex, P=1).  Reduced-CP OTFS prepends one cyclic prefix per
% frame, adding the intra-frame Doppler phase and the boundary-wrap phase to
% the delay-Doppler mean (T_S -> T).  Parameter ordering: eta = [aR,aI,tau,nu].
%
% Wrapper over the multi-target builder with a single path.

    J = getFIM_RCP_OTFS_multi(x, deltaf, Ts, sigma2, alpha, tau, nu, true);
end
