function [A, B] = getAB_CP_OTFS_multi_ICI(x, xhat, deltaf, Ts, sigma2, alpha, tau, nu, eta_eval)
% Multi-target realization MCRB component matrices A, B for CP-OTFS, with
% OFDM-side ICI (cp_otfs_mcrb.tex Sec. 1.2.2 / 2.2.2).  Stacking
% eta = [eta_1;...;eta_P], eta_p = [aR_p,aI_p,tau_p,nu_p].  Returns 4P x 4P.
%
% By (I)SFFT unitarity the per-realization CP-OTFS MCRB equals the OFDM ICI
% multi MCRB at X = ISFFT(x), Xhat = ISFFT(xhat).  Delegates accordingly.
% eta_eval: (optional) 4P x 1 EVALUATION point (e.g. pseudo-true); default = ground truth.

    if nargin < 9, eta_eval = []; end
    X    = isfft_u(x);
    Xhat = isfft_u(xhat);
    [A, B] = getAB_OFDM_multi_ICI(X, Xhat, deltaf, Ts, sigma2, alpha, tau, nu, eta_eval);
end
