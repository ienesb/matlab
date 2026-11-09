function [A, B] = getAB_CP_OTFS_ICI(x, xhat, deltaf, Ts, sigma2, alpha, tau, nu, eta_eval)
% Single-target realization MCRB component matrices A, B for CP-OTFS, with
% OFDM-side ICI (cp_otfs_mcrb.tex Sec. 1.1.2 / 2.1.2).  eta=[aR,aI,tau,nu].
%
% By (I)SFFT unitarity the per-realization CP-OTFS MCRB equals the OFDM ICI
% MCRB at X = ISFFT(x), Xhat = ISFFT(xhat).  Delegates accordingly.
%
% Inputs:
%   x, xhat : N x M delay-Doppler true / receiver-side symbol grids
%   deltaf, Ts, sigma2, alpha (scalar), tau, nu.
%   eta_eval: (optional) 4 x 1 EVALUATION point (e.g. pseudo-true); default = ground truth.

    if nargin < 9, eta_eval = []; end
    X    = isfft_u(x);
    Xhat = isfft_u(xhat);
    [A, B] = getAB_OFDM_ICI(X, Xhat, deltaf, Ts, sigma2, alpha, tau, nu, eta_eval);
end
