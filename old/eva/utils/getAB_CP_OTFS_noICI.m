function [A, B] = getAB_CP_OTFS_noICI(x, xhat, deltaf, Ts, sigma2, alpha, tau, nu, eta_eval)
% Single-target realization MCRB component matrices A, B for CP-OTFS, no
% OFDM-side ICI (cp_otfs_mcrb.tex Sec. 1.1.1 / 2.1).  Parameter order
% eta = [alphaR, alphaI, tau, nu].  Returns 4x4 A, B.
%
% By (I)SFFT unitarity all derivative/Hessian fields and the mismatch transform
% by the same unitary SFFT, so per realization the CP-OTFS MCRB equals the OFDM
% MCRB at X = ISFFT(x), Xhat = ISFFT(xhat).  Delegates accordingly.
%
% Inputs:
%   x, xhat : N x M delay-Doppler true / receiver-side symbol grids
%   deltaf, Ts, sigma2, alpha (scalar), tau, nu.

% eta_eval: (optional) 4 x 1 EVALUATION point (e.g. pseudo-true); default = ground truth.

    if nargin < 9, eta_eval = []; end
    X    = isfft_u(x);
    Xhat = isfft_u(xhat);
    [A, B] = getAB_OFDM_noICI(X, Xhat, deltaf, Ts, sigma2, alpha, tau, nu, eta_eval);
end
