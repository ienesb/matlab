function [A, B] = getAB_OFDM_ICI(X, Xhat, deltaf, Ts, sigma2, alpha, tau, nu, eta_eval)
% Single-target realization MCRB component matrices A and B for OFDM,
% ICI channel (ofdm_mcrb.tex Sec. 1.1.2 / 2.1.2).
% Parameter order: eta = [alphaR, alphaI, tau, nu].  Returns 4x4 A and B.
%
%   B = Jtilde + b b^T ,    A = -Jtilde + A_Delta ,
% with Jtilde the ICI FIM evaluated at Xhat and the bias/Hessian corrections
% driven by the mismatch E = X - Xhat.  Uses the unit-peak ICI kernel and the
% I' = (1/j2pi) dI/dnu, I'' = (1/j2pi) dI'/dnu convention.  Xhat = X gives the
% CRB (B = J, A = -J).
%
% Inputs:
%   X, Xhat : N x M true and receiver-side symbol grids
%   deltaf  : subcarrier spacing [Hz]
%   Ts      : OFDM symbol duration incl. CP [s]
%   sigma2  : noise variance
%   alpha   : complex path gain (scalar)        -- TRUE parameter eta_bar
%   tau, nu : path delay [s], Doppler [Hz]       -- TRUE parameter eta_bar
%   eta_eval: (optional) 4 x 1 EVALUATION point [aR;aI;tau;nu] at which A,B are
%             formed (e.g. the pseudo-true eta_0).  Default = ground truth.

    if nargin < 9 || isempty(eta_eval)
        aE = alpha;  tE = tau;  nE = nu;
    else
        [aE, tE, nE] = unpackEta(eta_eval);
    end
    [A, B] = getAB_OFDM_core(X, Xhat, deltaf, Ts, sigma2, aE, tE, nE, true, alpha, tau, nu);
end
