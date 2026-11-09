function [A, B] = getAB_OFDM_multi_noICI(X, Xhat, deltaf, Ts, sigma2, alpha, tau, nu, eta_eval)
% Multi-target realization MCRB component matrices A and B for OFDM,
% no-ICI channel (ofdm_mcrb.tex Sec. 1.2.1 / 2.2.1).
% Parameter stacking: eta = [eta_1; ...; eta_P], eta_p = [aR_p, aI_p, tau_p, nu_p].
% Returns 4P x 4P A and B.
%
%   B = Jtilde + b b^T ,    A = -Jtilde + A_Delta ,
% where Jtilde is the multi-target no-ICI FIM at Xhat, A_Delta is block-diagonal
% (cross-target Hessians vanish) while the off-diagonal blocks A_{p,q} = -Jtilde_{p,q}.
%
% Inputs:
%   X, Xhat : N x M true and receiver-side symbol grids
%   deltaf  : subcarrier spacing [Hz]
%   Ts      : OFDM symbol duration incl. CP [s]
%   sigma2  : noise variance
%   alpha   : P x 1 complex path gains
%   tau     : P x 1 path delays   [s]
%   nu      : P x 1 path Dopplers [Hz]

    % (alpha,tau,nu) are the TRUE parameter eta_bar; optional eta_eval (4P x 1
    % stacked [aR;aI;tau;nu]) is the EVALUATION point, default = ground truth.
    if nargin < 9 || isempty(eta_eval)
        aE = alpha;  tE = tau;  nE = nu;
    else
        [aE, tE, nE] = unpackEta(eta_eval);
    end
    [A, B] = getAB_OFDM_core(X, Xhat, deltaf, Ts, sigma2, aE, tE, nE, false, alpha, tau, nu);
end
