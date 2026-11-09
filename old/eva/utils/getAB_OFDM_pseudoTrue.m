function [A, B, eta0, alpha0, tau0, nu0] = ...
        getAB_OFDM_pseudoTrue(X, Xhat, deltaf, Ts, sigma2, alpha, tau, nu, use_ici)
% MCRB component matrices A and B for OFDM evaluated at the PSEUDO-TRUE
% parameter eta_0 (pseudo_true.tex), rather than at the ground-truth eta_bar.
%
% Standard MCRB theory evaluates A and B at the value to which the mismatched
% ML estimator converges -- the pseudo-true parameter, the least-squares
% projection of the true mean onto the assumed-mean manifold:
%   eta_0 = argmin_eta || mu(X, eta_bar) - mu_tilde(Xhat, eta) ||.
% This wrapper (i) solves for eta_0 with pseudoTrue_OFDM, then (ii) builds A, B
% with getAB_OFDM_core at the EVALUATION point eta_0 while keeping the TRUE mean
% in the mismatch Delta = mu(X, eta_bar) - mu_tilde(Xhat, eta_0).
%
% Handles all four cases via (numel(alpha), use_ici): single/multi target,
% no-ICI/ICI.  Parameter stacking eta = [eta_1;...;eta_P], eta_p=[aR;aI;tau;nu].
%
% Inputs:
%   X, Xhat   : N x M transmitted / detected symbol grids
%   deltaf,Ts : subcarrier spacing [Hz], OFDM symbol duration [s]
%   sigma2    : noise variance
%   alpha,tau,nu : P x 1 TRUE parameters eta_bar (gains/delays/Dopplers)
%   use_ici   : true -> ICI kernel, false -> no-ICI (I = delta)
%
% Outputs:
%   A, B            : 4P x 4P MCRB component matrices at eta_0
%   eta0            : 4P x 1 pseudo-true parameter [aR;aI;tau;nu] per target
%   alpha0,tau0,nu0 : P x 1 unpacked pseudo-true gains/delays/Dopplers
%
% With Xhat = X (no detection errors) eta_0 = eta_bar, Delta = 0, so A = -J,
% B = J and the MCRB collapses to the CRB.

    % (i) pseudo-true parameter = LS projection (pseudo_true.tex last equation).
    [eta0, alpha0, tau0, nu0] = ...
        pseudoTrue_OFDM(X, Xhat, deltaf, Ts, alpha, tau, nu, use_ici);

    % (ii) A, B at the evaluation point eta_0 with the true mean from eta_bar.
    [A, B] = getAB_OFDM_core(X, Xhat, deltaf, Ts, sigma2, ...
                             alpha0, tau0, nu0, use_ici, alpha, tau, nu);
end
