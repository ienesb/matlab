function [eta0, alpha0, tau0, nu0, resnorm, exitflag] = ...
        pseudoTrue_CP_OTFS(x, xhat, deltaf, Ts, alpha_bar, tau_bar, nu_bar, use_ici, eta_init, use_varpro)
% Pseudo-true parameter for the misspecified CP-OTFS ISAC model (pseudo_true.tex).
% Solves the delay-Doppler least-squares projection
%
%   eta0 = argmin_eta || mu_DD(x, eta_bar) - mu_tilde_DD(xhat, eta) ||,
%
% where mu_DD / mu_tilde_DD are the CP-OTFS noise-free means built from the
% transmitted / detected DELAY-DOPPLER symbols x, xhat.
%
% By (I)SFFT unitarity the DD-domain norm is preserved under the ISFFT, so the
% projection is identical to the OFDM (time-frequency) projection evaluated at
% the ISFFT-mapped grids:
%   || mu_DD(x,.) - mu_tilde_DD(xhat,.) || = || mu_TF(ISFFT x,.) - mu_tilde_TF(ISFFT xhat,.) ||.
% The channel parameter eta is the same object in both domains, so the
% pseudo-true parameter is obtained by delegating to pseudoTrue_OFDM at
% X = ISFFT(x), Xhat = ISFFT(xhat) -- the same equivalence used by the
% getAB_CP_OTFS_* and getFIM_CP_OTFS_* families.  The ISFFT is applied once
% here (not inside the optimizer), so this costs the same as the OFDM solve.
%
% Handles all four cases via (numel(alpha_bar), use_ici): single/multi target,
% no-ICI/ICI.  Parameter stacking eta = [eta_1;...;eta_P], eta_p=[aR;aI;tau;nu].
%
% Inputs:
%   x, xhat   : N x M transmitted / detected delay-Doppler symbol grids
%   deltaf,Ts : subcarrier spacing [Hz], OFDM symbol duration [s]
%   alpha_bar,tau_bar,nu_bar : P x 1 TRUE parameters eta_bar
%   use_ici   : true -> OFDM-side ICI kernel, false -> no-ICI
%   eta_init  : (optional) 4P x 1 initial guess; default = the true parameter
%   use_varpro: (optional) true -> variable-projection (gain-concentrated) solve
%               over the 2P delay/Doppler unknowns; false (default) -> full 4P.
%               Passed through to pseudoTrue_OFDM; by (I)SFFT unitarity the
%               projection is identical in the TF and DD domains (cp_otfs_mcrb.tex
%               Sec. "Pseudo-true parameter", equivalence note).
%
% Outputs:
%   eta0            : 4P x 1 pseudo-true parameter [aR;aI;tau;nu] per target
%   alpha0,tau0,nu0 : P x 1 unpacked pseudo-true gains/delays/Dopplers
%   resnorm,exitflag: lsqnonlin diagnostics

    X    = isfft_u(x);
    Xhat = isfft_u(xhat);
    if nargin < 9,  eta_init   = []; end
    if nargin < 10, use_varpro = false; end
    [eta0, alpha0, tau0, nu0, resnorm, exitflag] = ...
        pseudoTrue_OFDM(X, Xhat, deltaf, Ts, alpha_bar, tau_bar, nu_bar, use_ici, eta_init, use_varpro);
end
