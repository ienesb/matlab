function [eta_hat, alpha_hat, tau_hat, nu_hat] = ...
        mmle_CP_OTFS(y, xhat, deltaf, Ts, P, use_ici, Zpad)
% Misspecified maximum-likelihood estimator (MMLE) for the CP-OTFS ISAC model
% (cp_otfs_mmle.tex).  Estimates the P-target channel parameters
% eta = [aR_p; aI_p; tau_p; nu_p] from the noisy delay-Doppler observation y,
% using the DETECTED delay-Doppler symbols xhat as the despreading reference.
%
% By the (I)SFFT unitarity (cp_otfs_mmle.tex, "Equivalence with the OFDM matched
% filter", eq. 16), the DD observation and the DD steering vector are SFFT images
% of their time-frequency counterparts, y^DD = U Y^TF and c^DD = U c^TF with U
% the unitary SFFT.  Because U^H U = I, both the matched-filter correlation
% c^H y and the norm ||c||^2 are invariant under the transform, so the CP-OTFS
% MMLE yields the SAME metric as the OFDM MMLE and may be evaluated in the TF
% grid.  We therefore map the DD grids to TF,
%   Y^TF = ISFFT(y^DD),   Xhat^TF = ISFFT(xhat^DD),
% and delegate to mmle_OFDM -- exactly the wrapper pattern used by
% pseudoTrue_CP_OTFS / getAB_CP_OTFS_*.  The channel parameter eta is the same
% object in both domains.
%
% Inputs:
%   y       : N x M noisy delay-Doppler observation (true mean + noise)
%   xhat    : N x M detected delay-Doppler symbol grid (despreading reference)
%   deltaf  : subcarrier spacing [Hz]
%   Ts      : OFDM symbol duration incl. CP [s]
%   P       : number of targets to extract
%   use_ici : true -> OFDM-side ICI kernel; false -> no-ICI
%   Zpad    : (optional) FFT zero-pad factor for the coarse search; default 4
%
% Outputs:
%   eta_hat              : 4P x 1 stacked estimate [aR;aI;tau;nu] per target
%   alpha_hat,tau_hat,nu_hat : P x 1 gains / delays / Dopplers

    if nargin < 7, Zpad = []; end
    Y    = isfft_u(y);        % Y^TF = ISFFT(y^DD)      (U^{-1} = ISFFT)
    Xhat = isfft_u(xhat);     % Xhat^TF = ISFFT(xhat^DD)
    [eta_hat, alpha_hat, tau_hat, nu_hat] = ...
        mmle_OFDM(Y, Xhat, deltaf, Ts, P, use_ici, Zpad);
end
