function [alpha, tau, nu] = unpackEta(eta)
% Unpack a stacked parameter vector into per-target gains/delays/Dopplers.
%   eta = [aR_1; aI_1; tau_1; nu_1;  ...;  aR_P; aI_P; tau_P; nu_P]   (4P x 1)
% -> alpha = aR + j aI (P x 1 complex), tau (P x 1), nu (P x 1).
% Inverse of the stacking used by getAB_OFDM_core / pseudoTrue_OFDM (eta_0).

    eta = eta(:);
    if mod(numel(eta), 4) ~= 0
        error('unpackEta:size', 'eta length must be a multiple of 4 (got %d).', numel(eta));
    end
    alpha = eta(1:4:end) + 1j*eta(2:4:end);
    tau   = eta(3:4:end);
    nu    = eta(4:4:end);
end
