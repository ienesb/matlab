function [alpha, tau, nu] = eta_unpack(eta)
% ETA_UNPACK  Inverse of eta_pack: eta (4P x 1) -> alpha (P x 1 complex), tau, nu.
    eta = eta(:);
    assert(mod(numel(eta), 4) == 0, 'eta_unpack:size', 'eta length must be a multiple of 4 (got %d).', numel(eta));
    alpha = eta(1:4:end) + 1j*eta(2:4:end);
    tau   = eta(3:4:end);
    nu    = eta(4:4:end);
end
