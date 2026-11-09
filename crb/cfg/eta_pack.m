function eta = eta_pack(alpha, tau, nu)
% ETA_PACK  Stack per-target parameters into eta = [eta_1; ...; eta_P],
% eta_p = [alpha_R; alpha_I; tau; nu]  (4P x 1).   Inverse: eta_unpack.
    alpha = alpha(:);  tau = tau(:);  nu = nu(:);
    P = numel(alpha);
    assert(numel(tau) == P && numel(nu) == P, 'eta_pack:size', 'alpha, tau, nu must have equal length.');
    eta = reshape([real(alpha).'; imag(alpha).'; tau.'; nu.'], [], 1);
end
