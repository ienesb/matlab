function pc = mcrb_assemble(A, Jt, b, Jgenie, eta_bar, eta0, g, active)
% MCRB_ASSEMBLE  Affine-in-noise pieces of the mismatched bounds at sigma2 = 1,
% inverted in bin units (stab_quad_bin):
%   G  = A^-1 Jt A^-1     (noise term)          F  = A^-1 b b^T A^-1   (score-bias floor)
%   BB = (eta_bar - eta0)(eta_bar - eta0)^T     Gc = Jgenie^-1          (genie CRB)
% so that for noise N0:  MCRB0 = N0 G + F,  LB = MCRB0 + BB,  CRB = N0 Gc
% (mcrb_eval).  A, Jt, b must have been built at sigma2 = 1 (J and b scale as
% 1/sigma2, so A^-1 B A^-1 has exactly one sigma2^1 and one sigma2^0 piece;
% eta0 is a noiseless projection and does not depend on sigma2).
% active: logical 4P mask (default all; MT40 masks nu).  pc.kappa = cond of the
% bin-unit A, pc.valid = finite and cond within 1e12.
    P = numel(eta_bar)/4;
    D = eta_scale(g, P);
    if nargin < 8 || isempty(active), active = true(4*P, 1); end
    [pc.G, pc.kappa] = stab_quad_bin(A, Jt, D, active);
    pc.F  = stab_quad_bin(A, b*b.', D, active);
    dlt = eta_bar(:) - eta0(:);  dlt(~active) = 0;
    pc.BB = dlt*dlt.';
    pc.Gc = stab_quad_bin(-Jgenie, Jgenie, D, active);
    pc.active = active;  pc.P = P;
    pc.valid = all(isfinite(pc.G(:))) && all(isfinite(pc.F(:))) && pc.kappa < 1e12;
    pc.eta_bar = eta_bar(:);  pc.eta0 = eta0(:);
end
