function [Y, kappa] = stab_quad_bin(A, C, D, active)
% STAB_QUAD_BIN  Y = A^{-1} C A^{-1} evaluated in bin units, back in SI.
%   [Y, kappa] = stab_quad_bin(A, C, D, active)
% D = eta_scale(g, P); active = logical mask of estimated parameters (inactive
% rows/cols are structurally null and are excluded before inversion).  The
% scaled A_s = D A D is symmetrized and pseudo-inverted; kappa = cond(A_s).
% Replaces the five duplicated stab_quad copies of old/eva (Jacobi variant).
    n = size(A, 1);
    if nargin < 4 || isempty(active), active = true(n, 1); end
    As = D*A*D;  Cs = D*C*D;
    Aa = As(active, active);  Aa = (Aa + Aa.')/2;
    Ca = Cs(active, active);
    Ai = pinv(Aa);
    Ya = real(Ai * Ca * Ai);
    Yb = zeros(n);  Yb(active, active) = Ya;
    Y  = D*Yb*D;
    kappa = cond(Aa);
end
