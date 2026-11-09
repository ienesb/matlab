function D = eta_scale(g, P)
% ETA_SCALE  Diagonal SI -> bin-unit scaling for the stacked 4P parameter vector.
%
% Bin-unit parameters  eta_bin = [alpha_R; alpha_I; tau/Ts; nu*M*TS]  per target,
% i.e. eta = D * eta_bin with  D = diag([1 1 Ts 1/(M*TS)]) per target.
% Hence  J_bin = D*J*D  and  CRB = D * J_bin^{-1} * D  (exact change of
% variables, no approximation).  In SI units a joint (tau, nu) FIM has
% condition numbers ~1e20 (1/s^2 vs 1/Hz^2 entries); in bin units it is O(1).
    d = repmat([1; 1; g.bin_tau; g.bin_nu], P, 1);
    D = diag(d);
end
