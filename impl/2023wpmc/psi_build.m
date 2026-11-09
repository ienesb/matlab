function Psi = psi_build(scheme, h, tau, nu, sys)
% PSI_BUILD  Dispatch to the per-scheme closed-form channel matrix Psi^p.
%   scheme in {'cpfree_ofdm','cp_ofdm','zak_otfs','twosteps_otfs'}.
switch scheme
    case 'cpfree_ofdm',   Psi = psi_cpfree_ofdm(h,tau,nu,sys);
    case 'cp_ofdm',       Psi = psi_cp_ofdm(h,tau,nu,sys);
    case 'zak_otfs',      Psi = psi_zak(h,tau,nu,sys);
    case 'twosteps_otfs', Psi = psi_twosteps(h,tau,nu,sys);
    otherwise, error('psi_build:scheme','unknown scheme %s',scheme);
end
end
