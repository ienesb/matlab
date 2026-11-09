function J = getFIM_CP_OTFS_multi_noICI(x, deltaf, Ts, sigma2, alpha, tau, nu)
% Multi-target CP-OTFS Fisher information matrix, no OFDM-side ICI
% (Section 3.2.1 of cp_otfs_fim.tex).  Parameter ordering:
%   eta = [a1R, a1I, tau1, nu1,  a2R, a2I, tau2, nu2,  ...]
%
% CP-OTFS = OTFS with a cyclic prefix appended after every OFDM-like symbol
% (M CPs per frame).  A reduced-CP variant (RCP-OTFS, single CP per frame)
% will be implemented separately.
%
% Implementation: per-realization equivalent to OFDM no-ICI FIM at
% X = ISFFT(x) (unitarity argument).

    X = isfft_u(x);
    J = getFIM_OFDM_multi_noICI(X, deltaf, Ts, sigma2, alpha, tau, nu);
end
