function J = getFIM_CP_OTFS_noICI(x, deltaf, Ts, sigma2, alpha, tau, nu)
% Single-target CP-OTFS Fisher information matrix without OFDM-side ICI
% (Section 3.1.1 of cp_otfs_fim.tex).  Parameter ordering:
%   eta = [alphaR, alphaI, tau, nu].
%
% CP-OTFS = OTFS with a cyclic prefix appended after every OFDM-like symbol
% (M CPs per frame).  A reduced-CP variant (RCP-OTFS, single CP per frame)
% will be implemented separately.
%
% Same unitary argument as getFIM_CP_OTFS_ICI: J^CP-OTFS(x) at a realization
% equals J^OFDM_noICI(X) at X = ISFFT(x).

    X = isfft_u(x);
    J = getFIM_OFDM_noICI(X, deltaf, Ts, sigma2, alpha, tau, nu);
end
