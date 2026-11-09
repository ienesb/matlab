function J = getFIM_CP_OTFS_ICI(x, deltaf, Ts, sigma2, alpha, tau, nu)
% Single-target CP-OTFS Fisher information matrix with OFDM-side ICI
% (Section 3.1.2 of cp_otfs_fim.tex).  Parameter vector ordering:
%   eta = [alphaR, alphaI, tau, nu].
%
% CP-OTFS = OTFS with a cyclic prefix appended after every OFDM-like symbol
% (M CPs per frame).  A reduced-CP variant (RCP-OTFS, single CP per frame)
% will be implemented separately.
%
% Inputs:
%   x       : N x M delay-Doppler symbol grid (delay rows, Doppler cols)
%   deltaf  : subcarrier spacing [Hz]
%   Ts      : OFDM symbol duration including CP [s]
%   sigma2  : noise variance
%   alpha   : complex path gain
%   tau     : path delay [s]
%   nu      : path Doppler [Hz]
%
% Implementation: the (I)SFFT is unitary, so per realization the CP-OTFS
% gradient stack G_y equals (unitary operator) * G_Y.  Therefore
% G_y^H G_y = G_Y^H G_Y exactly, i.e. the CP-OTFS-with-OFDM-ICI FIM at a
% data realization x equals the OFDM-with-ICI FIM at X = ISFFT(x).
% This is the per-realization expression of the equivalence
% E{J^CP-OTFS} = E{J^OFDM} discussed in cp_otfs_fim.tex.

    X = isfft_u(x);
    J = getFIM_OFDM_ICI(X, deltaf, Ts, sigma2, alpha, tau, nu);
end
