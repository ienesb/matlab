function [A, B, Jt, b, Adelta, G, Delta] = mcrb_AB_cp_otfs(x, xhat, g, cr_eval, use_ici, sigma2, cr_true)
% MCRB_AB_CP_OTFS  CP-OTFS MCRB component matrices (latex/mcrb/cp_otfs_mcrb.tex):
% by (I)SFFT unitarity every ingredient equals the CP-OFDM one evaluated at
% X = ISFFT(x), Xhat = ISFFT(xhat), exactly per realization.  The BSC/detector
% errors live in the delay-Doppler grid (xhat), which is what separates the
% waveforms' fourth-order statistics.  RCP waveforms have no MCRB derivation.
    if nargin < 7, cr_true = []; end
    [A, B, Jt, b, Adelta, G, Delta] = mcrb_AB_cp_ofdm(isfft_dd(x), isfft_dd(xhat), g, cr_eval, use_ici, sigma2, cr_true);
end
