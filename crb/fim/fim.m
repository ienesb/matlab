function [J, G] = fim(S, wf, g, cr, use_ici, sigma2)
% FIM  Per-realization Fisher information matrix (latex/fim/fim.tex eq. 9):
%   J = (2/sigma2) Re{ G^H G },  G = grad_stack(S, wf, g, cr, use_ici)   (4P x 4P)
% S is the transmitted grid in the waveform's native domain (TF for OFDM,
% DD for OTFS).  Genie CRB = J^{-1}; use crb_from_fim for a conditioned inverse.
    if nargin < 6 || isempty(sigma2), sigma2 = 1; end
    G = grad_stack(S, wf, g, cr, use_ici);
    J = (2/sigma2) * real(G' * G);
    J = (J + J.')/2;
end
