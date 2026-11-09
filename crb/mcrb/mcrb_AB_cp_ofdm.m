function [A, B, Jt, b, Adelta, G, Delta] = mcrb_AB_cp_ofdm(X, Xhat, g, cr_eval, use_ici, sigma2, cr_true)
% MCRB_AB_CP_OFDM  Per-realization MCRB component matrices for CP-OFDM
% (latex/mcrb/MCRB derivations.tex; cp_ofdm_mcrb.tex Secs. B and A incl. the
% review Finding-5 gain-kinematic entries).
%
%   [A, B, Jt, b, Adelta] = mcrb_AB_cp_ofdm(X, Xhat, g, cr_eval, use_ici, sigma2, cr_true)
%
%   Delta   = mu(X, eta_true) - mu_tilde(Xhat, eta_eval)       (mean mismatch)
%   [Jt]_kl = (2/s2) Re{ g_k^H g_l }        g_k = d mu_tilde / d eta_k  (at Xhat, eta_eval)
%   [b]_k   = (2/s2) Re{ Delta^H g_k }
%   [A_D]_kl= (2/s2) Re{ Delta^H h_kl }     h_kl = d^2 mu_tilde / d eta_k d eta_l
%   B = Jt + b b^T,   A = -Jt + A_D         (A_D block-diagonal: cross-target Hessians vanish)
%
% cr_eval = evaluation point (normally the pseudo-true eta_0), cr_true = ground
% truth eta_bar with the transmitted X (defaults to cr_eval).  With Xhat = X and
% cr_eval = cr_true: b = 0, A = -J, B = J (matched CRB).  Hessian fields Xi,
% Omega1, Theta come from fields_cp_ofdm(order = 2) with the reduced I', I''.
% G (NM x 4P) and Delta (NM x 1) are returned for reuse.
% Ported from old/eva/utils/getAB_OFDM_core.m.
    if nargin < 7 || isempty(cr_true), cr_true = cr_eval; end
    [N, M] = size(X);
    P = numel(cr_eval.alpha);
    mu_true = io_cp_ofdm(X, g, cr_true, use_ici);
    mu_true = mu_true(:);
    G = zeros(N*M, 4*P);  mu_as = zeros(N*M, 1);
    Fc = cell(P, 1);
    for p = 1:P
        a = cr_eval.alpha(p);
        F = fields_cp_ofdm(Xhat, g, cr_eval.tau(p), cr_eval.nu(p), use_ici, 2);
        col = 4*(p-1);
        G(:, col+1) =             F.Phi(:);
        G(:, col+2) = 1j       *  F.Phi(:);
        G(:, col+3) = -1j*2*pi*a * F.Psi(:);
        G(:, col+4) =  1j*2*pi*a * F.Omega(:);
        mu_as = mu_as + a*F.Phi(:);
        Fc{p} = F;
    end
    Delta = mu_true - mu_as;
    Jt = (2/sigma2) * real(G'*G);  Jt = (Jt + Jt.')/2;
    b  = (2/sigma2) * real(G'*Delta);
    B  = Jt + b*b.';
    Adelta = zeros(4*P);
    contract = @(h) (2/sigma2) * real(Delta' * h(:));
    for p = 1:P
        a = cr_eval.alpha(p);  F = Fc{p};
        i1 = 4*(p-1)+1;  i2 = i1+1;  i3 = i1+2;  i4 = i1+3;
        Adelta(i1,i3) = contract(-1j*2*pi   * F.Psi);     Adelta(i3,i1) = Adelta(i1,i3);
        Adelta(i1,i4) = contract( 1j*2*pi   * F.Omega);   Adelta(i4,i1) = Adelta(i1,i4);
        Adelta(i2,i3) = contract( 2*pi      * F.Psi);     Adelta(i3,i2) = Adelta(i2,i3);
        Adelta(i2,i4) = contract(-2*pi      * F.Omega);   Adelta(i4,i2) = Adelta(i2,i4);
        Adelta(i3,i3) = contract(-4*pi^2*a  * F.Xi);
        Adelta(i4,i4) = contract(-4*pi^2*a  * F.Theta);
        Adelta(i3,i4) = contract( 4*pi^2*a  * F.Omega1);  Adelta(i4,i3) = Adelta(i3,i4);
    end
    A = -Jt + Adelta;
end
