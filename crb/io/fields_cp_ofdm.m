function F = fields_cp_ofdm(S, g, tau, nu, use_ici, order)
% FIELDS_CP_OFDM  Per-target derivative fields of the CP-OFDM mean
% (latex/fim/cp_ofdm_fim.tex eqs. 5-8, 60-62; latex/mcrb/cp_ofdm_mcrb.tex
% shorthand fields).  Single target (scalar tau, nu); callers loop over p.
%
%   F = fields_cp_ofdm(S, g, tau, nu, use_ici, order)
%
% With  A[n',m] = S[n',m] e^{-j2pi n' df tau} e^{+j2pi m TS nu}  and the
% 1/T-normalized ICI kernel I(n'-n) (delta[n'-n] when use_ici = false):
%   Phi    = sum_{n'} A I                          (mean field, d mu/d alpha_R)
%   Psi    = sum_{n'} (n' df)   A I                (d mu/d tau = -j2pi alpha Psi)
%   Omega  = m TS Phi + sum_{n'} A I'              (d mu/d nu  = +j2pi alpha Omega)
%   Xi     = sum_{n'} (n' df)^2 A I                (order 2, Hessian tau,tau)
%   Omega1 = m TS Psi + sum_{n'} (n' df) A I'      (order 2, Hessian tau,nu)
%   Theta  = (m TS)^2 Phi + 2 m TS sum A I' + sum A I''   (order 2, Hessian nu,nu)
% where I' := (1/j2pi) dI/dnu, I'' := (1/j2pi) dI'/dnu (ici_kernel).
%
% order = 0 returns only Phi (used by io_cp_ofdm); 1 adds Psi, Omega;
% 2 adds Xi, Omega1, Theta.  All fields are N x M.
%
% Ported from old/eva/utils/getAB_OFDM_core.m and getFIM_OFDM_multi_ICI.m
% (this is the single shared field builder; the four legacy ofdm_mean copies
% and the FIM/MCRB/pseudo-true builders all reduce to it).
    if nargin < 6, order = 1; end
    [N, M] = size(S);
    ns = (0:N-1).';
    ms = (0:M-1);
    ndf = ns*g.deltaf;              % n' df   (N x 1)
    mTS = ms*g.TS;                  % m TS    (1 x M)

    A = S .* (exp(-1j*2*pi*ndf*tau) * exp(1j*2*pi*mTS*nu));   % N x M

    if use_ici
        [NP, NN] = ndgrid(0:N-1, 0:N-1);
        [Imat, Ip, Ipp] = ici_kernel(NP - NN, nu, g.T);   % Imat(n'+1,n+1) = I(n'-n)
        F.Phi = Imat.' * A;
        if order >= 1
            F.Psi   = Imat.' * (ndf .* A);
            D1      = Ip.' * A;
            F.Omega = mTS .* F.Phi + D1;
        end
        if order >= 2
            F.Xi     = Imat.' * (ndf.^2 .* A);
            D1n      = Ip.'  * (ndf .* A);
            D2       = Ipp.' * A;
            F.Omega1 = mTS .* F.Psi + D1n;
            F.Theta  = mTS.^2 .* F.Phi + 2*mTS .* D1 + D2;
        end
    else
        F.Phi = A;
        if order >= 1
            F.Psi   = ndf .* A;
            F.Omega = mTS .* A;
        end
        if order >= 2
            F.Xi     = ndf.^2 .* A;
            F.Omega1 = (ndf .* mTS) .* A;
            F.Theta  = mTS.^2 .* A;
        end
    end
end
