function F = fields_rcp_ofdm(X, g, tau, nu, use_ici)
% FIELDS_RCP_OFDM  Per-target mean and derivative fields of the RCP-OFDM
% (single frame CP) relation, latex/io/rcp_ofdm_io.tex eq. RCP_OFDM_IO and
% latex/fim/rcp_ofdm_fim.tex (derivative fields Phi, Theta, Upsilon):
%
%   mu[n,m] = alpha e^{j2pi nu m T} sum_{n'} e^{-j2pi n' df tau}
%             ( X[n',m] I^S(n'-n) + X[n',(m-1)_M] I^ISI(n'-n) )
%
% with the split kernels I^ISI + I^S = I (isi_kernel / ici_kernel) and the
% useful duration T = 1/df in the symbol-time phase (no per-symbol CP).
%
% Outputs (N x M):
%   F.Phi   mean field                (d mu/d alpha_R = Phi, d mu/d alpha_I = j Phi)
%   F.Dtau  d mu/d tau  / alpha       includes the exact boundary term
%                                     (1/T) e^{j2pi((n'-n)df+nu)tau} (X_prev - X)
%   F.Dnu   d mu/d nu   / alpha       window phase term + split-kernel derivatives
%
% use_ici = false: kernels frozen at nu = 0, their nu-derivatives dropped,
% the delay-induced split RETAINED (README: "no ICI" removes only the
% Doppler-induced leakage), symbol-time phase e^{j2pi nu m T} kept.
%
% Ported from old/eva/utils/getFIM_RCP_OFDM_multi.m (kernel derivatives now
% come from ici_kernel/isi_kernel in reduced form; Dnu carries the j2pi).
    [N, M] = size(X);
    T  = g.T;
    ns = (0:N-1).';
    ms = (0:M-1);
    if use_ici, nuk = nu; else, nuk = 0; end

    [NP, NN] = ndgrid(0:N-1, 0:N-1);
    q = NP - NN;                                   % lag n' - n, rows = n'
    [I_full, Ip_full] = ici_kernel(q, nuk, T);     % I, reduced I'
    [I_isi, Ip_isi, Ebnd] = isi_kernel(q, nuk, T, tau);
    I_s = I_full - I_isi;

    Ap      = X .* exp(-1j*2*pi*g.deltaf*ns*tau);  % delay phase only
    Ap_prev = circshift(Ap, 1, 2);                 % column m <- column (m-1) mod M
    ph      = exp(1j*2*pi*T*ms*nu);                % per-window Doppler phase (1 x M)

    F.Phi  = ph .* (I_s.' * Ap + I_isi.' * Ap_prev);

    F.Dtau = ph .* ( -1j*2*pi*g.deltaf * (I_s.' * (ns .* Ap) + I_isi.' * (ns .* Ap_prev)) ...
                     + Ebnd.' * (Ap_prev - Ap) );

    F.Dnu  = (1j*2*pi*T*ms) .* F.Phi;
    if use_ici
        Ip_s  = Ip_full - Ip_isi;
        F.Dnu = F.Dnu + 1j*2*pi * ph .* (Ip_s.' * Ap + Ip_isi.' * Ap_prev);
    end
end
