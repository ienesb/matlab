function Y = io_rcp_ofdm(X, g, cr, use_ici)
% IO_RCP_OFDM  Noise-free RCP-OFDM (single frame CP) time-frequency relation,
% latex/io/rcp_ofdm_io.tex eq. RCP_OFDM_IO:
%   Y[n,m] = sum_p alpha_p e^{j2pi nu_p m T} sum_{n'} e^{-j2pi n' df tau_p}
%            ( X[n',m] I_p^S(n'-n) + X[n',(m-1)_M] I_p^ISI(n'-n) )
% Symbol period T (no per-symbol CP); ISI from the previous (circular) symbol.
% use_ici=false keeps the delay split but freezes the kernels at nu = 0.
    Y = zeros(size(X));
    for p = 1:numel(cr.alpha)
        F = fields_rcp_ofdm(X, g, cr.tau(p), cr.nu(p), use_ici);
        Y = Y + cr.alpha(p) * F.Phi;
    end
end
