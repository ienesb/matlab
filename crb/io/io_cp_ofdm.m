function Y = io_cp_ofdm(X, g, cr, use_ici)
% IO_CP_OFDM  Noise-free CP-OFDM time-frequency input/output relation
% (latex/io/cp_ofdm_io.tex, boxed ICI equation; no-ICI form when use_ici=false):
%
%   Y[n,m] = sum_p alpha_p sum_{n'} X[n',m] e^{-j2pi n' df tau_p}
%                                          e^{ j2pi m  TS nu_p } I_p(n'-n)
%
% X : N x M transmitted TF grid;  g : make_grid struct (TS = T + Tcp);
% cr: channel realization (alpha, tau [s], nu [Hz], P x 1).
% Add noise outside:  Y + sqrt(N0/2)*(randn + 1j randn).
    Y = zeros(size(X));
    for p = 1:numel(cr.alpha)
        F = fields_cp_ofdm(X, g, cr.tau(p), cr.nu(p), use_ici, 0);
        Y = Y + cr.alpha(p) * F.Phi;
    end
end
