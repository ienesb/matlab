function y = io_cp_otfs(x, g, cr, use_ici)
% IO_CP_OTFS  Noise-free CP-OTFS delay-Doppler input/output relation
% (latex/io/cp_otfs_io.tex): ISFFT -> per-symbol CP-OFDM channel -> SFFT.
% Equals (1/NM) sum_p alpha_p sum_{l,k} x[l,k] K_p(l',l) B_p(k-k') with the
% two-index kernel K_p (ICI) or A_p(l'-l) (no ICI), by unitarity of the
% (I)SFFT.  x : N x M DD grid (delay rows, Doppler cols).
    y = sfft_dd(io_cp_ofdm(isfft_dd(x), g, cr, use_ici));
end
