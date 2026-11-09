function y = io_rcp_otfs(x, g, cr, use_ici)
% IO_RCP_OTFS  Noise-free sampled-receiver RCP-OTFS delay-Doppler relation
% (latex/io/rcp_otfs_io.tex eq. RCP_IO_general; integer delays required).
% use_ici=true: exact (intra-frame Doppler phase kept); false: eq. RCP_IO_noICI.
    y = zeros(size(x));
    for p = 1:numel(cr.alpha)
        F = fields_rcp_otfs(x, g, cr.tau(p), cr.nu(p), use_ici);
        y = y + cr.alpha(p) * F.Phi;
    end
end
