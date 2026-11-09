function G = grad_stack(S, wf, g, cr, use_ici)
% GRAD_STACK  Gradient stack G (NM x 4P), column 4(p-1)+k = vec(d mu / d eta_{p,k}),
% eta_p = [alpha_R; alpha_I; tau; nu], for one data realization S (N x M).
%
%   G = grad_stack(S, wf, g, cr, use_ici)
%   wf: 'cp-ofdm' | 'cp-otfs' | 'rcp-ofdm' | 'rcp-otfs' | 'rcp-otfs-cont'
%
% cp-ofdm       [Phi, jPhi, -j2pi a Psi, j2pi a Omega]           (cp_ofdm_fim.tex eqs 5-8;
%               Psi carries n' df, Omega = m TS Phi + I' term)
% cp-otfs       grad_stack(isfft_dd(S), 'cp-ofdm')               (cp_otfs_fim.tex: exact by unitarity)
% rcp-ofdm      [Phi, jPhi, a Dtau, a Dnu]  with the exact boundary term (rcp_ofdm_fim.tex)
% rcp-otfs      sampled branch [Phi, jPhi, -j2pi df a Psi, j2pi a Obar] (rcp_otfs_fim.tex;
%               integer delays only)
% rcp-otfs-cont grad_stack(isfft_dd(S), 'rcp-ofdm')              (rcp_ofdm_fim.tex line 319,
%               J^cont_RCP-OTFS(x) = J_RCP-OFDM(ISFFT x), exact per realization; carries the
%               boundary field Xi exactly)
%
% Domain note: for 'cp-otfs' and 'rcp-otfs-cont' the columns are the TF-domain
% fields of X = ISFFT(x) (not their SFFT).  The FIM and every MCRB ingredient
% are unitarily invariant, so this is exact; anything that pairs G with a
% mean mismatch must build the mismatch in the same (TF) domain.
% FIM: J = (2/sigma2) real(G'*G)  (fim.m).  Ported from old/eva/utils/getFIM_*_multi_*.m.
    [N, M] = size(S);
    P = numel(cr.alpha);
    wf = lower(wf);
    switch wf
        case 'cp-otfs'
            G = grad_stack(isfft_dd(S), 'cp-ofdm', g, cr, use_ici);  return;
        case 'rcp-otfs-cont'
            G = grad_stack(isfft_dd(S), 'rcp-ofdm', g, cr, use_ici); return;
    end
    G = zeros(N*M, 4*P);
    for p = 1:P
        a = cr.alpha(p);  col = 4*(p-1);
        switch wf
            case 'cp-ofdm'
                F = fields_cp_ofdm(S, g, cr.tau(p), cr.nu(p), use_ici, 1);
                G(:, col+1) =              F.Phi(:);
                G(:, col+2) = 1j        *  F.Phi(:);
                G(:, col+3) = -1j*2*pi*a * F.Psi(:);
                G(:, col+4) =  1j*2*pi*a * F.Omega(:);
            case 'rcp-ofdm'
                F = fields_rcp_ofdm(S, g, cr.tau(p), cr.nu(p), use_ici);
                G(:, col+1) =      F.Phi(:);
                G(:, col+2) = 1j * F.Phi(:);
                G(:, col+3) = a  * F.Dtau(:);
                G(:, col+4) = a  * F.Dnu(:);
            case 'rcp-otfs'
                F = fields_rcp_otfs(S, g, cr.tau(p), cr.nu(p), use_ici);
                G(:, col+1) =                       F.Phi(:);
                G(:, col+2) = 1j                 *  F.Phi(:);
                G(:, col+3) = -1j*2*pi*g.deltaf*a * F.Psi(:);
                G(:, col+4) =  1j*2*pi*a          * F.Obar(:);
            otherwise
                error('grad_stack:wf', 'unknown waveform %s', wf);
        end
    end
end
