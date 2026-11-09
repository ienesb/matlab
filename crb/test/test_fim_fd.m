function r = test_fim_fd()
% Finite-difference check of every gradient column for every waveform (Phase 1).
% Perturbations in bin units (h*Ts for tau, h*bin_nu for nu), central differences.
    r = [];
    N = 16; M = 8; df = 15e3;  Ts = 1/(N*df);
    g  = make_grid(N, M, df, 3*Ts, 'cp');
    gr = make_grid(N, M, df, 3*Ts, 'rcp');
    rs = crn_stream(11, 0, 1);
    S  = gen_symbols(N, M, '4qam', rs);
    cr.P = 2;  cr.alpha = [0.9*exp(1j*0.3); 0.6*exp(-1j*1.2)];
    cr.nu = [0.3; -0.45]/(M*g.TS);
    h = 1e-4;
    cases = { 'cp-ofdm', g,  false, [1.3; 2.6]*Ts;   'cp-ofdm', g,  true, [1.3; 2.6]*Ts;
              'cp-otfs', g,  true,  [1.3; 2.6]*Ts;   'rcp-ofdm', gr, false, [1.3; 2.6]*Ts;
              'rcp-ofdm', gr, true, [1.3; 2.6]*Ts;   'rcp-otfs', gr, true,  [1; 3]*Ts;
              'rcp-otfs', gr, false, [1; 3]*Ts;      'rcp-otfs-cont', gr, true, [1.3; 2.6]*Ts };
    for ic = 1:size(cases, 1)
        wf = cases{ic,1};  gg = cases{ic,2};  ici = cases{ic,3};  cr.tau = cases{ic,4};
        cr.eta_bar = eta_pack(cr.alpha, cr.tau, cr.nu);
        G = grad_stack(S, wf, gg, cr, ici);
        mu = @(c) mean_of(S, wf, gg, c, ici);
        worst = zeros(1, 4);
        for p = 1:cr.P
            for k = 1:4
                cp = cr;  cm = cr;  scale = 1;
                switch k
                    case 1, cp.alpha(p) = cp.alpha(p) + h;     cm.alpha(p) = cm.alpha(p) - h;
                    case 2, cp.alpha(p) = cp.alpha(p) + 1j*h;  cm.alpha(p) = cm.alpha(p) - 1j*h;
                    case 3, cp.tau(p) = cp.tau(p) + h*Ts;      cm.tau(p) = cm.tau(p) - h*Ts;      scale = Ts;
                    case 4, cp.nu(p)  = cp.nu(p)  + h*gg.bin_nu; cm.nu(p) = cm.nu(p) - h*gg.bin_nu; scale = gg.bin_nu;
                end
                fd = (mu(cp) - mu(cm)) / (2*h);            % d mu / d eta_bin
                an = G(:, 4*(p-1)+k) * scale;
                worst(k) = max(worst(k), relerr(fd, an));
            end
        end
        r = tc(r, sprintf('FD %-14s ici=%d  [aR aI tau nu] max', wf, ici), max(worst), 1e-6);
        fprintf('      per column: %.1e %.1e %.1e %.1e\n', worst);
    end
end

function m = mean_of(S, wf, g, cr, ici)
% Mean in the domain of the gradient stack: the OTFS-by-unitarity waveforms
% ('cp-otfs', 'rcp-otfs-cont') keep their columns in the TF domain of
% ISFFT(x); the FIM is invariant, but a column-wise FD check must use the
% same domain.
    switch wf
        case 'cp-ofdm',       m = io_cp_ofdm(S, g, cr, ici);
        case 'cp-otfs',       m = io_cp_ofdm(isfft_dd(S), g, cr, ici);
        case 'rcp-ofdm',      m = io_rcp_ofdm(S, g, cr, ici);
        case 'rcp-otfs'
            m = zeros(size(S));
            for p = 1:numel(cr.alpha)
                F = fields_rcp_otfs(S, g, cr.tau(p), cr.nu(p), ici, false);   % off-grid allowed for FD
                m = m + cr.alpha(p)*F.Phi;
            end
        case 'rcp-otfs-cont', m = io_rcp_ofdm(isfft_dd(S), g, cr, ici);
    end
    m = m(:);
end
