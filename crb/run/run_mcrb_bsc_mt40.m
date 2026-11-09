function res = run_mcrb_bsc_mt40(varargin)
% RUN_MCRB_BSC_MT40  MCRB / LB vs error rate p (and vs SNR) on MT40 (plan F3.3),
% delay-only targets (nu masked), no-ICI model, 16-PSK (Gray-bit BSC) and
% 16-QAM (uniform symbol error).  Two framings:
%   '2d' (default): MT40-2D frame N = 128 x M = 16, Tcp = 127 Ts, {cp-ofdm, cp-otfs}
%   '1d':           MT40-1D block N = 128, M = 1, Tcp = 0 with Liu's unitary
%                   precoders X = Q^H s: {OFDM, OTFS (N1 = 32, N2 = 4), SC}
% L = o.L (10 by default) for the p sweep; when o.L ~= 40 the p = 0 point is ALSO
% computed at L = 40 (res.L40_p0) and cross-checked against the stored F1.1
% (results/crb_mt40_block.mat, 1-D) / F1.3 (results/crb_mt40_frame_*.mat, 2-D)
% genie CRBs.  Statistics: median (headline) and mean over valid draws of the
% PER-DRAW composite bounds (mcrb_stats); kappa filter; pilot-only CRB (20 %
% comb, omitted when 3L > 2 Npil) and the Jensen bound as references; the
% data-vs-pilot crossover is printed at 0 and 20 dB.
%
%   res = run_mcrb_bsc_mt40('tier', 'smoke', 'framing', '1d', ...)
% Options: tier smoke(20)|dev(200)|paper(2000) draws or a number; framing;
% L 10; p_list (0:0.025:0.1 smoke, 0:0.002:0.1 otherwise); mods {'16psk','16qam'};
% snr_db [-10:1:6 10 20]; seed 7000; kappa_thresh 1e4; Tcp_samples 127; M 16
% (2-D only); pilot_frac 0.20; tag ''.
    o = struct('tier', 'smoke', 'framing', '2d', 'L', 10, 'p_list', [], 'mods', {{'16psk', '16qam'}}, ...
        'snr_db', [-10:1:6 10 20], 'seed', 7000, 'kappa_thresh', 1e4, 'Tcp_samples', 127, 'M', 16, 'pilot_frac', 0.20, 'tag', '');
    o = parse_opts(o, varargin, 'run_mcrb_bsc_mt40');
    root = fileparts(fileparts(mfilename('fullpath')));
    nT = tier_trials(o.tier, [20 200 2000]);
    if isempty(o.p_list), if strcmpi(o.tier, 'smoke'), o.p_list = 0:0.025:0.1; else, o.p_list = 0:0.002:0.1; end, end
    if isnumeric(o.tier), o.tier = sprintf('n%d', nT); end
    N = 128;  B = 160e6;  Ts = 1/B;
    is1d = strcmpi(o.framing, '1d');
    if is1d
        M = 1;  g = make_grid(N, 1, B/N, 0, 'cp');  Tcp_samples = 0;
        Qh = liu_precoders(N, 32, 4, 1/16, 1/8);  wfs = {'OFDM', 'OTFS', 'SC'};  fr = 'MT40-1D';
    else
        M = o.M;  g = make_grid(N, M, B/N, o.Tcp_samples*Ts, 'cp');  Tcp_samples = o.Tcp_samples;
        Qh = [];  wfs = {'cp-ofdm', 'cp-otfs'};  fr = 'MT40-2D';
    end
    ch = chan_mt40(g, 'L', o.L);  cr = chan_realize(ch, g, []);
    idx = eta_index(cr.P);  active = true(4*cr.P, 1);  active(idx.nu) = false;
    nP = numel(o.p_list);  nW = numel(wfs);  nMod = numel(o.mods);
    Jj = fim(ones(N, M), 'cp-ofdm', g, cr, false, 1);  Sj = crb_from_fim(Jj, g, cr);
    res = struct();
    res.jensen_tr_tau = Sj.tr_tau;                                   % sigma2 = 1
    Sp = crb_pilot_only(g, cr, false, o.pilot_frac, 1, 'cp-ofdm');  res.pilot_tr_tau = Sp.tr_tau;
    % a comb of Npil identical pilot subcarriers on every symbol gives only 2*Npil
    % real equations (no Doppler): with 3L unknowns the pilot-only FIM is rank
    % deficient for 3L > 2*Npil and pinv would return a meaningless truncated number
    res.pilot_identifiable = 3*cr.P <= 2*numel(Sp.pilot_idx);
    if ~res.pilot_identifiable
        fprintf('  pilot-only CRB not identifiable: 3L = %d unknowns > 2 x %d comb pilots -> omitted\n', 3*cr.P, numel(Sp.pilot_idx));
        res.pilot_tr_tau = NaN;
    end
    fld = {'G_tau', 'F_tau', 'BB_tau', 'Gc_tau', 'Gt_tau', 'Ft_tau', 'kappa', 'pcvalid', 'pt_iters', 'pt_sec', 'pt_exitflag'};
    fprintf('%s MCRB vs p: %s, L = %d, %d draws, p = %s, %s, waveforms %s\n', fr, ch.name, cr.P, nT, mat2str(o.p_list), strjoin(o.mods, '/'), strjoin(wfs, ' '));
    mods = o.mods;  p_list = o.p_list;  seed = o.seed;
    t0 = tic;
    out = cell(nT, 1);
    parfor t = 1:nT
        loc = struct();  for k = 1:numel(fld), loc.(fld{k}) = nan(nP, nW, nMod); end
        tr = @(Xm, i) sum(real(diag(Xm(i, i))));
        for im = 1:nMod
            mdm = modem(mods{im});
            errmode = 'bit';  if strcmpi(mods{im}, '16qam'), errmode = 'symbol'; end
            Sg = gen_symbols(N, M, mods{im}, crn_stream(seed, im, t));
            for iw = 1:nW
                if is1d, Q = Qh.(wfs{iw});  X = Q*Sg;  wf = 'cp-ofdm';  dom = 'tf';
                else,    Q = [];  X = Sg;  wf = wfs{iw};  dom = 'tf';  if strcmp(wf, 'cp-otfs'), dom = 'dd'; end, end
                Jg = fim(X, wf, g, cr, false, 1);
                eta_prev = [];
                for ip = 1:nP
                    p = p_list(ip);
                    Sh = bsc_symbols(Sg, p, mdm, crn_stream(seed, 100*im + ip, t), errmode);   % flips paired across waveforms
                    if is1d, Xh = Q*Sh; else, Xh = Sh; end
                    if p == 0
                        eta0 = cr.eta_bar;  pinfo = struct('iterations', 0, 'seconds', 0, 'exitflag', 1);
                    else
                        [eta0, pinfo] = pseudo_true(X, Xh, g, cr, false, struct('domain', dom, 'eta_init', eta_prev, 'p', p, 'mask_nu', true));
                    end
                    eta_prev = eta0;
                    [a0, t0v, n0] = eta_unpack(eta0);
                    ce = cr;  ce.alpha = a0;  ce.tau = t0v;  ce.nu = n0;
                    if strcmp(wf, 'cp-ofdm')
                        [A, ~, Jt, b]    = mcrb_AB_cp_ofdm(X, Xh, g, ce, false, 1, cr);
                        [At, ~, Jtt, bt] = mcrb_AB_cp_ofdm(X, Xh, g, cr, false, 1, cr);
                    else
                        [A, ~, Jt, b]    = mcrb_AB_cp_otfs(X, Xh, g, ce, false, 1, cr);
                        [At, ~, Jtt, bt] = mcrb_AB_cp_otfs(X, Xh, g, cr, false, 1, cr);
                    end
                    pc  = mcrb_assemble(A, Jt, b, Jg, cr.eta_bar, eta0, g, active);
                    pct = mcrb_assemble(At, Jtt, bt, Jg, cr.eta_bar, cr.eta_bar, g, active);
                    loc.G_tau(ip, iw, im) = tr(pc.G, idx.tau);   loc.F_tau(ip, iw, im) = tr(pc.F, idx.tau);
                    loc.BB_tau(ip, iw, im) = tr(pc.BB, idx.tau); loc.Gc_tau(ip, iw, im) = tr(pc.Gc, idx.tau);
                    loc.Gt_tau(ip, iw, im) = tr(pct.G, idx.tau); loc.Ft_tau(ip, iw, im) = tr(pct.F, idx.tau);
                    loc.kappa(ip, iw, im) = pc.kappa;            loc.pcvalid(ip, iw, im) = double(pc.valid);
                    loc.pt_iters(ip, iw, im) = pinfo.iterations;  loc.pt_sec(ip, iw, im) = pinfo.seconds;  loc.pt_exitflag(ip, iw, im) = pinfo.exitflag;
                end
            end
        end
        out{t} = loc;
        fprintf('  draw %d/%d (%.0f s)\n', t, nT, toc(t0));
    end
    for k = 1:numel(fld)
        A4 = nan(nT, nP, nW, nMod);
        for t = 1:nT, A4(t, :, :, :) = out{t}.(fld{k}); end
        res.(fld{k}) = A4;
    end
    res.valid = res.pcvalid > 0 & res.pt_exitflag > 0 & isfinite(res.kappa) & res.kappa <= o.kappa_thresh;
    res.drop_frac = reshape(1 - mean(res.valid, 1), nP, nW, nMod);
    res.p_list = o.p_list;  res.snr_db = o.snr_db;  res.opts = o;  res.opts.Tcp_samples = Tcp_samples;  res.nT = nT;  res.g = g;  res.ch = ch;  res.cr = cr;
    res.waveforms = wfs;  res.mods = o.mods;  res.framing = fr;  res.range_m2_per_samp2 = ch.range_unit_m^2;  res.seconds = toc(t0);
    res.statistic = sprintf('median (headline) and mean over valid draws of the per-draw composite bounds (kappa <= %g, converged pseudo-true)', o.kappa_thresh);
    % ---- p = 0 at L = 40 (F1.1 / F1.3 cross-check) --------------------------------
    if cr.P ~= 40
        ch40 = chan_mt40(g, 'L', 40);  cr40 = chan_realize(ch40, g, []);
        Gc40 = nan(nT, nW, nMod);  k40 = Gc40;
        parfor t = 1:nT
            gl = nan(nW, nMod);  kl = gl;
            for im = 1:nMod
                Sg = gen_symbols(N, M, mods{im}, crn_stream(seed, im, t));
                for iw = 1:nW
                    if is1d, X = Qh.(wfs{iw})*Sg;  wf = 'cp-ofdm'; else, X = Sg;  wf = wfs{iw}; end
                    S = crb_from_fim(fim(X, wf, g, cr40, false, 1), g, cr40);
                    gl(iw, im) = S.tr_tau;  kl(iw, im) = S.kappa;
                end
            end
            Gc40(t, :, :) = gl;  k40(t, :, :) = kl;
        end
        Sj40 = crb_from_fim(fim(ones(N, M), 'cp-ofdm', g, cr40, false, 1), g, cr40);
        res.L40_p0 = struct('Gc_tau', Gc40, 'kappa', k40, 'cr', cr40, 'jensen_tr_tau', Sj40.tr_tau);
    elseif any(o.p_list == 0)
        i0 = find(o.p_list == 0, 1);
        res.L40_p0 = struct('Gc_tau', reshape(res.Gc_tau(:, i0, :, :), nT, nW, nMod), 'kappa', reshape(res.kappa(:, i0, :, :), nT, nW, nMod), 'cr', cr, 'jensen_tr_tau', res.jensen_tr_tau);
    else
        res.L40_p0 = [];
    end
    print_summary(res, root);
    if ~exist(fullfile(root, 'results'), 'dir'), mkdir(fullfile(root, 'results')); end
    save(fullfile(root, 'results', sprintf('mcrb_bsc_mt40_%s_L%d_%s%s.mat', lower(strrep(fr, 'MT40-', '')), cr.P, o.tier, o.tag)), '-struct', 'res');
    plot_mcrb_mt40(res);
end

function print_summary(res, root)
    o = res.opts;  p = res.p_list;  nP = numel(p);  Ts = res.g.Ts;  c = res.range_m2_per_samp2;  nW = numel(res.waveforms);  nMod = numel(res.mods);
    i0 = find(p == 0, 1);  pm = p > 0;
    fprintf('\nAcceptance (F3.3, %s): delay trace bounds [m^2, sum over %d targets], median over valid draws (drop max %.1f%%)', res.framing, res.cr.P, 100*max(res.drop_frac(:)));
    if any(pm), fprintf('; pseudo-true mean %.1f it, %.3f s, converged %.1f%%', mean(res.pt_iters(:, pm, :, :), 'all'), mean(res.pt_sec(:, pm, :, :), 'all'), 100*mean(res.pt_exitflag(:, pm, :, :) > 0, 'all')); end
    fprintf('\n');
    for im = 1:nMod
        for iw = 1:nW
            S0  = mcrb_stats(res.G_tau(:, :, iw, im), res.F_tau(:, :, iw, im), res.BB_tau(:, :, iw, im), res.Gc_tau(:, :, iw, im), res.valid(:, :, iw, im), 1, 'median');
            T0  = mcrb_stats(res.Gt_tau(:, :, iw, im), res.Ft_tau(:, :, iw, im), 0*res.BB_tau(:, :, iw, im), res.Gc_tau(:, :, iw, im), res.valid(:, :, iw, im), 1, 'median');
            mcrb = S0.mcrb/Ts^2*c;  lb = S0.lb/Ts^2*c;  crb = S0.crb/Ts^2*c;  mt = T0.mcrb/Ts^2*c;
            fprintf('  %-6s %-8s p=%s\n     MCRB0 %s\n     LB    %s\n     MCRB@eta_bar %s\n     CRB %.4g  Jensen %.4g  pilot-only %.4g  drop %s\n', ...
                res.mods{im}, res.waveforms{iw}, mat2str(p, 2), mat2str(mcrb, 3), mat2str(lb, 3), mat2str(mt, 3), crb(max(i0, 1)), res.jensen_tr_tau/Ts^2*c, res.pilot_tr_tau/Ts^2*c, mat2str(res.drop_frac(:, iw, im).', 2));
            if ~isempty(i0), fprintf('     p=0 MCRB0 == CRB: %s', pw(abs(mcrb(i0) - crb(i0)) < 1e-8*crb(i0))); else, fprintf('     (no p = 0 point)'); end
            if any(pm), fprintf(';  LB / MCRB@eta_bar for p > 0: %s (old: 1.2-3.1x)  %s\n', mat2str([min(lb(pm)./mt(pm)) max(lb(pm)./mt(pm))], 3), pw(all(lb(pm) >= mt(pm)*(1 - 1e-9)))); else, fprintf('\n'); end
            if isfinite(res.pilot_tr_tau) && nP > 1
                for sdb = [0 20]
                    N0 = 10^(-sdb/10);
                    Ss = mcrb_stats(res.G_tau(:, :, iw, im), res.F_tau(:, :, iw, im), res.BB_tau(:, :, iw, im), res.Gc_tau(:, :, iw, im), res.valid(:, :, iw, im), N0, 'median');
                    lbs = Ss.lb;  pil = N0*res.pilot_tr_tau;
                    k = find(lbs > pil, 1);
                    if isempty(k), fprintf('     %2d dB: data-aided LB stays below the pilot-only CRB for all p <= %g\n', sdb, p(end));
                    elseif k == 1, fprintf('     %2d dB: data-aided LB above the pilot-only CRB already at p = %g\n', sdb, p(1));
                    else
                        pstar = p(k-1) + (p(k) - p(k-1))*(pil - lbs(k-1))/(lbs(k) - lbs(k-1));
                        fprintf('     %2d dB: data-vs-pilot crossover at p = %.4f (old: ~1 %% at 20 dB)\n', sdb, pstar);
                    end
                end
            end
        end
        % waveform ordering at every p (median composite MCRB0 at 0 dB)
        m = zeros(nW, nP);
        for iw = 1:nW
            S0 = mcrb_stats(res.G_tau(:, :, iw, im), res.F_tau(:, :, iw, im), res.BB_tau(:, :, iw, im), res.Gc_tau(:, :, iw, im), res.valid(:, :, iw, im), 1, 'median');
            m(iw, :) = S0.mcrb;
        end
        okall = true;  s = '';
        for iw = 1:nW-1
            ok = all(m(iw, :) <= m(iw+1, :)*(1 + 1e-9));  okall = okall && ok;
            s = [s sprintf(' %s/%s %s', res.waveforms{iw}, res.waveforms{iw+1}, mat2str(m(iw, :)./m(iw+1, :), 3))]; %#ok<AGROW>
        end
        fprintf('  %s: MCRB0 ordering %s at every p: %s  (ratios%s)\n', res.mods{im}, strjoin(res.waveforms, ' < '), pw(okall), s);
    end
    % ---- p = 0, L = 40 cross-check against the stored F1.1 / F1.3 genie CRBs --------
    if isempty(res.L40_p0), fprintf('  (no p = 0 point and L ~= 40: F1.1/F1.3 cross-check skipped)\n');  return; end
    G40 = res.L40_p0.Gc_tau;  K40 = res.L40_p0.kappa;
    if strcmpi(res.framing, 'MT40-1D')
        f = fullfile(root, 'results', 'crb_mt40_block.mat');
        if exist(f, 'file')
            b = load(f);
            fprintf('  p = 0, L = 40 vs F1.1 (results/crb_mt40_block.mat, %d draws, mean over valid) [samples^2 at sigma2 = 1]:\n', b.nT);
            for im = 1:nMod
                jm = find(strcmpi(b.mods, res.mods{im}));
                for iw = 1:nW
                    jw = find(strcmpi(b.waveforms, res.waveforms{iw}));
                    if isempty(jm) || isempty(jw), continue; end
                    ref = mean(b.tr_samp2(b.valid(:, jw, jm), jw, jm));
                    here = mean(G40(:, iw, im))/Ts^2;
                    fprintf('     %-6s %-5s here %.4e  F1.1 %.4e  ratio %.3f  %s\n', res.mods{im}, res.waveforms{iw}, here, ref, here/ref, pw(abs(here/ref - 1) < 0.10));
                end
            end
            fprintf('     Jensen (L = 40) here %.4e vs F1.1 %.4e  %s\n', res.L40_p0.jensen_tr_tau/Ts^2, b.jensen_tr_samp2, pw(abs(res.L40_p0.jensen_tr_tau/Ts^2/b.jensen_tr_samp2 - 1) < 1e-8));
        else
            fprintf('  (F1.1 file results/crb_mt40_block.mat not found: cross-check skipped)\n');
        end
    else
        d = dir(fullfile(root, 'results', 'crb_mt40_frame_*.mat'));
        if ~isempty(d)
            f = fullfile(d(1).folder, d(1).name);  b = load(f);
            fprintf('  p = 0, L = 40 vs F1.3 (%s, %d draws, mean) [samples^2 at sigma2 = 1]:\n', d(1).name, b.nT);
            for im = 1:nMod
                key = sprintf('cont_%s', res.mods{im});
                if ~isfield(b, key), continue; end
                for iw = 1:nW
                    jw = find(strcmp(b.(key).waveforms, res.waveforms{iw}));
                    if isempty(jw), continue; end
                    ref = b.(key).tr_tau_samp2(jw);  here = mean(G40(:, iw, im))/Ts^2;
                    fprintf('     %-6s %-8s here %.4e  F1.3 %.4e  ratio %.3f  %s\n', res.mods{im}, res.waveforms{iw}, here, ref, here/ref, pw(abs(here/ref - 1) < 0.10));
                end
            end
        else
            fprintf('  (F1.3 file results/crb_mt40_frame_*.mat not found: cross-check skipped)\n');
        end
    end
    fprintf('  L = 40, p = 0: kappa median %.1e max %.1e over draws\n', median(K40(:)), max(K40(:)));
end

function s = pw(ok)
    if ok, s = 'PASS'; else, s = 'WARN'; end
end

function nT = tier_trials(tier, tab)
    if isnumeric(tier), nT = tier; return; end
    switch lower(tier)
        case 'smoke', nT = tab(1);
        case 'dev',   nT = tab(2);
        case 'paper', nT = tab(3);
        otherwise, error('tier must be smoke|dev|paper or a number');
    end
end
