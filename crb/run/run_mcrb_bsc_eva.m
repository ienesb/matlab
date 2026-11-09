function res = run_mcrb_bsc_eva(varargin)
% RUN_MCRB_BSC_EVA  MCRB / LB vs BSC error rate p on the EVA channel (plan F3.1)
% for {cp-ofdm, cp-otfs} x {no-ICI, ICI}, all-data grid, Gray-QPSK BSC.
% Per trial and p: pseudo-true eta0 (VARPRO, warm-started along p), A, B at
% eta0, genie J at the truth, and the affine pieces G, F, BB, Gc (sigma2 = 1),
% so every SNR follows by mcrb_eval; also the true-parameter MCRB (eta0 := eta_bar)
% and the full A(eta_bar) / J matrices for the mean-field check.
%
%   res = run_mcrb_bsc_eva('tier', 'smoke', 'reduced', true, ...)
% Options: tier smoke(20)|dev(200)|paper(2000) or a number; reduced (true:
% EVA-reduced, false: full EVA); p_list (0:0.05:0.5 smoke, 0:0.02:0.5 otherwise;
% zoom: 0:0.002:0.1); snr_db [0 10 20]; waveforms; ici_list [false true];
% seed 20260707; Tcp_samples 20; gain_model 'rayleigh'; kappa_thresh Inf;
% statistic 'auto' (median for Rayleigh gains, mean otherwise; both are always
% stored and plotted); tag '' (appended to the file names).
% Statistics are of the PER-TRIAL composite bounds (mcrb_stats), over trials
% with pc.valid, a converged pseudo-true solve and kappa <= kappa_thresh.
    o = struct('tier', 'smoke', 'reduced', true, 'p_list', [], 'snr_db', [0 10 20], ...
        'waveforms', {{'cp-ofdm', 'cp-otfs'}}, 'ici_list', [false true], 'seed', 20260707, ...
        'Tcp_samples', 20, 'gain_model', 'rayleigh', 'modtype', '4qam', 'kappa_thresh', Inf, ...
        'statistic', 'auto', 'tag', '');
    o = parse_opts(o, varargin, 'run_mcrb_bsc_eva');
    root = fileparts(fileparts(mfilename('fullpath')));
    nT = tier_trials(o.tier, [20 200 2000]);
    if isempty(o.p_list), if strcmpi(o.tier, 'smoke'), o.p_list = 0:0.05:0.5; else, o.p_list = 0:0.02:0.5; end, end
    if isnumeric(o.tier), o.tier = sprintf('n%d', nT); end
    if strcmpi(o.statistic, 'auto'), if strcmpi(o.gain_model, 'rayleigh'), o.statistic = 'median'; else, o.statistic = 'mean'; end, end
    N = 512; M = 128; df = 15e3;  Ts = 1/(N*df);
    g = make_grid(N, M, df, o.Tcp_samples*Ts, 'cp');
    ch = chan_eva(g, 'reduced', o.reduced, 'gain_model', o.gain_model);
    P = ch.P;  nP = numel(o.p_list);  nW = numel(o.waveforms);  nI = numel(o.ici_list);
    mdm = modem(o.modtype);
    fprintf('MCRB vs p (BSC): %s, P = %d, %d trials, %d p-points, %d waveforms x %d ICI, headline statistic: %s\n', ch.name, P, nT, nP, nW, nI, o.statistic);
    fld = {'G_tau', 'F_tau', 'BB_tau', 'Gc_tau', 'G_nu', 'F_nu', 'BB_nu', 'Gc_nu', 'Gt_tau', 'Ft_tau', 'Gt_nu', 'Ft_nu', ...
           'kappa', 'pcvalid', 'pt_iters', 'pt_sec', 'pt_exitflag', 'pt_resnorm'};
    wfs = o.waveforms;  ici_list = o.ici_list;  p_list = o.p_list;  seed = o.seed;  modtype = o.modtype;
    t0 = tic;
    out = cell(nT, 1);
    parfor t = 1:nT
        loc = struct();  for k = 1:numel(fld), loc.(fld{k}) = nan(nP, nW, nI); end
        loc.Afull = nan(nP, nW, nI, 4*P, 4*P);  loc.Jfull = nan(nP, nW, nI, 4*P, 4*P);
        rs = crn_stream(seed, 0, t);
        cr = chan_realize(ch, g, rs);
        Sg = gen_symbols(N, M, modtype, rs);
        idx = eta_index(cr.P);
        tr = @(Xm, i) sum(real(diag(Xm(i, i))));
        for ii = 1:nI
            ici = ici_list(ii);
            for iw = 1:nW
                wf = wfs{iw};
                if strcmp(wf, 'cp-ofdm'), dom = 'tf'; else, dom = 'dd'; end
                X = Sg;
                Jg = fim(X, wf, g, cr, ici, 1);
                eta_prev = [];
                for ip = 1:nP
                    p = p_list(ip);
                    Xh = bsc_symbols(X, p, mdm, crn_stream(seed, ip, t));      % flips paired across waveforms / ICI
                    if p == 0
                        eta0 = cr.eta_bar;  pinfo = struct('iterations', 0, 'seconds', 0, 'exitflag', 1, 'resnorm', 0);
                    else
                        [eta0, pinfo] = pseudo_true(X, Xh, g, cr, ici, struct('domain', dom, 'eta_init', eta_prev, 'p', p));
                    end
                    eta_prev = eta0;
                    [a0, t0v, n0] = eta_unpack(eta0);
                    ce = cr;  ce.alpha = a0;  ce.tau = t0v;  ce.nu = n0;
                    if strcmp(wf, 'cp-ofdm')
                        [A, ~, Jt, b]     = mcrb_AB_cp_ofdm(X, Xh, g, ce, ici, 1, cr);
                        [At, ~, Jtt, bt]  = mcrb_AB_cp_ofdm(X, Xh, g, cr, ici, 1, cr);   % eta0 := eta_bar
                    else
                        [A, ~, Jt, b]     = mcrb_AB_cp_otfs(X, Xh, g, ce, ici, 1, cr);
                        [At, ~, Jtt, bt]  = mcrb_AB_cp_otfs(X, Xh, g, cr, ici, 1, cr);
                    end
                    pc  = mcrb_assemble(A, Jt, b, Jg, cr.eta_bar, eta0, g);
                    pct = mcrb_assemble(At, Jtt, bt, Jg, cr.eta_bar, cr.eta_bar, g);
                    loc.G_tau(ip, iw, ii) = tr(pc.G, idx.tau);   loc.F_tau(ip, iw, ii) = tr(pc.F, idx.tau);
                    loc.BB_tau(ip, iw, ii) = tr(pc.BB, idx.tau); loc.Gc_tau(ip, iw, ii) = tr(pc.Gc, idx.tau);
                    loc.G_nu(ip, iw, ii) = tr(pc.G, idx.nu);     loc.F_nu(ip, iw, ii) = tr(pc.F, idx.nu);
                    loc.BB_nu(ip, iw, ii) = tr(pc.BB, idx.nu);   loc.Gc_nu(ip, iw, ii) = tr(pc.Gc, idx.nu);
                    loc.Gt_tau(ip, iw, ii) = tr(pct.G, idx.tau); loc.Ft_tau(ip, iw, ii) = tr(pct.F, idx.tau);
                    loc.Gt_nu(ip, iw, ii) = tr(pct.G, idx.nu);   loc.Ft_nu(ip, iw, ii) = tr(pct.F, idx.nu);
                    loc.kappa(ip, iw, ii) = pc.kappa;            loc.pcvalid(ip, iw, ii) = double(pc.valid);
                    loc.pt_iters(ip, iw, ii) = pinfo.iterations;  loc.pt_sec(ip, iw, ii) = pinfo.seconds;
                    loc.pt_exitflag(ip, iw, ii) = pinfo.exitflag; loc.pt_resnorm(ip, iw, ii) = pinfo.resnorm;
                    loc.Afull(ip, iw, ii, :, :) = At;  loc.Jfull(ip, iw, ii, :, :) = Jg;
                end
            end
        end
        out{t} = loc;
        fprintf('  trial %d/%d (%.0f s)\n', t, nT, toc(t0));
    end
    res = struct();
    for k = 1:numel(fld)
        A4 = nan(nT, nP, nW, nI);
        for t = 1:nT, A4(t, :, :, :) = out{t}.(fld{k}); end
        res.(fld{k}) = A4;
    end
    res.Afull = nan(nT, nP, nW, nI, 4*P, 4*P);  res.Jfull = res.Afull;
    for t = 1:nT, res.Afull(t, :, :, :, :, :) = out{t}.Afull;  res.Jfull(t, :, :, :, :, :) = out{t}.Jfull; end
    res.valid = res.pcvalid > 0 & res.pt_exitflag > 0 & isfinite(res.kappa) & res.kappa <= o.kappa_thresh;
    res.drop_frac = reshape(1 - mean(res.valid, 1), nP, nW, nI);
    res.statistic = sprintf('%s over valid trials of the per-trial composite bounds (valid: pc.valid, converged pseudo-true, kappa <= %g)', o.statistic, o.kappa_thresh);
    res.p_list = o.p_list;  res.snr_db = o.snr_db;  res.opts = o;  res.nT = nT;  res.g = g;  res.ch = ch;  res.P = P;
    res.waveforms = o.waveforms;  res.ici_list = o.ici_list;  res.seconds = toc(t0);
    print_acceptance(res);
    if ~exist(fullfile(root, 'results'), 'dir'), mkdir(fullfile(root, 'results')); end
    save(fullfile(root, 'results', sprintf('mcrb_bsc_%s_%s%s.mat', lower(strrep(ch.name, '-', '_')), o.tier, o.tag)), '-struct', 'res');
    plot_mcrb_vs_p(res, 'mean');
    plot_mcrb_vs_p(res, 'median');
end

function print_acceptance(res)
    o = res.opts;  p = res.p_list;  nP = numel(p);  nW = numel(res.waveforms);  nI = numel(res.ici_list);  P = res.P;
    i0 = find(p == 0, 1);  pm = p > 0;
    stats = {o.statistic, setdiff({'mean', 'median'}, {o.statistic})};  stats = [stats(1) stats{2}];
    fprintf('\nAcceptance (F3.1) -- statistics of the per-trial composite bounds over valid trials; drop fraction max %.1f%%\n', 100*max(res.drop_frac(:)));
    pw3 = p > 0 & p <= 0.3;                                          % mean-field regime (A stays well conditioned)
    if any(pm), fprintf('  pseudo-true: mean %.1f iterations, %.2f s; converged (exitflag > 0) in %.1f%% of solves; max F/(N0 G) at 20 dB in the eta0 arm = %.1e for p <= 0.3 (solver-residual witness, expect << 1; %.1e over all p incl. the p -> 0.5 wall)\n', ...
        mean(res.pt_iters(:, pm, :, :), 'all'), mean(res.pt_sec(:, pm, :, :), 'all'), 100*mean(res.pt_exitflag(:, pm, :, :) > 0, 'all'), ...
        max(res.F_tau(:, pw3, :, :) ./ (1e-2*res.G_tau(:, pw3, :, :)), [], 'all'), max(res.F_tau(:, pm, :, :) ./ (1e-2*res.G_tau(:, pm, :, :)), [], 'all')); end
    for is = 1:numel(stats)
        st = stats{is};
        fprintf(' [%s] delay RMSE bounds at 10 dB [ns]\n', st);
        S10 = mcrb_stats(res.G_tau, res.F_tau, res.BB_tau, res.Gc_tau, res.valid, 1e-1, st);
        T10 = mcrb_stats(res.Gt_tau, res.Ft_tau, 0*res.BB_tau, res.Gc_tau, res.valid, 1e-1, st);
        S20 = mcrb_stats(res.G_tau, res.F_tau, res.BB_tau, res.Gc_tau, res.valid, 1e-2, st);
        rv = @(x) reshape(x, 1, []);
        for iw = 1:nW, for ii = 1:nI
            mcrb = 1e9*rv(S10.rmse_mcrb(1, :, iw, ii));  lb = 1e9*rv(S10.rmse_lb(1, :, iw, ii));
            crb = 1e9*rv(S10.rmse_crb(1, :, iw, ii));    mt = 1e9*rv(T10.rmse_mcrb(1, :, iw, ii));
            fprintf('  %-8s ici=%d  p=%s\n     MCRB0 %s\n     LB    %s\n     MCRB@eta_bar %s\n     CRB   %s   drop %s\n', res.waveforms{iw}, res.ici_list(ii), ...
                mat2str(p, 2), mat2str(mcrb, 3), mat2str(lb, 3), mat2str(mt, 3), mat2str(crb, 3), mat2str(res.drop_frac(:, iw, ii).', 2));
            if ~isempty(i0)
                e0 = abs(mcrb(i0) - crb(i0))/crb(i0) + abs(lb(i0) - crb(i0))/crb(i0);
                fprintf('     p=0: MCRB0 = LB = CRB rel %.1e  %s', e0, pw(e0 < 1e-8));
            else
                fprintf('     (no p = 0 point: identity check skipped)');
            end
            mw = p <= 0.4;  dm = diff(mcrb(mw));
            fprintf(';  MCRB0 monotone in p (p <= 0.4): %s\n', pw(all(dm >= -1e-6*mcrb(1:nnz(mw)-1))));
            bbr = rv(S20.bb(1, :, iw, ii)) ./ rv(S20.mcrb(1, :, iw, ii));  mb = p >= 0.05;
            if any(mb), fprintf('     BB / MCRB0 at 20 dB for p >= 0.05: %s  (expect > 1)  %s\n', mat2str(bbr(mb), 2), pw(all(bbr(mb) > 1))); end
        end, end
        if nW == 2
            for ii = 1:nI
                m1 = rv(S10.mcrb(1, :, 1, ii));  m2 = rv(S10.mcrb(1, :, 2, ii));  mm = p <= 0.3;
                gap = max(abs(sqrt(m1(mm)./m2(mm)) - 1));
                fprintf('  %s vs %s MCRB0 (ici=%d, p <= 0.3): max rel gap %.3f (expect < 0.05)  %s\n', res.waveforms{1}, res.waveforms{2}, res.ici_list(ii), gap, pw(gap < 0.05));
            end
        end
    end
    % mean-field check on the within-target 4x4 blocks of E{A(eta_bar)} vs -(1-2p) E{J}
    % (cross-target blocks of A_Delta are zero, so those obey E{A} = -J with no (1-2p) factor)
    idx = eta_index(P);
    names = {'aR', 'aI', 'tau', 'nu'};
    ent_noici = [1 3; 1 4; 2 3; 2 4; 3 3; 3 4; 4 4];
    ent_ici_assert = [3 3; 3 4; 1 3; 2 3];  ent_ici_report = [1 4; 2 4; 4 4];
    for iw = 1:nW, for ii = 1:nI
        ici = res.ici_list(ii);
        if ici, ent_a = ent_ici_assert;  ent_r = ent_ici_report; else, ent_a = ent_noici;  ent_r = zeros(0, 2); end
        dev_a = nan(nP, size(ent_a, 1));  dev_r = nan(nP, size(ent_r, 1));  dev_x = nan(nP, 1);
        for ip = 1:nP
            if p(ip) == 0 || p(ip) > 0.2, continue; end
            v = res.valid(:, ip, iw, ii);
            if ~any(v), continue; end
            Ab = reshape(mean(res.Afull(v, ip, iw, ii, :, :), 1), 4*P, 4*P);  Jb = reshape(mean(res.Jfull(v, ip, iw, ii, :, :), 1), 4*P, 4*P);
            c = 1 - 2*p(ip);
            dev_a(ip, :) = 0;  dev_r(ip, :) = 0;
            for q = 1:P
                i4 = idx.target(q);  Aq = Ab(i4, i4);  Jq = Jb(i4, i4);
                dq = @(k, l) abs(Aq(k, l) + c*Jq(k, l)) / (c*sqrt(abs(Jq(k, k)*Jq(l, l))));
                for e = 1:size(ent_a, 1), dev_a(ip, e) = max(dev_a(ip, e), dq(ent_a(e,1), ent_a(e,2))); end
                for e = 1:size(ent_r, 1), dev_r(ip, e) = max(dev_r(ip, e), dq(ent_r(e,1), ent_r(e,2))); end
            end
            % cross-target blocks: E{A} = -J (info)
            if P > 1
                mx = 0;
                for q1 = 1:P, for q2 = 1:P
                    if q1 == q2, continue; end
                    i1 = idx.target(q1);  i2 = idx.target(q2);
                    mx = max(mx, norm(Ab(i1, i2) + Jb(i1, i2), 'fro') / max(norm(Jb(i1, i2), 'fro'), eps));
                end, end
                dev_x(ip) = mx;
            end
        end
        mp = ~isnan(dev_a(:, 1));
        if ~any(mp), continue; end
        lab = arrayfun(@(e) sprintf('(%s,%s)', names{ent_a(e,1)}, names{ent_a(e,2)}), 1:size(ent_a,1), 'UniformOutput', false);
        fprintf('  mean-field %s ici=%d (0 < p <= 0.2, within-target entries, worst target): max dev %s  %s\n', res.waveforms{iw}, ici, ...
            strjoin(arrayfun(@(e) sprintf('%s %.3f', lab{e}, max(dev_a(mp, e))), 1:numel(lab), 'UniformOutput', false), ', '), pw(max(dev_a(mp, :), [], 'all') < 0.05));
        if ~isempty(ent_r)
            labr = arrayfun(@(e) sprintf('(%s,%s)', names{ent_r(e,1)}, names{ent_r(e,2)}), 1:size(ent_r,1), 'UniformOutput', false);
            fprintf('     reported only (ICI W0 / R-vs-U terms): %s\n', strjoin(arrayfun(@(e) sprintf('%s %.3f', labr{e}, max(dev_r(mp, e))), 1:numel(labr), 'UniformOutput', false), ', '));
        end
        if P > 1, fprintf('     cross-target blocks E{A} = -J: max rel dev %.2e (info)\n', max(dev_x(mp))); end
    end, end
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
