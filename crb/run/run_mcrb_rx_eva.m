function res = run_mcrb_rx_eva(varargin)
% RUN_MCRB_RX_EVA  MCRB / LB vs SNR from the REAL detector's errors on the EVA
% channel (plan F3.2).  Per trial and SNR: ber_frame (OFDM dftdd CE + one-tap
% MMSE; OTFS embedded-pilot CE + MP) -> (X, Xhat) on the actually transmitted
% grid -> pseudo-true -> A, B at eta0 -> G, F, BB; genie CRB on the same grid;
% pilot-only CRB; and a randomized fixed-count symbol-error control arm.
% Its default uniform mode preserves SER, then selects positions without
% replacement and uniformly different symbols; it does not preserve BER or
% error magnitudes. Per-trial traces are kept; the headline statistic is the median
% for Rayleigh gains (the trial mean of J^-1 has no finite expectation) and
% the mean for constant-magnitude gains (the legacy reference).  Checkpoints
% the result struct after every (waveform, SNR).
%
%   res = run_mcrb_rx_eva('tier', 'smoke', 'reduced', true, 'snr_db', 0:5:25, ...)
% Options: tier smoke(20)|dev(200)|paper(2000) or a number; reduced (true
% default for speed; false = full 9-path EVA); snr_db (0:5:25 smoke, 0:1:25
% otherwise); waveforms {'cp-ofdm','cp-otfs'}; use_ici true; do_bsc true;
% bsc_mode 'uniform' | 'nearest' (bsc_match); seed 9100 (same as run_ber_eva,
% substream keyed on the SNR VALUE so the two drivers share draws at equal
% SNR); otfs_snr_cap 25; Tcp_samples 20; gain_model 'rayleigh'; E 10; mp_iter
% 20; mp_damp 0.7; mp_mode 'edge'; ce_mode 'dftdd'; n_dd 2; pilot_boost_dB 0;
% pilot_frac 0.20; kappa_thresh Inf; statistic 'auto'; tag ''; legacy false
% (true = the old/eva mcrb_ce_ici configuration: EVA-reduced, constmag gains,
% Tcp = 0, dftdd, posterior MP 0.3 x 30, seed 9000, mean statistic, and the
% regression print against old/eva/results/mcrb_ce_ici_results_2000.mat).
    o = struct('tier', 'smoke', 'reduced', true, 'snr_db', [], 'waveforms', {{'cp-ofdm', 'cp-otfs'}}, ...
        'use_ici', true, 'do_bsc', true, 'bsc_mode', 'uniform', 'seed', 9100, 'otfs_snr_cap', 25, 'Tcp_samples', 20, ...
        'gain_model', 'rayleigh', 'E', 10, 'mp_iter', 20, 'mp_damp', 0.7, 'mp_mode', 'edge', 'ce_mode', 'dftdd', 'n_dd', 2, ...
        'pilot_boost_dB', 0, 'pilot_frac', 0.20, 'kappa_thresh', Inf, 'statistic', 'auto', 'tag', '', 'legacy', false, ...
        'doppler_scale', 1, 'v_kmh', 120, 'merge_colliding', false, 'modtype', '4qam', ...
        'snap', true, 'taus_ns', [], 'pows_db', [], 'chan_name', 'EVA-custom', 'variant', '', ...
        'N', 512, 'M', 128, 'deltaf', 15e3);
    o = parse_opts(o, varargin, 'run_mcrb_rx_eva');
    if o.legacy
        o.reduced = true;  o.gain_model = 'constmag';  o.Tcp_samples = 0;  o.ce_mode = 'dftdd';
        o.mp_mode = 'posterior';  o.mp_iter = 30;  o.mp_damp = 0.3;  o.seed = 9000;  o.statistic = 'mean';
        if isempty(o.tag), o.tag = '_legacy'; end
    end
    root = fileparts(fileparts(mfilename('fullpath')));
    nT = tier_trials(o.tier, [20 200 2000]);
    if isempty(o.snr_db), if strcmpi(o.tier, 'smoke'), o.snr_db = 0:5:25; else, o.snr_db = 0:1:25; end, end
    if isnumeric(o.tier), o.tier = sprintf('n%d', nT); end
    if strcmpi(o.statistic, 'auto'), if strcmpi(o.gain_model, 'rayleigh'), o.statistic = 'median'; else, o.statistic = 'mean'; end, end
    N = o.N; M = o.M; df = o.deltaf;  Ts = 1/(N*df);
    g = make_grid(N, M, df, o.Tcp_samples*Ts, 'cp');
    ch = chan_eva(g, 'reduced', o.reduced, 'gain_model', o.gain_model, 'doppler_scale', o.doppler_scale, ...
        'v_kmh', o.v_kmh, 'merge_colliding', o.merge_colliding, 'snap', o.snap, ...
        'taus_ns', o.taus_ns, 'pows_db', o.pows_db, 'name', o.chan_name, 'variant', o.variant);
    mdm = modem(o.modtype);
    nS = numel(o.snr_db);  nW = numel(o.waveforms);
    trn = {'tr_tau_mcrb', 'tr_tau_bb', 'tr_tau_lb', 'tr_tau_crb', 'tr_nu_mcrb', 'tr_nu_bb', 'tr_nu_lb', 'tr_nu_crb', ...
           'tr_tau_lb_bsc', 'tr_nu_lb_bsc', 'tr_tau_pil', 'tr_nu_pil'};
    nQ = numel(trn);
    res = struct('schema_version', 'mcrb_rx_v2_arm_validity');
    res.trial = struct('tr', nan(nT, nQ, nS, nW), 'kappa', nan(nT, nS, nW), 'kappa_bsc', nan(nT, nS, nW), ...
        'exitflag', nan(nT, nS, nW), 'exitflag_bsc', nan(nT, nS, nW), 'nerr', nan(nT, nS, nW), 'nbits', nan(nT, nS, nW), ...
        'ser', nan(nT, nS, nW), 'pcvalid', nan(nT, nS, nW), 'pcvalid_bsc', nan(nT, nS, nW), 'F_over_G', nan(nT, nS, nW));
    res.opts = o;  res.snr_db = o.snr_db;  res.waveforms = o.waveforms;  res.nT = nT;  res.g = g;  res.ch = ch;  res.trn = trn;
    fprintf('MCRB from the real detector: %s, P = %d, %d trials, SNR %s, waveforms %s, MP %s %dx%.1f, CE %s, headline %s\n', ...
        ch.name, ch.P, nT, mat2str(o.snr_db), strjoin(o.waveforms, ' '), o.mp_mode, o.mp_iter, o.mp_damp, o.ce_mode, o.statistic);
    tAll = tic;
    for iw = 1:nW
        wf = o.waveforms{iw};
        for si = 1:nS
            snr = o.snr_db(si);  N0 = 10^(-snr/10);  key = 1000 + round(10*snr);   % same key as run_ber_eva
            if strcmp(wf, 'cp-otfs') && snr > o.otfs_snr_cap, continue; end
            [rx, dom] = rx_for(wf, o);
            acc = nan(nT, nQ);  kap = nan(nT, 2);  ef = nan(nT, 2);  bs = nan(nT, 5);  fg = nan(nT, 1);
            parfor t = 1:nT
                warning('off', 'MATLAB:nearlySingularMatrix');
                rs = crn_stream(o.seed, key, t);
                f = ber_frame(wf, g, ch, N0, rx, rs);
                cr = f.cr;  X = f.X;  Xh = f.Xhat;
                [pc, info] = mcrb_from_grids(wf, X, Xh, g, cr, o.use_ici, dom);
                lb = mcrb_eval(pc, N0);
                Sp = crb_pilot_only(g, cr, o.use_ici, o.pilot_frac, N0, 'cp-ofdm');
                row = [lb.tr_tau_mcrb, lb.tr_tau_bb, lb.tr_tau_lb, lb.tr_tau_crb, ...
                       lb.tr_nu_mcrb, lb.tr_nu_bb, lb.tr_nu_lb, lb.tr_nu_crb, NaN, NaN, Sp.tr_tau, Sp.tr_nu];
                k2 = [pc.kappa NaN];  e2 = [info.exitflag NaN];  vb = NaN;
                if o.do_bsc
                    Xb = bsc_match(X, Xh, f.masks.data, mdm, rs, o.bsc_mode);
                    [pcb, infob] = mcrb_from_grids(wf, X, Xb, g, cr, o.use_ici, dom);
                    lbb = mcrb_eval(pcb, N0);
                    row(9) = lbb.tr_tau_lb;  row(10) = lbb.tr_nu_lb;  k2(2) = pcb.kappa;  e2(2) = infob.exitflag;  vb = double(pcb.valid);
                end
                idx = eta_index(cr.P);
                acc(t, :) = row;  kap(t, :) = k2;  ef(t, :) = e2;
                bs(t, :) = [f.nerr, f.nbits, f.ser, double(pc.valid), vb];
                fg(t) = sum(real(diag(pc.F(idx.tau, idx.tau)))) / max(N0*sum(real(diag(pc.G(idx.tau, idx.tau)))), eps);
            end
            res.trial.tr(:, :, si, iw) = acc;  res.trial.kappa(:, si, iw) = kap(:, 1);  res.trial.kappa_bsc(:, si, iw) = kap(:, 2);
            res.trial.exitflag(:, si, iw) = ef(:, 1);  res.trial.exitflag_bsc(:, si, iw) = ef(:, 2);
            res.trial.nerr(:, si, iw) = bs(:, 1);  res.trial.nbits(:, si, iw) = bs(:, 2);  res.trial.ser(:, si, iw) = bs(:, 3);
            res.trial.pcvalid(:, si, iw) = bs(:, 4);  res.trial.pcvalid_bsc(:, si, iw) = bs(:, 5);  res.trial.F_over_G(:, si, iw) = fg;
            res = reduce_mcrb_rx_stats(res, o);
            fprintf('  %-8s %2d dB | BER %.2e SER %.2e | delay RMSE [ns] (%s): LB %.3f  MCRB0 %.3f  CRB %.3f  BSC-LB %.3f  pilot %.3f | Doppler [Hz]: LB %.3f CRB %.3f | drop %.0f%% kappa<= %.1e (%.0f s)\n', ...
                wf, snr, res.ber(si, iw), res.ser(si, iw), o.statistic, 1e9*res.rmse_tau_lb(si, iw), 1e9*res.rmse_tau_mcrb(si, iw), 1e9*res.rmse_tau_crb(si, iw), ...
                1e9*res.rmse_tau_lb_bsc(si, iw), 1e9*res.rmse_tau_pil(si, iw), res.rmse_nu_lb(si, iw), res.rmse_nu_crb(si, iw), 100*res.drop_frac(si, iw), res.kappa_max(si, iw), toc(tAll));
            if ~exist(fullfile(root, 'results'), 'dir'), mkdir(fullfile(root, 'results')); end
            save(fullfile(root, 'results', sprintf('CKPT_mcrb_rx_%s_%s%s.mat', lower(strrep(ch.name, '-', '_')), o.tier, o.tag)), '-struct', 'res');
        end
    end
    res.seconds = toc(tAll);
    print_acceptance(res);
    save(fullfile(root, 'results', sprintf('mcrb_rx_%s_%s%s.mat', lower(strrep(ch.name, '-', '_')), o.tier, o.tag)), '-struct', 'res');
    plot_mcrb_vs_snr(res);
    plot_ber_mcrb_panels(res);
end

function print_acceptance(res)
    o = res.opts;  snr = res.snr_db(:);  nW = numel(res.waveforms);
    fprintf('\nAcceptance (F3.2), headline statistic: %s\n', o.statistic);
    io = find(strcmp(res.waveforms, 'cp-ofdm'));  it = find(strcmp(res.waveforms, 'cp-otfs'));
    nerr = reshape(sum(res.trial.nerr, 1), numel(snr), nW);
    well = nerr >= 100;                                              % SNRs where the error rate is measured
    if ~isempty(io) && ~isempty(it)
        m = isfinite(res.rmse_tau_crb(:, io)) & isfinite(res.rmse_tau_crb(:, it));
        d = max(abs(res.rmse_tau_crb(m, io)./res.rmse_tau_crb(m, it) - 1));
        fprintf('  genie delay CRB OFDM vs OTFS (same draws): max rel dev %.3f (expect < 0.01 old; < 0.03 with pilot frames)  %s\n', d, pw(d < 0.03));
    end
    for iw = 1:nW
        m = snr > 18 & isfinite(res.rmse_tau_mcrb(:, iw));
        if any(m)
            d = max(abs(res.rmse_tau_mcrb(m, iw)./res.rmse_tau_crb(m, iw) - 1));
            fprintf('  %-8s MCRB0 vs CRB above 18 dB (delay): max rel dev %.3f (expect < 0.02)  %s\n', res.waveforms{iw}, d, pw(d < 0.02));
        end
        k = find(isfinite(res.rmse_tau_lb(:, iw)), 1, 'last');
        if ~isempty(k), fprintf('  %-8s LB / CRB at %g dB (SER %.1e): delay %.3f, Doppler %.3f (LB -> CRB as SER -> 0)\n', res.waveforms{iw}, snr(k), res.ser(k, iw), ...
                res.rmse_tau_lb(k, iw)/res.rmse_tau_crb(k, iw), res.rmse_nu_lb(k, iw)/res.rmse_nu_crb(k, iw)); end
        if o.do_bsc
            ps = res.(['stat_' o.statistic '_paired']);
            r = ps.rmse_tau_lb(:, iw) ./ ps.rmse_tau_lb_bsc(:, iw);
            fprintf('  %-8s real / randomized-control delay LB ratio on jointly valid trials: %s (paired n %s)\n', res.waveforms{iw}, mat2str(r.', 2), mat2str(res.n_valid_paired(:, iw).'));
            if o.legacy && iw == io
                lo = snr <= 8 & isfinite(r) & well(:, iw);  hi = snr >= 13 & isfinite(r) & well(:, iw);
                if any(lo), fprintf('     legacy ref 0-8 dB: 2.5-3.7x  -> here %s  %s\n', mat2str(r(lo).', 2), pw(all(r(lo) > 1.8 & r(lo) < 5))); end
                if any(hi), fprintf('     legacy ref 13-25 dB: 0.6-0.8x -> here %s  %s\n', mat2str(r(hi).', 2), pw(all(r(hi) > 0.45 & r(hi) < 1.05))); end
            end
        end
        fprintf('  %-8s max F/(N0 G) over trials (eta0-arm solver witness): %.1e; converged %.1f%%; drop %s\n', res.waveforms{iw}, ...
            max(res.trial.F_over_G(:, :, iw), [], 'all'), 100*mean(res.trial.exitflag(:, :, iw) > 0, 'all', 'omitnan'), mat2str(res.drop_frac(:, iw).', 2));
    end
    if o.legacy
        % old/eva mcrb_ce_ici (reduced EVA, constmag, Tcp = 0, dftdd, posterior MP, 2000 trials, trial-MEAN):
        % rmse_tau_mcrb there is the full LB (N0 G + F + BB); rmse_tau_crb the genie CRB.
        ref.snr = [0 10 12 15 20 25];
        ref.ofdm_lb  = [5.825e-09 4.513e-10 2.8e-10 1.369e-10 6.698e-11 3.634e-11];
        ref.otfs_lb  = [3.519e-09 2.354e-10 1.64e-10 1.145e-10 6.48e-11 3.598e-11];
        ref.ofdm_crb = [6.445e-10 2.043e-10 1.623e-10 1.145e-10 6.48e-11 3.597e-11];
        ref.ofdm_ber = [0.229 0.03465 0.0211 0.01007 0.003308 0.001326];
        ref.otfs_ber = [0.2116 0.003747 0.000282 1.419e-06 2.321e-08 0];
        [~, ia, ib] = intersect(snr, ref.snr);
        fprintf('  Legacy regression vs old/eva mcrb_ce_ici_results_2000.mat (trial mean; 25 %% on LB, 5 %% on CRB):\n');
        for k = 1:numel(ia)
            s = '';
            if ~isempty(io), s = [s sprintf('  OFDM LB %.3f vs %.3f ns (%s), CRB %.3f vs %.3f (%s), BER %.2e vs %.2e', 1e9*res.stat_mean.rmse_tau_lb(ia(k), io), 1e9*ref.ofdm_lb(ib(k)), ...
                    pw(abs(res.stat_mean.rmse_tau_lb(ia(k), io)/ref.ofdm_lb(ib(k)) - 1) < 0.25), 1e9*res.stat_mean.rmse_tau_crb(ia(k), io), 1e9*ref.ofdm_crb(ib(k)), ...
                    pw(abs(res.stat_mean.rmse_tau_crb(ia(k), io)/ref.ofdm_crb(ib(k)) - 1) < 0.05), res.ber(ia(k), io), ref.ofdm_ber(ib(k)))]; end %#ok<AGROW>
            if ~isempty(it) && isfinite(res.stat_mean.rmse_tau_lb(ia(k), it)), s = [s sprintf(' | OTFS LB %.3f vs %.3f ns (%s), BER %.2e vs %.2e', 1e9*res.stat_mean.rmse_tau_lb(ia(k), it), 1e9*ref.otfs_lb(ib(k)), ...
                    pw(abs(res.stat_mean.rmse_tau_lb(ia(k), it)/ref.otfs_lb(ib(k)) - 1) < 0.25), res.ber(ia(k), it), ref.otfs_ber(ib(k)))]; end %#ok<AGROW>
            fprintf('   %2d dB:%s\n', snr(ia(k)), s);
        end
    end
end

function [rx, dom] = rx_for(wf, o)
    rx = struct('E', o.E, 'mp_iter', o.mp_iter, 'mp_damp', o.mp_damp, 'mp_mode', o.mp_mode, 'n_dd', o.n_dd, ...
        'pilot_boost_dB', o.pilot_boost_dB, 'pilot', 'matched', 'use_ici', o.use_ici, 'ce_mode', o.ce_mode, ...
        'modtype', o.modtype);
    if strcmp(wf, 'cp-ofdm'), rx.ofdm_csi = 'ce';  dom = 'tf'; else, rx.otfs_csi = 'ce';  dom = 'dd'; end
end

function [pc, info] = mcrb_from_grids(wf, X, Xh, g, cr, use_ici, dom)
    [eta0, info] = pseudo_true(X, Xh, g, cr, use_ici, struct('domain', dom));
    [a0, t0, n0] = eta_unpack(eta0);
    ce = cr;  ce.alpha = a0;  ce.tau = t0;  ce.nu = n0;
    if strcmp(wf, 'cp-ofdm')
        [A, ~, Jt, b] = mcrb_AB_cp_ofdm(X, Xh, g, ce, use_ici, 1, cr);
    else
        [A, ~, Jt, b] = mcrb_AB_cp_otfs(X, Xh, g, ce, use_ici, 1, cr);
    end
    Jg = fim(X, wf, g, cr, use_ici, 1);           % genie CRB on the actual transmitted grid
    pc = mcrb_assemble(A, Jt, b, Jg, cr.eta_bar, eta0, g);
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
