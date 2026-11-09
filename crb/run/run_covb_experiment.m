function res = run_covb_experiment(nT)
% RUN_COVB_EXPERIMENT  Settle the Cov(b) question numerically (plan F3.0).
% latex/mcrb/cp_otfs_mcrb.tex Sec. 3-4 derives the CP-OTFS bias covariance as
% Kronecker (delay Gram x Doppler Gram) traces and says it does not reduce to
% the OFDM one; latex/mcrb/mcrb_expectation_comparison.tex Sec. 7.3 argues the
% no-ICI collapse Cov(b) = 4p(1-p)|alpha|^2 J/sigma2 extends to CP-OTFS.
% Monte-Carlo Cov(b) on the legacy 16 x 8 preset, p = 0.1, no ICI, for
%   case 1: single target (P = 1)  -> MC vs closed forms vs the collapse
%   case 2: three targets (P = 3)  -> MC vs closed forms, OFDM vs OTFS gap
% with batch-split standard errors.  Saves results/covb_experiment.mat and
% figures/fig_covb_experiment.png; the outcome is recorded in README.md.
% Runtime code never depends on it (it always uses realization B).
    if nargin < 1, nT = 40000; end
    root = fileparts(fileparts(mfilename('fullpath')));
    N = 16; M = 8; df = 15e3;  sigma2 = 1.3;  p = 0.1;
    g0 = make_grid(N, M, df, 0, 'cp');
    ch = chan_eva(g0, 'legacy_preset', true);  cr3 = chan_realize(ch, g0, []);
    cr1 = struct('P', 1, 'alpha', cr3.alpha(2), 'tau', cr3.tau(2), 'nu', cr3.nu(2));
    cr1.eta_bar = eta_pack(cr1.alpha, cr1.tau, cr1.nu);
    cases = {cr1, cr3};
    mdm = modem('4qam');
    nB = 10;  nb = floor(nT/nB);
    assert(nb >= 1, 'run_covb_experiment:nT', 'nT must be >= %d', nB);
    nT = nb*nB;                                   % trials actually averaged
    wfs = {'cp-ofdm', 'cp-otfs'};
    res = struct('p', p, 'nT', nT, 'sigma2', sigma2, 'nB', nB);
    fprintf('Cov(b) experiment: legacy 16x8, p = %.2f, no ICI, %d trials (%d batches)\n', p, nT, nB);
    fig = new_fig([100 100 1100 380]);
    for ic = 1:2
        cr = cases{ic};  nq = 4*cr.P;
        Sb = zeros(nq, nq, nB, 2);  sb = zeros(nq, nB, 2);  batches = zeros(nq, nq, nB, 2);
        parfor ib = 1:nB
            rs = crn_stream(96 + ic, ib, 1);
            Sbb = zeros(nq, nq, 2);  sbb = zeros(nq, 2);
            for t = 1:nb
                x = gen_symbols(N, M, '4qam', rs);  xh = bsc_symbols(x, p, mdm, rs);
                [~, ~, ~, b1] = mcrb_AB_cp_ofdm(x, xh, g0, cr, false, sigma2);
                [~, ~, ~, b2] = mcrb_AB_cp_otfs(x, xh, g0, cr, false, sigma2);
                Sbb(:, :, 1) = Sbb(:, :, 1) + b1*b1.';  sbb(:, 1) = sbb(:, 1) + b1;
                Sbb(:, :, 2) = Sbb(:, :, 2) + b2*b2.';  sbb(:, 2) = sbb(:, 2) + b2;
            end
            bl = zeros(nq, nq, 2);
            for w = 1:2, bl(:, :, w) = Sbb(:, :, w)/nb - (sbb(:, w)/nb)*(sbb(:, w)/nb).'; end
            batches(:, :, ib, :) = reshape(bl, nq, nq, 1, 2);
            Sb(:, :, ib, :) = reshape(Sbb, nq, nq, 1, 2);  sb(:, ib, :) = reshape(sbb, nq, 1, 2);
        end
        c = struct('P', cr.P);
        for w = 1:2
            wn = strrep(wfs{w}, '-', '_');
            St = sum(Sb(:, :, :, w), 3);  st = sum(sb(:, :, w), 2);
            c.(['mc_' wn]) = St/nT - (st/nT)*(st/nT).';
            c.(['batch_' wn]) = batches(:, :, :, w);
        end
        c.cf_cp_ofdm = covb_cp_ofdm(g0, cr, false, sigma2, p);
        c.cf_cp_otfs = covb_cp_otfs(g0, cr, false, sigma2, p);
        J = fim_expected('cp-ofdm', g0, cr, false, sigma2);
        c.J = J;
        idx = eta_index(cr.P);
        blocks = {'gain', idx.gain; 'kin', idx.kin};
        fprintf('\n-- case %d: P = %d --\n', ic, cr.P);
        fprintf('%-8s %-6s %-16s %-16s %-16s\n', 'wf', 'block', 'MC vs cf', 'MC OFDM vs OTFS', 'rel SE');
        for w = 1:2
            wn = strrep(wfs{w}, '-', '_');
            for k = 1:2
                ii = blocks{k, 2};
                mc = c.(['mc_' wn])(ii, ii);  cf = c.(['cf_' wn])(ii, ii);
                other = c.(['mc_' strrep(wfs{3-w}, '-', '_')])(ii, ii);
                bt = reshape(c.(['batch_' wn])(ii, ii, :), [], nB);
                se = norm(std(bt, 0, 2))/sqrt(nB)/max(norm(mc(:)), eps);
                st_ = struct('mc_vs_cf', relerr(mc, cf), 'mc_vs_other', relerr(mc, other), 'rel_se', se);
                c.([wn '_' blocks{k, 1}]) = st_;
                fprintf('%-8s %-6s %-16.3e %-16.3e %-16.3e\n', wn, blocks{k, 1}, st_.mc_vs_cf, st_.mc_vs_other, se);
            end
        end
        c.cf_gap_kin = relerr(c.cf_cp_otfs(idx.kin, idx.kin), c.cf_cp_ofdm(idx.kin, idx.kin));
        c.cf_gap_all = relerr(c.cf_cp_otfs, c.cf_cp_ofdm);
        fprintf('closed-form OTFS vs OFDM Cov(b): kinematic block %.3e, full %.3e\n', c.cf_gap_kin, c.cf_gap_all);
        % PAIRED difference (same (x, xhat) draws for both waveforms -> common random
        % numbers): far sharper than comparing the two MC estimates separately.
        % Scalar test statistic: the kinematic-block TRACE of the per-batch paired
        % difference, t = mean/SE against a t distribution with nB-1 dof; plus the
        % per-entry z-scores of the unique kinematic entries.
        ii = idx.kin;  nk = numel(ii);
        dB = c.batch_cp_otfs(ii, ii, :) - c.batch_cp_ofdm(ii, ii, :);
        dmc = mean(dB, 3);  dse_ent = std(dB, 0, 3)/sqrt(nB);
        dcf = c.cf_cp_otfs(ii, ii) - c.cf_cp_ofdm(ii, ii);
        dtr = arrayfun(@(ib) trace(dB(:, :, ib)), 1:nB);
        t_tr = mean(dtr) / max(std(dtr)/sqrt(nB), eps);
        t_crit = tinv_local(0.995, nB - 1);
        up = triu(true(nk));
        z = dmc(up) ./ max(dse_ent(up), eps);
        c.paired = struct('dmc', dmc, 'dcf', dcf, 'dse_entry', dse_ent, 't_trace', t_tr, 't_crit', t_crit, ...
            'trace_mc', mean(dtr), 'trace_cf', trace(dcf), 'trace_se', std(dtr)/sqrt(nB), 'z_max', max(abs(z)), 'chi2', sum(z.^2), 'n_entries', nnz(up), ...
            'rel_err', norm(dmc(:) - dcf(:))/max(norm(dcf(:)), eps), 'rel_se', norm(dse_ent(:))/max(norm(dcf(:)), eps), ...
            'gap_rel_to_covb', norm(dcf(:))/norm(reshape(c.cf_cp_ofdm(ii, ii), [], 1)));
        fprintf(['paired MC (OTFS - OFDM) kinematic block: trace t = %.1f (crit %.2f at 99%%, %d dof), MC trace %.3e +/- %.1e vs closed-form %.3e; ' ...
                 'max |z| over %d entries %.1f, chi2 %.1f; full-block rel err vs closed form %.3e (rel SE %.3e); closed-form gap = %.3e of Cov(b)\n'], ...
            c.paired.t_trace, c.paired.t_crit, nB - 1, c.paired.trace_mc, c.paired.trace_se, c.paired.trace_cf, c.paired.n_entries, c.paired.z_max, c.paired.chi2, ...
            c.paired.rel_err, c.paired.rel_se, c.paired.gap_rel_to_covb);
        if cr.P == 1
            collapse = 4*p*(1-p)*abs(cr.alpha)^2/sigma2 * J;
            c.collapse = collapse;
            c.cf_ofdm_vs_collapse = relerr(c.cf_cp_ofdm, collapse);
            c.cf_otfs_vs_collapse = relerr(c.cf_cp_otfs, collapse);
            c.mc_otfs_vs_collapse_kin = relerr(c.mc_cp_otfs(idx.kin, idx.kin), collapse(idx.kin, idx.kin));
            fprintf('P = 1 collapse 4p(1-p)|a|^2 J/sigma2: closed-form OFDM %.3e, closed-form OTFS %.3e, MC OTFS kin %.3e (SE %.3e)\n', ...
                c.cf_ofdm_vs_collapse, c.cf_otfs_vs_collapse, c.mc_otfs_vs_collapse_kin, c.cp_otfs_kin.rel_se);
        end
        res.(sprintf('P%d', cr.P)) = c;
        subplot(1, 2, ic);
        if cr.P == 1, vals = [c.mc_cp_ofdm(:), c.mc_cp_otfs(:), c.cf_cp_otfs(:), c.collapse(:)];  lg = {'MC OFDM', 'MC OTFS', 'cf OTFS', 'collapse'};
        else,         vals = [c.mc_cp_ofdm(:), c.mc_cp_otfs(:), c.cf_cp_ofdm(:), c.cf_cp_otfs(:)];  lg = {'MC OFDM', 'MC OTFS', 'cf OFDM', 'cf OTFS'}; end
        bar(vals ./ max(abs(vals(:))));  grid on;  legend(lg, 'Location', 'best');
        title(sprintf('P = %d: Cov(b) entries (column-major %dx%d), p = %.2f, no ICI', cr.P, nq, nq, p));
    end
    % verdicts
    c1 = res.P1;  c3 = res.P3;
    tol1 = 5*c1.cp_otfs_kin.rel_se + 0.01;
    if c1.cf_otfs_vs_collapse < 1e-10 && c1.mc_otfs_vs_collapse_kin < tol1
        res.verdict_P1 = 'single target, no ICI: the cp_otfs_mcrb.tex Sec. 4 Kronecker form evaluates EXACTLY to the collapse 4p(1-p)|a|^2 J/sigma2 and MC agrees -> no contradiction for P = 1';
    else
        res.verdict_P1 = sprintf('single target: closed-form OTFS vs collapse %.2e, MC vs collapse %.2e (tol %.2e) -> inspect', c1.cf_otfs_vs_collapse, c1.mc_otfs_vs_collapse_kin, tol1);
    end
    tol3 = 5*max(c3.cp_otfs_kin.rel_se, c3.cp_ofdm_kin.rel_se) + 0.01;
    pr = c3.paired;
    if c3.cp_otfs_kin.mc_vs_cf < tol3 && c3.cp_ofdm_kin.mc_vs_cf < tol3
        tr_ok = abs(pr.trace_mc - pr.trace_cf) < 3*pr.trace_se + 0.05*abs(pr.trace_cf);
        if abs(pr.t_trace) > pr.t_crit && tr_ok
            res.verdict_P3 = sprintf(['three targets: both closed forms match MC (OFDM %.2e, OTFS %.2e); the PAIRED MC difference OTFS-OFDM is ' ...
                'non-zero (trace t = %.1f > %.2f) and its trace matches the closed-form difference (%.3e +/- %.1e vs %.3e) -> the OTFS Cov(b) does NOT ' ...
                'reduce to the OFDM one for P > 1 (gap %.2e of Cov(b); cp_otfs_mcrb.tex Sec. 3-4 supported, the Sec. 7.3 collapse holds only at P = 1 / no ICI)'], ...
                c3.cp_ofdm_kin.mc_vs_cf, c3.cp_otfs_kin.mc_vs_cf, pr.t_trace, pr.t_crit, pr.trace_mc, pr.trace_se, pr.trace_cf, pr.gap_rel_to_covb);
        elseif abs(pr.t_trace) <= pr.t_crit
            res.verdict_P3 = sprintf('three targets: both closed forms match MC; the paired OTFS-OFDM MC difference is not significant (trace t = %.1f, crit %.2f; closed-form gap %.2e of Cov(b)) -> not resolved at %d trials', pr.t_trace, pr.t_crit, pr.gap_rel_to_covb, nT);
        else
            res.verdict_P3 = sprintf('three targets: paired MC difference is significant (t = %.1f) but its trace %.3e +/- %.1e disagrees with the closed-form difference %.3e -> inspect', pr.t_trace, pr.trace_mc, pr.trace_se, pr.trace_cf);
        end
    else
        res.verdict_P3 = sprintf('three targets: MC vs closed form OFDM %.2e, OTFS %.2e (tol %.2e) -> inconclusive / inspect', c3.cp_ofdm_kin.mc_vs_cf, c3.cp_otfs_kin.mc_vs_cf, tol3);
    end
    fprintf('\nVERDICT P=1: %s\nVERDICT P=3: %s\n', res.verdict_P1, res.verdict_P3);
    if ~exist(fullfile(root, 'results'), 'dir'), mkdir(fullfile(root, 'results')); end
    if ~exist(fullfile(root, 'figures'), 'dir'), mkdir(fullfile(root, 'figures')); end
    save(fullfile(root, 'results', 'covb_experiment.mat'), '-struct', 'res');
    print(fig, fullfile(root, 'figures', 'fig_covb_experiment.png'), '-dpng', '-r120');
end

function q = tinv_local(pr, dof)
% Student-t quantile without the Statistics Toolbox (bisection on the CDF via betainc).
    cdf = @(t) 1 - 0.5*betainc(dof./(dof + t.^2), dof/2, 0.5);
    lo = 0;  hi = 100;
    for it = 1:200
        mid = 0.5*(lo + hi);
        if cdf(mid) < pr, lo = mid; else, hi = mid; end
    end
    q = 0.5*(lo + hi);
end
