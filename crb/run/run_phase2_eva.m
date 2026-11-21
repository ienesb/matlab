function res = run_phase2_eva(varargin)
% RUN_PHASE2_EVA  Submission-plan Phase 2 driver: run_mcrb_rx_eva's frame loop
% (identical draws and identical legacy outputs for the same seed, SNR and
% trial) plus what the reviewer-gap experiments need, in one resumable run:
%   * several randomized fixed-count control arms on the SAME error positions
%     ('uniform': random different symbol, the paper's control; 'nearest': one
%     Gray bit flipped, so only the placement differs from the real errors);
%   * the mismatched ML estimator (mmle_cp 'lsq') on the detected reference,
%     end to end (the observation that produced the decisions) and with a fresh
%     noise draw (the conditional setting the bound is derived for);
%   * per-trial vectors and kinematic-block matrices (eta_bar, eta0, G, Gc), so
%     later statistics need no new simulation;
%   * a speed option, a genie-CSI detector option and the full 9-path EVA.
%
%   res = run_phase2_eva('tier', 500, 'v_kmh', 240, 'snr_db', 0:2:24, ...)
% Options (defaults): tier 500 (trials, or smoke|dev|paper); v_kmh 120; reduced
% true (false = full 9-path EVA); snr_db 0:2:24 (evaluated in the given order);
% waveforms {'cp-ofdm','cp-otfs'}; seed 9100; controls {'uniform','nearest'}
% (the first one fills the legacy *_bsc slots); mmle_snr [] (SNR values that
% get the estimator arms); mmle_arms, default all four of e2e, fresh, e2e_truth,
% fresh_truth ('e2e': the
% observation that produced the decisions; 'fresh': the same decisions with an
% independent noise draw; suffix '_truth': joint fit started at the true
% parameter instead of the periodogram/CLEAN start, i.e. the local optimum the
% bound describes); csi 'ce' | 'genie' (detector only: OFDM banded MMSE
% with the true channel, OTFS MP on the true DD response; the sensing side never
% sees it); chunk 40 (trials per checkpoint); resume true; pass [] (SNR values to
% evaluate in THIS call, in that order; default all of snr_db; the result file
% is written once every SNR is finished, so a campaign can take the most useful
% points of every run first); frames [] (trials per cell to reach in THIS call,
% default all; a call with 'frames', 250 runs trials 1-250 of every cell of the
% pass and keeps the checkpoint, a later call without it tops the same run up to
% the full count, and trials 1-250 are the same either way); ofdm_lattice []
% ([dN dM] pins the OFDM pilot spacing, file suffix _lat<dN>x<dM>; default: the
% search of pilot_spacing_match, which gives (37, 2) at 120-360 km/h but (4, 19)
% at 60 km/h on the reduced EVA); tag ''; and the
% receiver knobs of run_mcrb_rx_eva (use_ici, Tcp_samples, E, mp_*, ce_mode,
% n_dd, pilot_boost_dB, pilot_frac, modtype, kappa_thresh, N, M, deltaf).
%
% Output: results/p2_<channel>_v<speed>_n<trials>_s<seed><tag>.mat with the
% run_mcrb_rx_eva schema (trial.tr, kappa, exitflag, ..., reduce_mcrb_rx_stats
% summaries) and the struct x:
%   x.eta_bar, x.eta0          nT x 4P x nS x nW   true / pseudo-true parameter
%   x.Gk, x.Gck                nT x 2P x 2P x nS x nW  noise term G and genie
%                              CRB at sigma2 = 1 on the [tau; nu] block (SI units)
%   x.err_energy, x.nsym_err   nT x nS x nW  sum |x - xhat|^2 and symbol errors (data cells)
%   x.ctl(c).name, .tr (nT x 6 x nS x nW: tau mcrb, bb, lb, nu mcrb, bb, lb),
%            .kappa, .exitflag, .pcvalid, .eta0, .err_energy
%   x.mmle.arms, .eta_hat (nT x 4P x nA x nS x nW, paths matched to the true
%            ones), .exitflag, .resnorm, .iters, .snr_db
%   x.sec (nS x nW wall seconds), ndone (nS x nW trials finished)
% A killed run restarts from CKPT_<same name> (checkpoints are written
% atomically after every chunk); a finished run is not repeated; an existing
% result is never written over.  Creating results/STOP_PHASE2 stops the run
% with an error after the current chunk.
    o = struct('tier', 500, 'v_kmh', 120, 'reduced', true, 'snr_db', 0:2:24, 'waveforms', {{'cp-ofdm', 'cp-otfs'}}, ...
        'seed', 9100, 'controls', {{'uniform', 'nearest'}}, 'mmle_snr', [], 'mmle_arms', {{'e2e', 'fresh', 'e2e_truth', 'fresh_truth'}}, 'csi', 'ce', ...
        'chunk', 40, 'resume', true, 'pass', [], 'frames', [], 'ofdm_lattice', [], 'tag', '', 'use_ici', true, 'Tcp_samples', 20, 'E', 10, 'mp_iter', 20, ...
        'mp_damp', 0.7, 'mp_mode', 'edge', 'ce_mode', 'dftdd', 'n_dd', 2, 'pilot_boost_dB', 0, 'pilot_frac', 0.20, ...
        'modtype', '4qam', 'kappa_thresh', Inf, 'N', 512, 'M', 128, 'deltaf', 15e3);
    o = parse_opts(o, varargin, 'run_phase2_eva');
    assert(any(strcmpi(o.csi, {'ce', 'genie'})), 'run_phase2_eva:csi', 'csi must be ce or genie');
    assert(all(ismember(o.controls, {'uniform', 'nearest'})), 'run_phase2_eva:controls', 'controls must be uniform and/or nearest');
    nT = tier_trials(o.tier);  o.tier = sprintf('n%d', nT);
    root = fileparts(fileparts(mfilename('fullpath')));
    N = o.N;  M = o.M;  Ts = 1/(N*o.deltaf);
    g = make_grid(N, M, o.deltaf, o.Tcp_samples*Ts, 'cp');
    ch = chan_eva(g, 'reduced', o.reduced, 'v_kmh', o.v_kmh);
    P = ch.P;  nS = numel(o.snr_db);  nW = numel(o.waveforms);  nC = numel(o.controls);
    assert(all(ismember(o.mmle_arms, {'e2e', 'fresh', 'e2e_truth', 'fresh_truth'})), 'run_phase2_eva:arms', 'unknown estimator arm');
    arms = {};  if ~isempty(o.mmle_snr), arms = o.mmle_arms; end
    nA = numel(arms);
    csi_sfx = '';  if strcmpi(o.csi, 'genie'), csi_sfx = '_geniecsi'; end
    if ~isempty(o.ofdm_lattice), csi_sfx = sprintf('%s_lat%dx%d', csi_sfx, o.ofdm_lattice(1), o.ofdm_lattice(2)); end
    stem = sprintf('p2_%s_v%g_%s_s%d%s%s', lower(strrep(ch.name, '-', '_')), o.v_kmh, o.tier, o.seed, csi_sfx, o.tag);
    f_out = fullfile(root, 'results', [stem '.mat']);
    f_ck  = fullfile(root, 'results', ['CKPT_' stem '.mat']);
    f_stop = fullfile(root, 'results', 'STOP_PHASE2');           % touch this file to stop cleanly after the current chunk
    sig = rmfield(o, {'chunk', 'resume', 'pass', 'frames'});
    if isempty(o.ofdm_lattice), sig = rmfield(sig, 'ofdm_lattice'); end     % checkpoints made before this option stay resumable
    nGoal = nT;  if ~isempty(o.frames), nGoal = min(nT, o.frames); end
    if isempty(o.pass), o.pass = o.snr_db; end
    [okp, pass_idx] = ismember(round(10*o.pass), round(10*o.snr_db));
    assert(all(okp), 'run_phase2_eva:pass', 'pass must be a subset of snr_db');
    if isfile(f_out)
        fprintf('run_phase2_eva: %s exists, nothing to do\n', f_out);
        if nargout, res = load(f_out); end
        return;
    end

    trn = {'tr_tau_mcrb', 'tr_tau_bb', 'tr_tau_lb', 'tr_tau_crb', 'tr_nu_mcrb', 'tr_nu_bb', 'tr_nu_lb', 'tr_nu_crb', ...
           'tr_tau_lb_bsc', 'tr_nu_lb_bsc', 'tr_tau_pil', 'tr_nu_pil'};
    nQ = numel(trn);
    if isfile(f_ck)
        assert(o.resume, 'run_phase2_eva:ckpt', 'Checkpoint %s exists and resume is false', f_ck);
        res = load(f_ck);
        assert(isequaln(res.sig, sig), 'run_phase2_eva:sig', 'Checkpoint %s was made with different options', f_ck);
        fprintf('run_phase2_eva: resuming %s (%d of %d trial-cells done)\n', stem, sum(res.x.ndone(:)), nT*nS*nW);
    else
        res = struct('schema_version', 'phase2_v1 (mcrb_rx_v2_arm_validity + x)');
        z3 = nan(nT, nS, nW);
        res.trial = struct('tr', nan(nT, nQ, nS, nW), 'kappa', z3, 'kappa_bsc', z3, 'exitflag', z3, 'exitflag_bsc', z3, ...
            'nerr', z3, 'nbits', z3, 'ser', z3, 'pcvalid', z3, 'pcvalid_bsc', z3, 'F_over_G', z3);
        lo = o;  lo.do_bsc = nC > 0;  lo.bsc_mode = '';  if nC > 0, lo.bsc_mode = o.controls{1}; end
        lo.statistic = 'median';  lo.gain_model = 'rayleigh';  lo.legacy = false;  lo.otfs_snr_cap = Inf;
        res.opts = lo;  res.sig = sig;  res.snr_db = o.snr_db;  res.waveforms = o.waveforms;  res.nT = nT;  res.g = g;  res.ch = ch;  res.trn = trn;
        x = struct('eta_bar', nan(nT, 4*P, nS, nW), 'eta0', nan(nT, 4*P, nS, nW), 'Gk', nan(nT, 2*P, 2*P, nS, nW), ...
            'Gck', nan(nT, 2*P, 2*P, nS, nW), 'err_energy', z3, 'nsym_err', z3, 'pt_iters', z3, ...
            'sec', zeros(nS, nW), 'ndone', zeros(nS, nW));
        x.kin_index = 'rows/columns of Gk, Gck: [tau_1..tau_P, nu_1..nu_P]';
        x.ctl = struct('name', {}, 'tr', {}, 'kappa', {}, 'exitflag', {}, 'pcvalid', {}, 'eta0', {}, 'err_energy', {});
        for c = 1:nC
            x.ctl(c) = struct('name', o.controls{c}, 'tr', nan(nT, 6, nS, nW), 'kappa', z3, 'exitflag', z3, 'pcvalid', z3, ...
                'eta0', nan(nT, 4*P, nS, nW), 'err_energy', z3);
        end
        x.ctl_trn = {'tr_tau_mcrb', 'tr_tau_bb', 'tr_tau_lb', 'tr_nu_mcrb', 'tr_nu_bb', 'tr_nu_lb'};
        x.mmle = struct('arms', {arms}, 'snr_db', o.mmle_snr, 'eta_hat', nan(nT, 4*P, max(nA, 1), nS, nW), ...
            'exitflag', nan(nT, max(nA, 1), nS, nW), 'resnorm', nan(nT, max(nA, 1), nS, nW), 'iters', nan(nT, max(nA, 1), nS, nW));
        res.x = x;
    end
    fprintf('Phase 2: %s at %g km/h, P = %d, %d trials, SNR %s, seed %d, csi %s, controls {%s}, MMLE at %s {%s}\n', ...
        ch.name, o.v_kmh, P, nT, mat2str(o.snr_db), o.seed, o.csi, strjoin(o.controls, ' '), mat2str(o.mmle_snr), strjoin(arms, ' '));

    tAll = tic;
    for si = pass_idx(:).'
        snr = o.snr_db(si);  N0 = 10^(-snr/10);  key = 1000 + round(10*snr);     % same key as run_mcrb_rx_eva / run_ber_eva
        do_mmle = nA > 0 && any(abs(o.mmle_snr - snr) < 1e-9);
        for iw = 1:nW
            wf = o.waveforms{iw};
            [rx, dom] = rx_for(wf, o);
            if res.x.ndone(si, iw) >= nGoal, continue; end
            while res.x.ndone(si, iw) < nGoal
                tt = res.x.ndone(si, iw) + (1:min(o.chunk, nGoal - res.x.ndone(si, iw)));
                out = cell(1, numel(tt));  tC = tic;
                parfor k = 1:numel(tt)
                    warning('off', 'MATLAB:nearlySingularMatrix');
                    out{k} = frame_eval(wf, dom, g, ch, N0, rx, o, key, tt(k), do_mmle, arms);
                end
                for k = 1:numel(tt)
                    t = tt(k);  q = out{k};
                    res.trial.tr(t, :, si, iw) = q.row;  res.trial.kappa(t, si, iw) = q.kappa;  res.trial.exitflag(t, si, iw) = q.exitflag;
                    res.trial.nerr(t, si, iw) = q.nerr;  res.trial.nbits(t, si, iw) = q.nbits;  res.trial.ser(t, si, iw) = q.ser;
                    res.trial.pcvalid(t, si, iw) = q.pcvalid;  res.trial.F_over_G(t, si, iw) = q.F_over_G;
                    res.x.eta_bar(t, :, si, iw) = q.eta_bar;  res.x.eta0(t, :, si, iw) = q.eta0;
                    res.x.Gk(t, :, :, si, iw) = q.Gk;  res.x.Gck(t, :, :, si, iw) = q.Gck;
                    res.x.err_energy(t, si, iw) = q.err_energy;  res.x.nsym_err(t, si, iw) = q.nsym_err;  res.x.pt_iters(t, si, iw) = q.pt_iters;
                    for c = 1:nC
                        res.x.ctl(c).tr(t, :, si, iw) = q.ctl(c).tr;  res.x.ctl(c).kappa(t, si, iw) = q.ctl(c).kappa;
                        res.x.ctl(c).exitflag(t, si, iw) = q.ctl(c).exitflag;  res.x.ctl(c).pcvalid(t, si, iw) = q.ctl(c).pcvalid;
                        res.x.ctl(c).eta0(t, :, si, iw) = q.ctl(c).eta0;  res.x.ctl(c).err_energy(t, si, iw) = q.ctl(c).err_energy;
                    end
                    if nC > 0
                        res.trial.kappa_bsc(t, si, iw) = q.ctl(1).kappa;  res.trial.exitflag_bsc(t, si, iw) = q.ctl(1).exitflag;
                        res.trial.pcvalid_bsc(t, si, iw) = q.ctl(1).pcvalid;
                    end
                    if do_mmle
                        res.x.mmle.eta_hat(t, :, :, si, iw) = q.mmle_eta;  res.x.mmle.exitflag(t, :, si, iw) = q.mmle_flag;
                        res.x.mmle.resnorm(t, :, si, iw) = q.mmle_res;  res.x.mmle.iters(t, :, si, iw) = q.mmle_it;
                    end
                end
                res.x.ndone(si, iw) = tt(end);  res.x.sec(si, iw) = res.x.sec(si, iw) + toc(tC);
                if res.x.ndone(si, iw) < nGoal, save_atomic(f_ck, res); end
                if isfile(f_stop)
                    save_atomic(f_ck, res);
                    error('run_phase2_eva:stop', 'Stop file %s found: checkpoint saved after %s %g dB trial %d; delete the file and call again to resume.', f_stop, wf, snr, tt(end));
                end
            end
            res = reduce_mcrb_rx_stats(res, res.opts);
            fprintf('  %-8s %2d dB %3d fr | BER %.2e SER %.2e | delay RMSE [ns] (median): LB %.3f  MCRB0 %.3f  CRB %.3f  ctl-LB %.3f  pilot %.3f | Doppler [Hz]: LB %.3f CRB %.3f | drop %.0f%% kappa<= %.1e%s (cell %.0f s, run %.0f s)\n', ...
                wf, snr, res.x.ndone(si, iw), res.ber(si, iw), res.ser(si, iw), 1e9*res.rmse_tau_lb(si, iw), 1e9*res.rmse_tau_mcrb(si, iw), 1e9*res.rmse_tau_crb(si, iw), ...
                1e9*res.rmse_tau_lb_bsc(si, iw), 1e9*res.rmse_tau_pil(si, iw), res.rmse_nu_lb(si, iw), res.rmse_nu_crb(si, iw), ...
                100*res.drop_frac(si, iw), res.kappa_max(si, iw), mmle_line(res, si, iw, do_mmle), res.x.sec(si, iw), toc(tAll));
            save_atomic(f_ck, res);
        end
    end
    if any(res.x.ndone(:) < nT)
        fprintf('PHASE2_PASS_COMPLETED %s (SNR %s, %d frames); cells: %d of %d complete, %d partly done, checkpoint kept\n', ...
            stem, mat2str(o.pass), nGoal, nnz(res.x.ndone == nT), nS*nW, nnz(res.x.ndone > 0 & res.x.ndone < nT));
        return;
    end
    res.seconds = sum(res.x.sec(:));
    assert(~isfile(f_out), 'run_phase2_eva:exists', 'Refusing to overwrite %s', f_out);
    save_atomic(f_out, res);
    delete(f_ck);
    fprintf('PHASE2_RUN_COMPLETED %s seconds=%.0f\n', stem, res.seconds);
end

% ---- one frame: every arm of one (waveform, SNR, trial) ------------------------
function q = frame_eval(wf, dom, g, ch, N0, rx, o, key, t, do_mmle, arms)
    mdm = modem(o.modtype);
    rs = crn_stream(o.seed, key, t);
    f = ber_frame(wf, g, ch, N0, rx, rs);
    cr = f.cr;  X = f.X;  Xh = f.Xhat;  P = cr.P;
    idx = eta_index(P);  kin = [idx.tau(:); idx.nu(:)];
    Jg = fim(X, wf, g, cr, o.use_ici, 1);                      % genie FIM on the transmitted grid
    [pc, info, eta0] = mcrb_from_grids(wf, X, Xh, g, cr, o.use_ici, dom, Jg);
    lb = mcrb_eval(pc, N0);
    Sp = crb_pilot_only(g, cr, o.use_ici, o.pilot_frac, N0, 'cp-ofdm');
    q.row = [lb.tr_tau_mcrb, lb.tr_tau_bb, lb.tr_tau_lb, lb.tr_tau_crb, ...
             lb.tr_nu_mcrb, lb.tr_nu_bb, lb.tr_nu_lb, lb.tr_nu_crb, NaN, NaN, Sp.tr_tau, Sp.tr_nu];
    q.kappa = pc.kappa;  q.exitflag = info.exitflag;  q.pcvalid = double(pc.valid);  q.pt_iters = info.iterations;
    q.nerr = f.nerr;  q.nbits = f.nbits;  q.ser = f.ser;
    q.F_over_G = sum(real(diag(pc.F(idx.tau, idx.tau)))) / max(N0*sum(real(diag(pc.G(idx.tau, idx.tau)))), eps);
    q.eta_bar = cr.eta_bar(:).';  q.eta0 = eta0(:).';
    q.Gk = reshape(pc.G(kin, kin), 1, 2*P, 2*P);  q.Gck = reshape(pc.Gc(kin, kin), 1, 2*P, 2*P);
    dm = f.masks.data;
    q.err_energy = sum(abs(X(dm) - Xh(dm)).^2);  q.nsym_err = nnz(X(dm) ~= Xh(dm));
    % Control arms: the first one consumes the stream exactly as run_mcrb_rx_eva
    % does; the stream state is rewound before each further arm, so all arms move
    % the errors to the SAME cells and differ only in the wrong symbol's value.
    st = rs.State;
    q.ctl = struct('tr', {}, 'kappa', {}, 'exitflag', {}, 'pcvalid', {}, 'eta0', {}, 'err_energy', {});
    for c = 1:numel(o.controls)
        rs.State = st;
        Xb = bsc_match(X, Xh, dm, mdm, rs, o.controls{c});
        [pcb, infob, eta0b] = mcrb_from_grids(wf, X, Xb, g, cr, o.use_ici, dom, Jg);
        lbb = mcrb_eval(pcb, N0);
        q.ctl(c) = struct('tr', [lbb.tr_tau_mcrb, lbb.tr_tau_bb, lbb.tr_tau_lb, lbb.tr_nu_mcrb, lbb.tr_nu_bb, lbb.tr_nu_lb], ...
            'kappa', pcb.kappa, 'exitflag', infob.exitflag, 'pcvalid', double(pcb.valid), 'eta0', eta0b(:).', ...
            'err_energy', sum(abs(X(dm) - Xb(dm)).^2));
        if c == 1, q.row(9) = lbb.tr_tau_lb;  q.row(10) = lbb.tr_nu_lb; end
    end
    % Mismatched ML estimator on the detected reference.
    nA = numel(arms);
    q.mmle_eta = nan(1, 4*P, max(nA, 1));  q.mmle_flag = nan(1, max(nA, 1));  q.mmle_res = q.mmle_flag;  q.mmle_it = q.mmle_flag;
    if do_mmle
        rs2 = crn_stream(o.seed + 104729, key, t);              % fresh noise, shared by the waveforms
        Z2 = sqrt(N0/2) * (randn(rs2, g.N, g.M) + 1j*randn(rs2, g.N, g.M));
        if strcmp(wf, 'cp-ofdm'), Yf = io_cp_ofdm(X, g, cr, o.use_ici) + Z2; else, Yf = io_cp_otfs(X, g, cr, o.use_ici) + Z2; end
        for a = 1:nA
            mo = struct('domain', dom, 'method', 'lsq');
            if endsWith(arms{a}, '_truth'), mo.eta_init = cr.eta_bar; end
            if startsWith(arms{a}, 'e2e'), Yobs = f.Y; else, Yobs = Yf; end
            [eh, mi] = mmle_cp(Yobs, Xh, g, P, o.use_ici, mo);
            q.mmle_eta(1, :, a) = match_paths(eh, cr, g);
            q.mmle_flag(a) = mi.exitflag;  q.mmle_res(a) = mi.resnorm;  q.mmle_it(a) = mi.iterations;
        end
    end
end

function [pc, info, eta0] = mcrb_from_grids(wf, X, Xh, g, cr, use_ici, dom, Jg)
    [eta0, info] = pseudo_true(X, Xh, g, cr, use_ici, struct('domain', dom));
    [a0, t0, n0] = eta_unpack(eta0);
    ce = cr;  ce.alpha = a0;  ce.tau = t0;  ce.nu = n0;
    if strcmp(wf, 'cp-ofdm')
        [A, ~, Jt, b] = mcrb_AB_cp_ofdm(X, Xh, g, ce, use_ici, 1, cr);
    else
        [A, ~, Jt, b] = mcrb_AB_cp_otfs(X, Xh, g, ce, use_ici, 1, cr);
    end
    pc = mcrb_assemble(A, Jt, b, Jg, cr.eta_bar, eta0, g);
end

function eh = match_paths(eh, cr, g)
% Reorder the estimated paths to the true ones (minimum total squared distance
% in resolution units), so that eta_hat - eta_bar is an error vector.
    [a, t, n] = eta_unpack(eh);
    C = ((t(:) - cr.tau(:).')/g.bin_tau).^2 + ((n(:) - cr.nu(:).')/g.bin_nu).^2;   % rows: estimates, columns: true
    C(~isfinite(C)) = 1e11;
    mp = matchpairs(C, 1e12);
    perm = zeros(numel(t), 1);  perm(mp(:, 2)) = mp(:, 1);
    eh = eta_pack(a(perm), t(perm), n(perm)).';
end

function s = mmle_line(res, si, iw, do_mmle)
% Normalized squared error of the estimator arms: median over valid trials of
% |theta_hat - theta_bar|^2 / tr(LB) per block (1 = the bound predicts the error).
    s = '';
    if ~do_mmle, return; end
    P = res.ch.P;  it = 3:4:4*P;  in = 4:4:4*P;
    v = res.trial.valid(:, si, iw);
    lbt = res.trial.tr(v, 3, si, iw);  lbn = res.trial.tr(v, 7, si, iw);
    eb = res.x.eta_bar(v, :, si, iw);
    s = ' | MMLE err^2/LB median (tau, nu):';
    for a = 1:numel(res.x.mmle.arms)
        d = res.x.mmle.eta_hat(v, :, a, si, iw) - eb;
        s = [s sprintf(' %s %.2f %.2f', res.x.mmle.arms{a}, median(sum(d(:, it).^2, 2)./lbt, 'omitnan'), median(sum(d(:, in).^2, 2)./lbn, 'omitnan'))]; %#ok<AGROW>
    end
end

function [rx, dom] = rx_for(wf, o)
    rx = struct('E', o.E, 'mp_iter', o.mp_iter, 'mp_damp', o.mp_damp, 'mp_mode', o.mp_mode, 'n_dd', o.n_dd, ...
        'pilot_boost_dB', o.pilot_boost_dB, 'pilot', 'matched', 'use_ici', o.use_ici, 'ce_mode', o.ce_mode, ...
        'modtype', o.modtype, 'ofdm_lattice', o.ofdm_lattice);
    genie = strcmpi(o.csi, 'genie');
    if strcmp(wf, 'cp-ofdm')
        dom = 'tf';  if genie, rx.ofdm_csi = 'idealbanded'; else, rx.ofdm_csi = 'ce'; end
    else
        dom = 'dd';  if genie, rx.otfs_csi = 'ideal'; else, rx.otfs_csi = 'ce'; end
    end
end

function save_atomic(fn, res)
    tmp = [fn '.tmp.mat'];
    save(tmp, '-struct', 'res');                    % v7: readable by scipy.io.loadmat (the exporter)
    movefile(tmp, fn, 'f');
end

function nT = tier_trials(tier)
    if isnumeric(tier), nT = tier; return; end
    switch lower(tier)
        case 'smoke', nT = 20;
        case 'dev',   nT = 200;
        case 'paper', nT = 2000;
        otherwise, error('tier must be smoke|dev|paper or a number');
    end
end
