function res = run_phase2_mmle(varargin)
% RUN_PHASE2_MMLE  Estimator-only supplement to a run_phase2_eva result: the
% mismatched ML estimator (mmle_cp 'lsq') at SNR points of that run which have
% bounds but no estimator arms.  The frames are regenerated from the run's own
% options (its stored signature) and the same common-random-number streams, and
% every regenerated frame is checked against the stored true parameter and bit
% error count, so each estimate belongs to the frame whose bounds the run holds.
% The source file is only read; the estimates go to a file of their own.
%
%   run_phase2_mmle('source', 'CKPT_p2_eva_reduced_v120_n500_s9100', 'snr_db', [2 6 10 24])
%
% Options (defaults): source (file in results/, with or without .mat); snr_db
% (SNR points of the source to estimate at); frames 250 (frames 1..frames of
% every point); arms {'e2e','fresh'} (as in run_phase2_eva: 'e2e' = the
% observation that produced the decisions, 'fresh' = the same decisions with an
% independent noise draw, suffix '_truth' = fit started at the true parameter;
% the default arms start from the periodogram with successive cancellation);
% check_frames 8 (frames per waveform recomputed at check_snr and compared with
% the estimates stored in the source); check_snr [] (default: the first SNR of
% the source that has estimator arms; the check is skipped when there is none);
% chunk 40; resume true; tag ''.
%
% Output: results/p2mmle_<source stem><tag>.mat (v7, readable by scipy) with
%   snr_db, waveforms, arms, frames
%   eta_hat   frames x 4P x nA x nS x nW   estimates, paths matched to the true ones
%   exitflag, resnorm, iters   frames x nA x nS x nW
%   eta_bar   frames x 4P x nS x nW   true parameter of the regenerated frame
%   nerr      frames x nS x nW        bit errors of the regenerated frame
%   ndone, sec (nS x nW), check (the comparison with the source), sig, source
% A killed run restarts from CKPT_<same name>; a finished run is not repeated;
% an existing result is never written over.  Creating results/STOP_PHASE2 stops
% the run with an error after the current chunk.
%
% The frame and estimator code below repeats frame_eval, rx_for and match_paths
% of run_phase2_eva (local functions there); the check against the stored
% estimates guards the copy.
    o = struct('source', '', 'snr_db', [], 'frames', 250, 'arms', {{'e2e', 'fresh'}}, 'check_frames', 8, ...
        'check_snr', [], 'chunk', 40, 'resume', true, 'tag', '');
    o = parse_opts(o, varargin, 'run_phase2_mmle');
    assert(~isempty(o.source) && ~isempty(o.snr_db), 'run_phase2_mmle:args', 'source and snr_db are required');
    assert(all(ismember(o.arms, {'e2e', 'fresh', 'e2e_truth', 'fresh_truth'})), 'run_phase2_mmle:arms', 'unknown estimator arm');
    root = fileparts(fileparts(mfilename('fullpath')));
    [~, sname] = fileparts(char(o.source));
    f_src = fullfile(root, 'results', [sname '.mat']);
    assert(isfile(f_src), 'run_phase2_mmle:source', 'Source %s not found', f_src);
    stem = ['p2mmle_' regexprep(sname, '^(CKPT_)?p2_', '') o.tag];
    f_out = fullfile(root, 'results', [stem '.mat']);
    f_ck  = fullfile(root, 'results', ['CKPT_' stem '.mat']);
    f_stop = fullfile(root, 'results', 'STOP_PHASE2');
    if isfile(f_out)
        fprintf('run_phase2_mmle: %s exists, nothing to do\n', f_out);
        if nargout, res = load(f_out); end
        return;
    end

    src = load(f_src, 'sig', 'snr_db', 'waveforms', 'nT', 'g', 'ch', 'x', 'trial');
    so = src.sig;  if ~isfield(so, 'ofdm_lattice'), so.ofdm_lattice = []; end
    Ts = 1/(so.N*so.deltaf);
    g = make_grid(so.N, so.M, so.deltaf, so.Tcp_samples*Ts, 'cp');
    ch = chan_eva(g, 'reduced', so.reduced, 'v_kmh', so.v_kmh);
    assert(isequaln(g, src.g) && isequaln(ch, src.ch), 'run_phase2_mmle:grid', 'The grid or channel rebuilt from the source options differs from the stored one');
    P = ch.P;  wfs = src.waveforms;  nW = numel(wfs);  nS = numel(o.snr_db);  nA = numel(o.arms);  nF = o.frames;
    [ok, s_src] = ismember(round(10*o.snr_db), round(10*src.snr_db));
    assert(all(ok), 'run_phase2_mmle:snr', 'snr_db must be SNR points of the source');
    assert(all(all(src.x.ndone(s_src, :) >= nF)), 'run_phase2_mmle:frames', 'The source has fewer than %d frames at one of these SNR points', nF);
    sig = struct('source', sname, 'snr_db', o.snr_db, 'frames', nF, 'arms', {o.arms}, 'tag', o.tag, 'src_sig', src.sig);

    if isfile(f_ck)
        assert(o.resume, 'run_phase2_mmle:ckpt', 'Checkpoint %s exists and resume is false', f_ck);
        res = load(f_ck);
        assert(isequaln(res.sig, sig), 'run_phase2_mmle:sig', 'Checkpoint %s was made with different options', f_ck);
        fprintf('run_phase2_mmle: resuming %s (%d of %d frame-cells done)\n', stem, sum(res.ndone(:)), nF*nS*nW);
    else
        d = dir(f_src);
        res = struct('schema_version', 'phase2_mmle_v1', 'sig', sig, 'source', struct('file', [sname '.mat'], 'bytes', d.bytes, 'date', d.date), ...
            'snr_db', o.snr_db, 'waveforms', {wfs}, 'arms', {o.arms}, 'frames', nF, ...
            'eta_hat', nan(nF, 4*P, nA, nS, nW), 'exitflag', nan(nF, nA, nS, nW), 'resnorm', nan(nF, nA, nS, nW), 'iters', nan(nF, nA, nS, nW), ...
            'eta_bar', nan(nF, 4*P, nS, nW), 'nerr', nan(nF, nS, nW), 'ndone', zeros(nS, nW), 'sec', zeros(nS, nW), 'check', struct([]));
    end
    fprintf('Phase 2 estimator supplement: %s, %s at %g km/h, P = %d, %d frames, SNR %s, arms {%s}\n', ...
        sname, ch.name, so.v_kmh, P, nF, mat2str(o.snr_db), strjoin(o.arms, ' '));

    % ---- the copy of the frame code must reproduce the estimates the source holds ----
    if isempty(res.check)
        res.check = check_against_source(src, so, g, ch, wfs, o);
        save_atomic(f_ck, res);
    end

    tAll = tic;
    for si = 1:nS
        snr = o.snr_db(si);  N0 = 10^(-snr/10);  key = 1000 + round(10*snr);
        for iw = 1:nW
            wf = wfs{iw};
            while res.ndone(si, iw) < nF
                tt = res.ndone(si, iw) + (1:min(o.chunk, nF - res.ndone(si, iw)));
                out = cell(1, numel(tt));  tC = tic;
                parfor k = 1:numel(tt)
                    warning('off', 'MATLAB:nearlySingularMatrix');
                    out{k} = frame_mmle(wf, g, ch, N0, so, key, tt(k), o.arms);
                end
                for k = 1:numel(tt)
                    t = tt(k);  q = out{k};
                    same = isequal(q.eta_bar, src.x.eta_bar(t, :, s_src(si), iw)) && q.nerr == src.trial.nerr(t, s_src(si), iw);
                    assert(same, 'run_phase2_mmle:frame', '%s %g dB frame %d is not the frame of the source (true parameter or bit errors differ)', wf, snr, t);
                    res.eta_hat(t, :, :, si, iw) = q.eta;  res.exitflag(t, :, si, iw) = q.flag;
                    res.resnorm(t, :, si, iw) = q.res;  res.iters(t, :, si, iw) = q.it;
                    res.eta_bar(t, :, si, iw) = q.eta_bar;  res.nerr(t, si, iw) = q.nerr;
                end
                res.ndone(si, iw) = tt(end);  res.sec(si, iw) = res.sec(si, iw) + toc(tC);
                save_atomic(f_ck, res);
                if isfile(f_stop)
                    error('run_phase2_mmle:stop', 'Stop file %s found: checkpoint saved after %s %g dB frame %d; delete the file and call again to resume.', f_stop, wf, snr, tt(end));
                end
            end
            fprintf('  %-8s %2d dB %3d fr | estimator err^2/LB median (tau, nu):%s (cell %.0f s, run %.0f s)\n', ...
                wf, snr, res.ndone(si, iw), ratio_line(res, src, si, s_src(si), iw, P), res.sec(si, iw), toc(tAll));
        end
    end
    res.seconds = sum(res.sec(:));
    assert(~isfile(f_out), 'run_phase2_mmle:exists', 'Refusing to overwrite %s', f_out);
    save_atomic(f_out, res);
    delete(f_ck);
    fprintf('PHASE2_MMLE_COMPLETED %s seconds=%.0f\n', stem, res.seconds);
end

% ---- one frame: detection, then the estimator arms (copy of frame_eval's estimator part) ----
function q = frame_mmle(wf, g, ch, N0, so, key, t, arms)
    [rx, dom] = rx_for(wf, so);
    rs = crn_stream(so.seed, key, t);
    f = ber_frame(wf, g, ch, N0, rx, rs);
    cr = f.cr;  X = f.X;  Xh = f.Xhat;  P = cr.P;  nA = numel(arms);
    q.eta_bar = cr.eta_bar(:).';  q.nerr = f.nerr;
    q.eta = nan(1, 4*P, nA);  q.flag = nan(1, nA);  q.res = q.flag;  q.it = q.flag;
    rs2 = crn_stream(so.seed + 104729, key, t);                 % fresh noise, shared by the waveforms
    Z2 = sqrt(N0/2) * (randn(rs2, g.N, g.M) + 1j*randn(rs2, g.N, g.M));
    if strcmp(wf, 'cp-ofdm'), Yf = io_cp_ofdm(X, g, cr, so.use_ici) + Z2; else, Yf = io_cp_otfs(X, g, cr, so.use_ici) + Z2; end
    for a = 1:nA
        mo = struct('domain', dom, 'method', 'lsq');
        if endsWith(arms{a}, '_truth'), mo.eta_init = cr.eta_bar; end
        if startsWith(arms{a}, 'e2e'), Yobs = f.Y; else, Yobs = Yf; end
        [eh, mi] = mmle_cp(Yobs, Xh, g, P, so.use_ici, mo);
        q.eta(1, :, a) = match_paths(eh, cr, g);
        q.flag(a) = mi.exitflag;  q.res(a) = mi.resnorm;  q.it(a) = mi.iterations;
    end
end

function c = check_against_source(src, so, g, ch, wfs, o)
% Recompute the estimator arms that the source also holds, on its first frames
% at one SNR point, and compare (relative difference per frame and arm).
    c = struct('snr_db', NaN, 'frames', 0, 'arms', {{}}, 'max_rel_diff', NaN, 'skipped', true);
    if o.check_frames < 1 || ~isfield(src.x, 'mmle') || isempty(src.x.mmle.arms), return; end
    arms = intersect(o.arms, src.x.mmle.arms, 'stable');
    snr = o.check_snr;  if isempty(snr), snr = src.x.mmle.snr_db(1); end
    si = find(abs(src.snr_db - snr) < 1e-9, 1);
    if isempty(arms) || isempty(si) || ~any(abs(src.x.mmle.snr_db - snr) < 1e-9), return; end
    [~, ia] = ismember(arms, src.x.mmle.arms);
    N0 = 10^(-snr/10);  key = 1000 + round(10*snr);  nC = o.check_frames;  nW = numel(wfs);
    dmax = zeros(nC, nW);
    for iw = 1:nW
        wf = wfs{iw};  out = cell(1, nC);
        parfor t = 1:nC
            warning('off', 'MATLAB:nearlySingularMatrix');
            out{t} = frame_mmle(wf, g, ch, N0, so, key, t, arms);
        end
        for t = 1:nC
            ref = squeeze(src.x.mmle.eta_hat(t, :, ia, si, iw));  got = squeeze(out{t}.eta);
            assert(isequal(out{t}.eta_bar, src.x.eta_bar(t, :, si, iw)), 'run_phase2_mmle:check', 'check frame %d (%s) is not the frame of the source', t, wf);
            dmax(t, iw) = max(abs(got(:) - ref(:)) ./ max(abs(ref(:)), realmin));
        end
    end
    c = struct('snr_db', snr, 'frames', nC, 'arms', {arms}, 'max_rel_diff', max(dmax(:)), 'skipped', false);
    fprintf('  check at %g dB: %d frames per waveform, arms {%s}, largest relative difference to the stored estimates %.2e\n', ...
        snr, nC, strjoin(arms, ' '), c.max_rel_diff);
    assert(c.max_rel_diff <= 1e-6, 'run_phase2_mmle:check', 'The regenerated estimates differ from those stored in the source (%.2e)', c.max_rel_diff);
end

function s = ratio_line(res, src, si, ss, iw, P)
% Median over the source's valid frames of |theta_hat - theta_bar|^2 / tr(LB) per block.
    it = 3:4:4*P;  in = 4:4:4*P;  nF = res.frames;
    v = src.trial.pcvalid(1:nF, ss, iw) > 0 & src.trial.exitflag(1:nF, ss, iw) > 0;
    lbt = src.trial.tr(1:nF, 3, ss, iw);  lbn = src.trial.tr(1:nF, 7, ss, iw);
    s = '';
    for a = 1:numel(res.arms)
        d = res.eta_hat(:, :, a, si, iw) - res.eta_bar(:, :, si, iw);
        s = [s sprintf(' %s %.2f %.2f', res.arms{a}, median(sum(d(v, it).^2, 2)./lbt(v), 'omitnan'), median(sum(d(v, in).^2, 2)./lbn(v), 'omitnan'))]; %#ok<AGROW>
    end
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
