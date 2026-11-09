function res = run_ber_fig9(varargin)
% RUN_BER_FIG9  Raviteja et al. TVT'19 Fig. 9 on the EVA channel (plan F2.1):
% BER vs SNR_d for OTFS with embedded-pilot threshold CE + MP at
% SNR_p = {40, 45, 50} dB, OTFS with ideal CSI + MP, and OFDM with ideal CSI
% (banded ICI MMSE).  Raviteja's noise convention: sigma^2 = 1, data symbols
% scaled by sqrt(SNR_d), pilot |x_p|^2 = SNR_p, threshold 3 sigma on |y|.
% This is the acceptance test of the MP detector / damping.
%
%   res = run_ber_fig9('tier', 'smoke', ...)
% Options: tier 'smoke'(10)|'dev'(50)|'paper'(200) trials, snr_d_db 10:2:18,
% snr_p_db [40 45 50], E 10, mp_iter 20, mp_damp 0.7, seed 2019, Tcp_samples 20.
    o = struct('tier', 'smoke', 'snr_d_db', 10:2:18, 'snr_p_db', [40 45 50], 'E', 10, ...
        'mp_iter', 20, 'mp_damp', 0.7, 'seed', 2019, 'Tcp_samples', 20, 'mp_mode', 'edge');
    for i = 1:2:numel(varargin), o.(varargin{i}) = varargin{i+1}; end
    root = fileparts(fileparts(mfilename('fullpath')));
    nT = tier_trials(o.tier, [10 50 200]);
    if isnumeric(o.tier), o.tier = sprintf('n%d', nT); end
    N = 512; M = 128; df = 15e3;  Ts = 1/(N*df);
    g  = make_grid(N, M, df, o.Tcp_samples*Ts, 'cp');
    ch = chan_eva(g);                                   % full EVA, Rayleigh, snapped
    nS = numel(o.snr_d_db);  nP = numel(o.snr_p_db);
    curves = [{'OTFS ideal'}, arrayfun(@(p) sprintf('OTFS CE SNR_p=%d', p), o.snr_p_db, 'UniformOutput', false), {'OFDM ideal'}];
    nC = numel(curves);
    nerr = zeros(nT, nS, nC);  nbits = zeros(nT, nS, nC);
    fprintf('Fig. 9: %s, %d trials, SNR_d %s dB, SNR_p %s dB, E=%d, MP %d iters damp %.1f (%s)\n', ...
        ch.name, nT, mat2str(o.snr_d_db), mat2str(o.snr_p_db), o.E, o.mp_iter, o.mp_damp, o.mp_mode);
    t0 = tic;
    parfor t = 1:nT
        le = zeros(nS, nC);  lb = zeros(nS, nC);
        for si = 1:nS
            amp = sqrt(10^(o.snr_d_db(si)/10));
            base = struct('data_amp', amp, 'E', o.E, 'mp_iter', o.mp_iter, 'mp_damp', o.mp_damp, 'mp_mode', o.mp_mode);
            % OTFS ideal CSI, full data frame
            rx = base;  rx.pilot = 'none';  rx.otfs_csi = 'ideal';
            f = ber_frame('cp-otfs', g, ch, 1, rx, crn_stream(o.seed, si, t));
            le(si, 1) = f.nerr;  lb(si, 1) = f.nbits;
            % OTFS CE at each pilot SNR
            for ip = 1:nP
                rx = base;  rx.pilot = 'explicit';  rx.x_p = sqrt(10^(o.snr_p_db(ip)/10));  rx.otfs_csi = 'ce';  rx.thr = 3;
                f = ber_frame('cp-otfs', g, ch, 1, rx, crn_stream(o.seed, si, t));
                le(si, 1+ip) = f.nerr;  lb(si, 1+ip) = f.nbits;
            end
            % OFDM ideal CSI (known ICI, banded MMSE), no pilots
            rx = base;  rx.pilot = 'none';  rx.ofdm_csi = 'idealbanded';
            f = ber_frame('cp-ofdm', g, ch, 1, rx, crn_stream(o.seed, si, t));
            le(si, nC) = f.nerr;  lb(si, nC) = f.nbits;
        end
        nerr(t, :, :) = le;  nbits(t, :, :) = lb;
        fprintf('  trial %d/%d (%.0f s)\n', t, nT, toc(t0));
    end
    ber = reshape(sum(nerr, 1) ./ max(sum(nbits, 1), 1), nS, nC);   % nS x nC (reshape, not squeeze: nS may be 1)
    res = struct('ber', ber, 'curves', {curves}, 'snr_d_db', o.snr_d_db, 'snr_p_db', o.snr_p_db, ...
        'nT', nT, 'opts', o, 'g', g, 'ch', ch, 'nerr', nerr, 'nbits', nbits, 'seconds', toc(t0));
    % reference (impl/otfs_ber/README.md, 48 trials) for the shared SNR points
    ref.snr = 10:2:18;
    ref.otfs_ideal = [1.8e-2 3.2e-3 9.7e-4 2.9e-4 1.2e-4];
    ref.ce50 = [1.9e-2 3.3e-3 1.9e-3 4.4e-4 9.0e-5];
    ref.ce45 = [2.5e-2 1.1e-2 2.0e-3 9.5e-4 3.4e-4];
    ref.ce40 = [3.2e-2 1.6e-2 1.1e-2 4.8e-3 3.4e-3];
    ref.ofdm  = [4.5e-2 3.0e-2 2.0e-2 1.3e-2 8.1e-3];
    res.ref = ref;
    fprintf('\nBER table (rows: SNR_d dB):\n  %-8s', 'SNR_d');  fprintf('%-20s', curves{:});  fprintf('\n');
    for si = 1:nS
        fprintf('  %-8d', o.snr_d_db(si));  fprintf('%-20.2e', ber(si, :));  fprintf('\n');
    end
    fprintf('Acceptance (F2.1): ordering CE40 > CE45 > CE50 >= ideal at 18 dB: %s; OFDM top: %s\n', ...
        pw(all(diff(ber(end, 2:1+nP)) < 0) && ber(end, 1+nP) >= 0.8*ber(end, 1)), pw(all(ber(:, nC) >= ber(:, 1))));
    % Compare with the old reproduction only where BOTH runs are well measured
    % (BER > 20x the measurement floor 1/bits): the OTFS tail is dominated by
    % rare deep-fade channels and needs >> 50 trials to resolve.  On identical
    % channel draws the two detectors are bit-identical and the two channel
    % models agree within ~2x (scratch check, 2026-09-07).
    [~, ia, ib] = intersect(o.snr_d_db, ref.snr);
    if ~isempty(ia)
        fl = 20 ./ reshape(sum(nbits(:, ia, [1 nC]), 1), numel(ia), 2);   % 20 x floor
        refv = [ref.otfs_ideal(ib).' ref.ofdm(ib).'];
        rat = ber(ia, [1 nC]) ./ refv;
        ok = ber(ia, [1 nC]) > fl & refv > 20/(48*2*N*M);
        fprintf('  vs impl/otfs_ber reference (well-measured points only): OTFS-ideal ratio %s, OFDM ratio %s (expect within x2)  %s\n', ...
            mat2str(rat(:,1), 2), mat2str(rat(:,2), 2), pw(all(rat(ok) > 0.5 & rat(ok) < 2)));
    end
    if ~exist(fullfile(root, 'results'), 'dir'), mkdir(fullfile(root, 'results')); end
    save(fullfile(root, 'results', sprintf('ber_fig9_%s.mat', o.tier)), '-struct', 'res');
    plot_fig9(res);
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
