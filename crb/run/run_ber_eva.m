function res = run_ber_eva(varargin)
% RUN_BER_EVA  Main BER vs SNR on the EVA channel (plan F2.2) with matched,
% data-referenced pilot budgets: OFDM {dftdd, dft, interp} CE + one-tap MMSE,
% OFDM ideal-CSI one-tap and ideal banded-ICI MMSE, OTFS embedded-pilot CE +
% MP and OTFS ideal-taps + MP.  N0 = 1/SNR, unit-power data.
%
%   res = run_ber_eva('tier', 'smoke', 'channel', 'eva', ...)
% Options (defaults): tier 'smoke'(20)|'dev'(200)|'paper'(2000) trials;
% snr_db 0:1:25; channel 'eva' | 'eva-reduced' | 'awgn' | 'flat_rayleigh';
% gain_model 'rayleigh'; doppler_scale 1; use_ici true; Tcp_samples 20;
% configs (see below); seed 9100; pilot_boost_dB 0; E 10; mp_iter 20;
% mp_damp 0.7; mp_mode 'edge'; n_dd 2; otfs_snr_cap Inf (skip OTFS above it).
% Saves results/ber_<channel>_<tier>.mat and plots with plot_ber_vs_snr.
    o = struct('tier', 'smoke', 'snr_db', 0:1:25, 'channel', 'eva', 'gain_model', 'rayleigh', ...
        'doppler_scale', 1, 'use_ici', true, 'Tcp_samples', 20, ...
        'configs', {{'ofdm-dftdd', 'ofdm-dft', 'ofdm-interp', 'ofdm-ideal1tap', 'ofdm-idealbanded', 'otfs-ce', 'otfs-ideal'}}, ...
        'seed', 9100, 'pilot_boost_dB', 0, 'E', 10, 'mp_iter', 20, 'mp_damp', 0.7, 'mp_mode', 'edge', 'n_dd', 2, ...
        'otfs_snr_cap', Inf, 'tag', '');
    o = parse_opts(o, varargin, 'run_ber_eva');
    root = fileparts(fileparts(mfilename('fullpath')));
    nT = tier_trials(o.tier, [20 200 2000]);
    if isnumeric(o.tier), o.tier = sprintf('n%d', nT); end
    N = 512; M = 128; df = 15e3;  Ts = 1/(N*df);
    g = make_grid(N, M, df, o.Tcp_samples*Ts, 'cp');
    switch lower(o.channel)
        case 'eva',          ch = chan_eva(g, 'gain_model', o.gain_model, 'doppler_scale', o.doppler_scale);
        case 'eva-reduced',  ch = chan_eva(g, 'reduced', true, 'gain_model', o.gain_model, 'doppler_scale', o.doppler_scale);
        case 'eva-reduced-fast', ch = chan_eva(g, 'variant', 'reduced-fast', 'gain_model', o.gain_model, 'doppler_scale', o.doppler_scale);
        case 'awgn',         ch = chan_controls(g, 'awgn');
        case 'flat_rayleigh', ch = chan_controls(g, 'flat_rayleigh');
        otherwise, error('run_ber_eva:channel', 'unknown channel %s', o.channel);
    end
    if strcmp(ch.family, 'control')
        % a flat channel has l_tau = 0: give the pilot frames a minimal guard so
        % the OTFS pilot band and the OFDM lattice still exist
        ch.l_tau = 1;  ch.k_nu = max(ch.k_nu, 1);
    end
    nS = numel(o.snr_db);  nC = numel(o.configs);
    nerr = zeros(nT, nS, nC);  nbits = zeros(nT, nS, nC);  nserr = zeros(nT, nS, nC);  nsym = zeros(nT, nS, nC);
    fprintf('BER run: %s (%s gains, doppler x%g, ici=%d), %d trials, %d SNR points, configs: %s\n', ...
        ch.name, ch.gain_model, ch.doppler_scale, o.use_ici, nT, nS, strjoin(o.configs, ' '));
    cfgs = o.configs;  snr_db = o.snr_db;
    t0 = tic;
    parfor t = 1:nT
        le = zeros(nS, nC);  lb = zeros(nS, nC);  lse = zeros(nS, nC);  lsy = zeros(nS, nC);
        for si = 1:nS
            N0 = 10^(-snr_db(si)/10);
            for ci = 1:nC
                [wf, rx] = config_to_rx(cfgs{ci}, o);
                if strcmp(wf, 'cp-otfs') && snr_db(si) > o.otfs_snr_cap, continue; end
                % substream keyed on the SNR VALUE (not its grid index) so that
                % run_mcrb_rx_eva with the same seed sees the same draws at equal SNR
                f = ber_frame(wf, g, ch, N0, rx, crn_stream(o.seed, 1000 + round(10*snr_db(si)), t));
                le(si, ci) = f.nerr;  lb(si, ci) = f.nbits;
                lse(si, ci) = f.ser * nnz(f.masks.data);  lsy(si, ci) = nnz(f.masks.data);
            end
        end
        nerr(t, :, :) = le;  nbits(t, :, :) = lb;  nserr(t, :, :) = lse;  nsym(t, :, :) = lsy;
        fprintf('  trial %d/%d (%.0f s)\n', t, nT, toc(t0));
    end
    ber = squeeze(sum(nerr, 1) ./ max(sum(nbits, 1), 1));  ber(squeeze(sum(nbits, 1)) == 0) = NaN;
    ser = squeeze(sum(nserr, 1) ./ max(sum(nsym, 1), 1));  ser(squeeze(sum(nsym, 1)) == 0) = NaN;
    res = struct('ber', ber, 'ser', ser, 'configs', {o.configs}, 'snr_db', o.snr_db, 'nT', nT, 'opts', o, ...
        'g', g, 'ch', ch, 'nerr', nerr, 'nbits', nbits, 'seconds', toc(t0), ...
        'floor', 1 ./ max(squeeze(sum(nbits, 1)), 1));
    print_acceptance(res);
    if ~exist(fullfile(root, 'results'), 'dir'), mkdir(fullfile(root, 'results')); end
    save(fullfile(root, 'results', sprintf('ber_%s_%s%s.mat', strrep(o.channel, '-', '_'), o.tier, o.tag)), '-struct', 'res');
    plot_ber_vs_snr(res);
end

function [wf, rx] = config_to_rx(name, o)
    rx = struct('E', o.E, 'mp_iter', o.mp_iter, 'mp_damp', o.mp_damp, 'mp_mode', o.mp_mode, 'n_dd', o.n_dd, ...
        'pilot', 'matched', 'pilot_boost_dB', o.pilot_boost_dB, 'use_ici', o.use_ici);
    switch lower(name)
        case 'ofdm-dftdd',       wf = 'cp-ofdm';  rx.ofdm_csi = 'ce';  rx.ce_mode = 'dftdd';
        case 'ofdm-dft',         wf = 'cp-ofdm';  rx.ofdm_csi = 'ce';  rx.ce_mode = 'dft';
        case 'ofdm-interp',      wf = 'cp-ofdm';  rx.ofdm_csi = 'ce';  rx.ce_mode = 'interp';
        case 'ofdm-ideal1tap',   wf = 'cp-ofdm';  rx.ofdm_csi = 'ideal1tap';
        case 'ofdm-idealbanded', wf = 'cp-ofdm';  rx.ofdm_csi = 'idealbanded';
        case 'otfs-ce',          wf = 'cp-otfs';  rx.otfs_csi = 'ce';
        case 'otfs-ideal',       wf = 'cp-otfs';  rx.otfs_csi = 'ideal';
        otherwise, error('run_ber_eva:config', 'unknown config %s', name);
    end
end

function print_acceptance(res)
    o = res.opts;  snr = res.snr_db;  lin = 10.^(snr/10);
    fprintf('\nBER table (rows SNR dB):\n  %-6s', 'SNR');  fprintf('%-18s', res.configs{:});  fprintf('\n');
    for si = 1:numel(snr), fprintf('  %-6g', snr(si));  fprintf('%-18.2e', res.ber(si, :));  fprintf('\n'); end
    gb = lin/2;                                                  % per-bit SNR for unit-power 4-QAM
    th_awgn = 0.5*erfc(sqrt(gb));                                 % Q(sqrt(2 gamma_b))
    th_ray  = 0.5*(1 - sqrt(gb ./ (1 + gb)));
    ci = @(n) find(strcmp(res.configs, n));
    fprintf('Acceptance (F2.2):\n');
    switch lower(o.channel)
        case 'awgn'
            for c = 1:numel(res.configs)
                m = res.ber(:, c) > 20*res.floor(:, c);          % only well-measured points
                if any(m), fprintf('  %-16s vs Q(sqrt(2 gamma_b)): max rel dev %.2f  %s\n', res.configs{c}, max(abs(res.ber(m, c)./th_awgn(m).' - 1)), pw(max(abs(res.ber(m, c)./th_awgn(m).' - 1)) < 0.10)); end
            end
        case 'flat_rayleigh'
            % One fade per frame: the ensemble average needs >> 1000 draws, so the
            % sharp check is PAIRED: per trial the conditional AWGN theory
            % Q(sqrt(2 gamma_b |alpha|^2)) at that trial's own gain (regenerated from
            % the same stream) vs the simulated errors, summed over trials.
            for c = 1:numel(res.configs)
                if ~any(strcmp(res.configs{c}, {'ofdm-ideal1tap', 'otfs-ideal'})), continue; end
                nerr_sim = zeros(numel(snr), 1);  nerr_th = zeros(numel(snr), 1);
                for si = 1:numel(snr)
                    for t = 1:res.nT
                        cr = chan_realize(res.ch, res.g, crn_stream(o.seed, 1000 + round(10*snr(si)), t));
                        gb_t = gb(si)*abs(cr.alpha)^2;
                        nerr_th(si) = nerr_th(si) + res.nbits(t, si, c)*0.5*erfc(sqrt(gb_t));
                        nerr_sim(si) = nerr_sim(si) + res.nerr(t, si, c);
                    end
                end
                m = nerr_sim >= 100;
                r = nerr_sim ./ max(nerr_th, eps);
                fprintf('  %-16s paired conditional Rayleigh check (sim errors / theory at each trial''s own gain): %s  %s\n', res.configs{c}, mat2str(r.', 3), pw(all(abs(r(m) - 1) < 0.10)));
            end
            c = ci('ofdm-ideal1tap');
            if ~isempty(c), fprintf('  OFDM ideal 1-tap vs the order-1 Rayleigh ENSEMBLE curve (needs >> 1000 fades): max rel dev %.2f (informational at %d trials)\n', max(abs(res.ber(:, c)./th_ray.' - 1)), res.nT); end
        otherwise
            c1 = ci('ofdm-dftdd');  c2 = ci('ofdm-ideal1tap');  c3 = ci('otfs-ce');  c4 = ci('ofdm-interp');
            if ~isempty(c1) && ~isempty(c2)
                m = snr >= 10;  d = max(abs(res.ber(m, c1)./res.ber(m, c2) - 1));
                fprintf('  dftdd vs ideal 1-tap above 10 dB: max rel dev %.2f (expect < 0.10; old: < 0.01)  %s\n', d, pw(d < 0.10));
            end
            if ~isempty(c2) && o.doppler_scale > 0
                fprintf('  OFDM ideal 1-tap at %g dB = %.2e (with Doppler an ICI floor ~4e-4 on reduced EVA is expected)\n', snr(end), res.ber(end, c2));
            end
            if ~isempty(c3) && ~isempty(c1)
                k = find(snr == 18);
                if ~isempty(k), fprintf('  OTFS-CE vs OFDM-dftdd at 18 dB: %.2e vs %.2e (expect >= 2 decades gap)  %s\n', res.ber(k, c3), res.ber(k, c1), pw(res.ber(k, c3) <= 1e-2*res.ber(k, c1) || res.ber(k, c3) == 0)); end
            end
            if ~isempty(c4) && ~isempty(c1)
                fprintf('  interp floor vs dftdd at %g dB: %.2e vs %.2e\n', snr(end), res.ber(end, c4), res.ber(end, c1));
            end
    end
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
