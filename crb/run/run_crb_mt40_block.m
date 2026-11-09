function res = run_crb_mt40_block(varargin)
% RUN_CRB_MT40_BLOCK  Liu et al. Fig. 1 on the MT40-1D framing (plan F1.1):
% ranging CRB [m^2] (summed over the 40 targets) vs SNR for the SC / AFDM /
% OTFS / OFDM unitary bases and the Jensen bound, 16-QAM and 16-PSK panels.
%
%   res = run_crb_mt40_block('tier', 'dev', 'snr_db', -10:1:6, ...)
% Options: tier ('smoke' 20 | 'dev' 500 | 'paper' 5000 data draws), snr_db,
% mods ({'16qam','16psk'}), waveforms, kappa_thresh (1e4), seed (7000), L (40).
% The FIM is built once per draw at sigma2 = 1 and scaled by 1/SNR (exact).
% Saves results/crb_mt40_block.mat; plot with plot_liu_fig1.
    o = struct('tier', 'dev', 'snr_db', -10:1:6, 'mods', {{'16qam', '16psk'}}, ...
        'waveforms', {{'SC', 'AFDM', 'OTFS', 'OFDM'}}, 'kappa_thresh', Inf, 'seed', 7000, 'L', 40);
    % kappa_thresh: with 3L = 120 parameters in N = 128 samples the bin-unit FIM
    % has cond > 1e4 on EVERY draw (identifiability edge), so no Omega_N cut is
    % applied by default at p = 0; the condition numbers are stored and printed.
    for i = 1:2:numel(varargin), o.(varargin{i}) = varargin{i+1}; end
    nT = tier_trials(o.tier, [20 500 5000]);
    root = fileparts(fileparts(mfilename('fullpath')));

    N = 128;  B = 160e6;
    g  = make_grid(N, 1, B/N, 0, 'cp');
    ch = chan_mt40(g, 'L', o.L);  cr = chan_realize(ch, g, []);
    Qh = liu_precoders(N, 32, 4, 1/16, 1/8);
    m2_per_samp2 = ch.range_unit_m^2;           % 0.87890625 m^2 per sample^2
    snr_lin = 10.^(o.snr_db/10);

    % Jensen bound: |X|^2 = 1 deterministic -> waveform/constellation independent
    Jb = fim(ones(N,1), 'cp-ofdm', g, cr, false, 1);
    Sb = crb_from_fim(Jb, g, cr);
    res.jensen_tr_samp2 = Sb.tr_tau / g.Ts^2;      % trace at sigma2 = 1, in samples^2
    res.jensen_m2 = res.jensen_tr_samp2 * m2_per_samp2 ./ snr_lin;

    nW = numel(o.waveforms);  nMod = numel(o.mods);
    tr = nan(nT, nW, nMod);  kap = tr;
    t0 = tic;
    for im = 1:nMod
        mdm = modem(o.mods{im});
        for iw = 1:nW
            Q = Qh.(o.waveforms{iw});
            trw = nan(nT, 1);  kw = nan(nT, 1);
            parfor t = 1:nT
                s = mdm.random(N, crn_stream(o.seed, im, t)).';        % CRN: same s for every waveform
                X = Q * s;
                J = fim(X, 'cp-ofdm', g, cr, false, 1);
                S = crb_from_fim(J, g, cr);
                trw(t) = S.tr_tau / g.Ts^2;  kw(t) = S.kappa;
            end
            tr(:, iw, im) = trw;  kap(:, iw, im) = kw;
            fprintf('  %-6s %-6s  median tr %.4e samp^2  mean %.4e  kappa median %.1e max %.1e  drop(kappa>%g) %.1f%%   [%.0f s]\n', ...
                o.mods{im}, o.waveforms{iw}, median(trw), mean(trw(kw <= o.kappa_thresh)), median(kw), max(kw), o.kappa_thresh, 100*mean(kw > o.kappa_thresh), toc(t0));
        end
    end
    valid = kap <= o.kappa_thresh;
    res.tr_samp2 = tr;  res.kappa = kap;  res.valid = valid;
    res.mean_m2   = zeros(numel(o.snr_db), nW, nMod);  res.median_m2 = res.mean_m2;
    for im = 1:nMod, for iw = 1:nW
        v = valid(:, iw, im);
        res.mean_m2(:, iw, im)   = mean(tr(v, iw, im))   * m2_per_samp2 ./ snr_lin(:);
        res.median_m2(:, iw, im) = median(tr(v, iw, im)) * m2_per_samp2 ./ snr_lin(:);
    end, end
    res.drop_frac = squeeze(1 - mean(valid, 1));
    res.opts = o;  res.nT = nT;  res.snr_db = o.snr_db;  res.g = g;  res.ch = ch;
    res.waveforms = o.waveforms;  res.mods = o.mods;  res.statistic = 'mean over valid draws (Liu); median stored';

    % acceptance table
    fprintf('\nAcceptance (F1.1):\n');
    io = find(strcmp(o.waveforms, 'OFDM'));  ip = find(strcmp(o.mods, '16psk'));  iq = find(strcmp(o.mods, '16qam'));
    if ~isempty(io) && ~isempty(ip)
        e = max(abs(tr(:, io, ip) - res.jensen_tr_samp2)) / res.jensen_tr_samp2;
        fprintf('  OFDM+16-PSK per-draw CRB == Jensen: max rel %.2e  %s\n', e, pw(e < 1e-10));
    end
    k10 = find(o.snr_db == -10);
    if ~isempty(k10)
        for im = 1:nMod
            fprintf('  %s at -10 dB [m^2] (mean | median):', o.mods{im});
            for iw = 1:nW, fprintf('  %s %.3f|%.3f', o.waveforms{iw}, res.mean_m2(k10, iw, im), res.median_m2(k10, iw, im)); end
            fprintf('  Jensen %.3f\n', res.jensen_m2(k10));
        end
        fprintf('  (paper read-off at -10 dB: cluster ~6.2, OFDM ~4.85 (16-QAM), Jensen ~3.9-4.0; absolute levels depend on the\n');
        fprintf('   random delay geometry and on beta = 1, so only the ratios below are asserted)\n');
        if ~isempty(io) && ~isempty(iq)
            others = setdiff(1:nW, io);
            gain_db = 10*log10(mean(res.mean_m2(k10, others, iq)) / res.mean_m2(k10, io, iq));
            fprintf('  16-QAM: OFDM gain vs the SC/AFDM/OTFS cluster = %.2f dB (paper: ~1 dB)  %s\n', gain_db, pw(gain_db > 0.5 && gain_db < 1.8));
        end
        if ~isempty(io) && ~isempty(ip)
            others = setdiff(1:nW, io);
            gain_db = 10*log10(mean(res.mean_m2(k10, others, ip)) / res.mean_m2(k10, io, ip));
            fprintf('  16-PSK: OFDM gain vs the cluster = %.2f dB (paper: ~2 dB)  %s\n', gain_db, pw(gain_db > 1.3 && gain_db < 2.7));
        end
    end
    if ~exist(fullfile(root, 'results'), 'dir'), mkdir(fullfile(root, 'results')); end
    save(fullfile(root, 'results', 'crb_mt40_block.mat'), '-struct', 'res');
    plot_liu_fig1(res);
end

function s = pw(ok)
    if ok, s = 'PASS'; else, s = 'WARN'; end
end

function nT = tier_trials(tier, tab)
    switch lower(tier)
        case 'smoke', nT = tab(1);
        case 'dev',   nT = tab(2);
        case 'paper', nT = tab(3);
        otherwise, nT = tier;   % numeric
    end
end
