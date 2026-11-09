function res = run_crb_mt40_frame(varargin)
% RUN_CRB_MT40_FRAME  Genie CRB vs SNR on the MT40-2D framing (plan F1.3):
% N = 128 x M = 16, per-symbol CP Tcp = 127 Ts (CP assumption satisfied),
% delay-only targets (nu block masked), every waveform of the derivations.
%
%   res = run_crb_mt40_frame('tier', 'dev', 'geometry', 'both', ...)
% Options (defaults): tier 'smoke'(20)|'dev'(200)|'paper'(2000) data draws;
% M 16; Tcp_samples 127; mods {'16psk','16qam','4qam'}; geometry
% 'cont' (continuous delays; the sampled RCP-OTFS branch is skipped) |
% 'int' (MT40-int, all 5 waveforms) | 'both'; L 40; doppler 'none' | 'mt40d'
% (MT40-D adds ICI and Doppler); snr_db -10:1:6; seed 7000; kappa_thresh Inf.
% At nu = 0 the ICI and no-ICI models coincide on the estimated parameters,
% so only the no-ICI case is run unless doppler = 'mt40d'.
    o = struct('tier', 'smoke', 'M', 16, 'Tcp_samples', 127, 'mods', {{'16psk', '16qam', '4qam'}}, ...
        'geometry', 'both', 'L', 40, 'doppler', 'none', 'snr_db', -10:1:6, 'seed', 7000, 'kappa_thresh', Inf);
    for i = 1:2:numel(varargin), o.(varargin{i}) = varargin{i+1}; end
    root = fileparts(fileparts(mfilename('fullpath')));
    nT = tier_trials(o.tier, [20 200 2000]);
    N = 128;  B = 160e6;  df = B/N;  Ts = 1/B;
    gcp  = make_grid(N, o.M, df, o.Tcp_samples*Ts, 'cp');
    grcp = make_grid(N, o.M, df, o.Tcp_samples*Ts, 'rcp');
    geoms = {};  if any(strcmp(o.geometry, {'cont', 'both'})), geoms{end+1} = 'cont'; end
    if any(strcmp(o.geometry, {'int', 'both'})), geoms{end+1} = 'int'; end
    if strcmp(o.doppler, 'mt40d'), ici_list = [true false]; else, ici_list = false; end
    res = struct();
    for ig = 1:numel(geoms)
        ch = chan_mt40(gcp, 'L', o.L, 'integer', strcmp(geoms{ig}, 'int'), 'doppler', o.doppler);
        if strcmp(geoms{ig}, 'cont')
            wfs = {'cp-ofdm', 'cp-otfs', 'rcp-ofdm', 'rcp-otfs-cont'};
        else
            wfs = {'cp-ofdm', 'cp-otfs', 'rcp-ofdm', 'rcp-otfs', 'rcp-otfs-cont'};
        end
        for im = 1:numel(o.mods)
            oc = struct('nT', nT, 'waveforms', {wfs}, 'ici_list', ici_list, 'modtype', o.mods{im}, ...
                'seed', o.seed + im, 'kappa_thresh', o.kappa_thresh);
            fprintf('MT40 frame: %s geometry, %s, %d draws\n', geoms{ig}, o.mods{im}, nT);
            r = crb_frame_core(ch, gcp, grcp, oc);
            r.snr_db = o.snr_db;  r.opts = o;  r.opts.tier = sprintf('%s_%s_%s', o.tier, geoms{ig}, o.mods{im});
            r.opts.gain_model = ch.gain_model;  r.opts.Tcp_samples = o.Tcp_samples;
            r.range_convention = 'monostatic c/(2B)';
            r.tr_tau_samp2 = squeeze(mean(sum(r.crb_tau, 2), 1)) / Ts^2;    % trace, samples^2, sigma2 = 1
            r.med_tau_samp2 = squeeze(median(sum(r.crb_tau, 2), 1)) / Ts^2;
            fprintf('  delay trace CRB at sigma2 = 1 [samples^2] (mean | median):');
            for iw = 1:numel(wfs), fprintf('  %s %.4f|%.4f', wfs{iw}, r.tr_tau_samp2(iw), r.med_tau_samp2(iw)); end
            fprintf('\n');
            res.(sprintf('%s_%s', geoms{ig}, o.mods{im})) = r;
            plot_crb_vs_snr(r, 'mean');
        end
    end
    res.opts = o;  res.nT = nT;
    if ~exist(fullfile(root, 'results'), 'dir'), mkdir(fullfile(root, 'results')); end
    save(fullfile(root, 'results', sprintf('crb_mt40_frame_%s.mat', o.tier)), '-struct', 'res');
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
