function res = run_crb_eva(varargin)
% RUN_CRB_EVA  Genie CRB vs SNR on the EVA channel for every waveform / ICI
% case of the derivations (plan F1.2):
%   {cp-ofdm, cp-otfs, rcp-ofdm, rcp-otfs (sampled), rcp-otfs-cont} x {ICI, no-ICI}
% Monte Carlo over channel draws (Rayleigh gains, Jakes Doppler) and data.
%
%   res = run_crb_eva('tier', 'smoke', 'reduced', false, ...)
% Options (defaults): tier 'smoke'(20)|'dev'(200)|'paper'(2000) trials;
% snr_db 0:1:25; reduced (auto from results/diag_cond_eva.mat, else false);
% gain_model 'rayleigh'; Tcp_samples 20; doppler_scale 1; snap true;
% modtype '4qam'; waveforms; ici_list [true false]; seed 9000; kappa_thresh 1e10.
% Every FIM is built at sigma2 = 1 and scaled by 1/SNR (exact).  Saves
% results/crb_eva_<tier>.mat and plots with plot_crb_vs_snr.
    o = struct('tier', 'smoke', 'snr_db', 0:1:25, 'reduced', [], 'gain_model', 'rayleigh', ...
        'Tcp_samples', 20, 'doppler_scale', 1, 'snap', true, 'modtype', '4qam', ...
        'waveforms', {{'cp-ofdm', 'cp-otfs', 'rcp-ofdm', 'rcp-otfs', 'rcp-otfs-cont'}}, ...
        'ici_list', [true false], 'seed', 9000, 'kappa_thresh', 1e10, 'tag', '', 'variant', '');
    o = parse_opts(o, varargin, 'run_crb_eva');
    root = fileparts(fileparts(mfilename('fullpath')));
    if isempty(o.reduced)
        f = fullfile(root, 'results', 'diag_cond_eva.mat');
        o.reduced = false;
        if exist(f, 'file')
            d = load(f, 'decision');
            o.reduced = contains(d.decision, 'reduced');
            fprintf('conditioning gate decision: %s\n', d.decision);
        else
            fprintf('no diag_cond_eva.mat found -> full EVA (run diag_cond_eva first)\n');
        end
    end
    o.nT = tier_trials(o.tier, [20 200 2000]);
    N = 512; M = 128; df = 15e3;  Ts = 1/(N*df);
    gcp  = make_grid(N, M, df, o.Tcp_samples*Ts, 'cp');
    grcp = make_grid(N, M, df, o.Tcp_samples*Ts, 'rcp');
    ch = chan_eva(gcp, 'reduced', o.reduced, 'gain_model', o.gain_model, 'snap', o.snap, ...
        'doppler_scale', o.doppler_scale, 'variant', o.variant);
    fprintf('EVA CRB run: %s, P = %d, %s gains, Tcp = %d Ts, %d trials, %d waveforms x %d ICI cases\n', ...
        ch.name, ch.P, o.gain_model, o.Tcp_samples, o.nT, numel(o.waveforms), numel(o.ici_list));
    res = crb_frame_core(ch, gcp, grcp, o);
    res.snr_db = o.snr_db;  res.opts = o;  res.range_convention = 'bistatic c*tau';
    print_acceptance(res);
    if ~exist(fullfile(root, 'results'), 'dir'), mkdir(fullfile(root, 'results')); end
    % file name carries the channel variant: crb_eva_<tier>.mat / crb_eva_reduced_<tier>.mat
    save(fullfile(root, 'results', sprintf('crb_%s_%s%s.mat', lower(strrep(ch.name, '-', '_')), o.tier, o.tag)), '-struct', 'res');
    plot_crb_vs_snr(res);
end

function print_acceptance(res)
    o = res.opts;  nW = numel(o.waveforms);  nI = numel(o.ici_list);
    fprintf('\nAcceptance (F1.2), traces at sigma2 = 1, mean over valid trials:\n');
    tr = @(A, iw, ii) mean(sum(A(res.valid(:, iw, ii), :, iw, ii), 2), 1);
    for iw = 1:nW
        fprintf('  %-14s', o.waveforms{iw});
        for ii = 1:nI
            fprintf('  ici=%d: tau %.3e s^2  nu %.3e Hz^2  kappa med %.1e  drop %.1f%%', o.ici_list(ii), ...
                tr(res.crb_tau, iw, ii), tr(res.crb_nu, iw, ii), median(res.kappa(:, iw, ii)), 100*res.drop_frac(iw, ii));
        end
        fprintf('\n');
    end
    io = find(strcmp(o.waveforms, 'cp-ofdm'));  ic1 = find(o.ici_list);  ic0 = find(~o.ici_list);
    if ~isempty(io) && ~isempty(ic1) && ~isempty(ic0)
        d = abs(tr(res.crb_tau, io, ic1) - tr(res.crb_tau, io, ic0)) / tr(res.crb_tau, io, ic0);
        fprintf('  CP-OFDM delay CRB ICI vs no-ICI: %.2f%% (expect < 0.5%%)  %s\n', 100*d, pw(d < 0.005));
        d = abs(tr(res.crb_nu, io, ic1) - tr(res.crb_nu, io, ic0)) / tr(res.crb_nu, io, ic0);
        fprintf('  CP-OFDM Doppler CRB ICI vs no-ICI: %.2f%% (expect < 2%%)  %s\n', 100*d, pw(d < 0.02));
        e = mean(sum(res.exp_tau(:, :, io, ic0), 2)) ;
        fprintf('  CP-OFDM no-ICI: mean of J^-1 vs (E J)^-1 delay trace: %.2f%% (Jensen gap; 0 for constant-modulus data)\n', ...
            100*abs(tr(res.crb_tau, io, ic0) - e)/e);
    end
    it = find(strcmp(o.waveforms, 'cp-otfs'));
    if ~isempty(io) && ~isempty(it) && ~isempty(ic1)
        d = abs(tr(res.crb_tau, it, ic1) - tr(res.crb_tau, io, ic1)) / tr(res.crb_tau, io, ic1);
        fprintf('  CP-OTFS vs CP-OFDM delay CRB (E[J^-1], differs only through data randomness): %.2f%%\n', 100*d);
    end
    ir = find(strcmp(o.waveforms, 'rcp-ofdm'));
    if ~isempty(io) && ~isempty(ir) && ~isempty(ic1)
        d = abs(tr(res.crb_tau, ir, ic1) - tr(res.crb_tau, io, ic1)) / tr(res.crb_tau, io, ic1);
        fprintf('  RCP-OFDM vs CP-OFDM delay CRB departure: %.2f%% (expect O(1/M)+O(tau/T), < 10%%)  %s\n', 100*d, pw(d < 0.10));
    end
    is = find(strcmp(o.waveforms, 'rcp-otfs'));  ic = find(strcmp(o.waveforms, 'rcp-otfs-cont'));
    if ~isempty(is) && ~isempty(ic) && ~isempty(ic1)
        d = abs(tr(res.crb_nu, is, ic1) - tr(res.crb_nu, ic, ic1)) / tr(res.crb_nu, ic, ic1);
        fprintf('  RCP-OTFS sampled vs continuous Doppler CRB gap: %.2f%% (different models, reported only)\n', 100*d);
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
