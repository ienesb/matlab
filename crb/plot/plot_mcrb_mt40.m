function fig = plot_mcrb_mt40(res)
% PLOT_MCRB_MT40  MT40 ranging bounds [m^2, sum over targets]: vs p at 0 dB and
% vs SNR at the largest p, per constellation; CRB, MCRB0, LB, MCRB@eta_bar,
% Jensen and pilot-only references.  Median over valid draws of the PER-DRAW
% composite bounds (mcrb_stats).  The pilot-only line is omitted when the comb
% cannot identify the targets (res.pilot_tr_tau = NaN).  Works for the 1-D
% (Liu precoders) and 2-D framings (res.framing).
    root = fileparts(fileparts(mfilename('fullpath')));
    p = res.p_list;  Ts = res.g.Ts;  c = res.range_m2_per_samp2;  nW = numel(res.waveforms);  nMod = numel(res.mods);
    snr = res.snr_db(:);  lin = 10.^(snr/10);  nS = numel(snr);
    cols = lines(nW);
    mk = '';  if numel(p) == 1, mk = 'o'; end                      % single p: show markers
    fig = new_fig([60 60 480*nMod 760]);
    for im = 1:nMod
        if strcmpi(res.mods{im}, '16qam'), xl = 'symbol error rate p'; else, xl = 'Gray-bit flip rate p (BSC)'; end
        subplot(2, nMod, im);  hold on;
        for iw = 1:nW
            v = res.valid(:, :, iw, im);
            S = mcrb_stats(res.G_tau(:, :, iw, im), res.F_tau(:, :, iw, im), res.BB_tau(:, :, iw, im), res.Gc_tau(:, :, iw, im), v, 1, 'median');
            T = mcrb_stats(res.Gt_tau(:, :, iw, im), res.Ft_tau(:, :, iw, im), 0*res.BB_tau(:, :, iw, im), res.Gc_tau(:, :, iw, im), v, 1, 'median');
            semilogy(p, S.crb/Ts^2*c, [':' mk], 'Color', cols(iw,:), 'DisplayName', [res.waveforms{iw} ' CRB']);
            semilogy(p, S.mcrb/Ts^2*c, ['-' mk], 'Color', cols(iw,:), 'LineWidth', 1.2, 'DisplayName', [res.waveforms{iw} ' MCRB0']);
            semilogy(p, S.lb/Ts^2*c, ['-' mk], 'Color', cols(iw,:), 'LineWidth', 2.2, 'DisplayName', [res.waveforms{iw} ' LB']);
            semilogy(p, T.mcrb/Ts^2*c, ['--' mk], 'Color', cols(iw,:), 'DisplayName', [res.waveforms{iw} ' MCRB@eta\_bar']);
        end
        yline(res.jensen_tr_tau/Ts^2*c, 'k-', 'DisplayName', 'Jensen (all-data)');
        if isfinite(res.pilot_tr_tau), yline(res.pilot_tr_tau/Ts^2*c, 'k-.', 'DisplayName', sprintf('pilot-only CRB (%.0f %% comb)', 100*res.opts.pilot_frac)); end
        set(gca, 'YScale', 'log');  grid on;  xlabel(xl);  ylabel('ranging bound (m^2), SNR = 0 dB');
        title(sprintf('%s, L = %d', upper(res.mods{im}), res.cr.P));  if im == 1, legend('Location', 'northwest', 'FontSize', 6); end
        subplot(2, nMod, nMod + im);  hold on;  ip = numel(p);
        for iw = 1:nW
            v = res.valid(:, ip, iw, im);
            y_mcrb = zeros(nS, 1);  y_lb = y_mcrb;  y_crb = y_mcrb;
            for k = 1:nS
                S = mcrb_stats(res.G_tau(:, ip, iw, im), res.F_tau(:, ip, iw, im), res.BB_tau(:, ip, iw, im), res.Gc_tau(:, ip, iw, im), v, 1/lin(k), 'median');
                y_mcrb(k) = S.mcrb;  y_lb(k) = S.lb;  y_crb(k) = S.crb;
            end
            semilogy(snr, y_crb/Ts^2*c, ':', 'Color', cols(iw,:), 'DisplayName', [res.waveforms{iw} ' CRB']);
            semilogy(snr, y_mcrb/Ts^2*c, '-', 'Color', cols(iw,:), 'LineWidth', 1.2, 'DisplayName', [res.waveforms{iw} ' MCRB0']);
            semilogy(snr, y_lb/Ts^2*c, '-', 'Color', cols(iw,:), 'LineWidth', 2.2, 'DisplayName', [res.waveforms{iw} ' LB']);
        end
        semilogy(snr, (res.jensen_tr_tau./lin)/Ts^2*c, 'k-', 'DisplayName', 'Jensen (all-data)');
        if isfinite(res.pilot_tr_tau), semilogy(snr, (res.pilot_tr_tau./lin)/Ts^2*c, 'k-.', 'DisplayName', 'pilot-only CRB'); end
        set(gca, 'YScale', 'log');  grid on;  xlabel('SNR (dB)');  ylabel(sprintf('ranging bound (m^2), %s = %g', xl, p(ip)));
        if im == 1, legend('Location', 'southwest', 'FontSize', 6); end
    end
    sgtitle(sprintf('%s %dx%d, Tcp = %d Ts, %d draws, median over valid (kappa <= %g), c/(2B); PSK: bit BSC, QAM: symbol errors', ...
        res.framing, res.g.N, res.g.M, res.opts.Tcp_samples, res.nT, res.opts.kappa_thresh), 'FontSize', 10);
    if ~exist(fullfile(root, 'figures'), 'dir'), mkdir(fullfile(root, 'figures')); end
    print(fig, fullfile(root, 'figures', sprintf('fig_mcrb_mt40_%s_L%d_%s%s.png', lower(strrep(res.framing, 'MT40-', '')), res.cr.P, res.opts.tier, res.opts.tag)), '-dpng', '-r120');
end
