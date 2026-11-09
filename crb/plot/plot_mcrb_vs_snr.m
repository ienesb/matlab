function fig = plot_mcrb_vs_snr(res)
% PLOT_MCRB_VS_SNR  Delay and Doppler RMSE bounds vs SNR from the real detector
% (run_mcrb_rx_eva): LB (bold), MCRB0 (solid), genie CRB (dotted), BSC-matched LB
% (dashed, when computed), pilot-only CRB (dash-dot), per waveform; headline
% statistic of the driver (res.opts.statistic).
    root = fileparts(fileparts(mfilename('fullpath')));
    o = res.opts;  snr = res.snr_db(:);  nW = numel(res.waveforms);  cols = lines(nW);
    fig = new_fig([80 80 1000 420]);
    for blk = 1:2
        if blk == 1, tag = 'tau'; sc = 1e9; lab = 'delay RMSE bound (ns)'; else, tag = 'nu'; sc = 1; lab = 'Doppler RMSE bound (Hz)'; end
        subplot(1, 2, blk);  hold on;
        for iw = 1:nW
            c = cols(iw, :);  nm = res.waveforms{iw};
            semilogy(snr, sc*res.(['rmse_' tag '_lb'])(:, iw),   '-',  'Color', c, 'LineWidth', 2.2, 'DisplayName', [nm ' LB']);
            semilogy(snr, sc*res.(['rmse_' tag '_mcrb'])(:, iw), '-',  'Color', c, 'LineWidth', 1.1, 'DisplayName', [nm ' MCRB0']);
            semilogy(snr, sc*res.(['rmse_' tag '_crb'])(:, iw),  ':',  'Color', c, 'LineWidth', 1.2, 'DisplayName', [nm ' genie CRB']);
            y = sc*res.(['rmse_' tag '_lb_bsc'])(:, iw);
            if o.do_bsc && any(isfinite(y)), semilogy(snr, y, '--', 'Color', c, 'LineWidth', 1.0, 'DisplayName', sprintf('%s BSC-matched LB (%s)', nm, o.bsc_mode)); end
            if iw == 1, semilogy(snr, sc*res.(['rmse_' tag '_pil'])(:, iw), '-.', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.0, 'DisplayName', sprintf('pilot-only CRB (%.0f %% comb)', 100*o.pilot_frac)); end
        end
        set(gca, 'YScale', 'log');  grid on;  xlabel('SNR (dB)');  ylabel(lab);
        if blk == 1, legend('Location', 'southwest', 'FontSize', 7); end
    end
    sgtitle(sprintf('%s real-detector MCRB: %d trials, %s (drop <= %.0f%%), %s gains, ICI=%d, Tcp=%d Ts, MP %s', ...
        res.ch.name, res.nT, o.statistic, 100*max(res.drop_frac(:), [], 'omitnan'), res.ch.gain_model, o.use_ici, o.Tcp_samples, o.mp_mode), 'Interpreter', 'none', 'FontSize', 10);
    if ~exist(fullfile(root, 'figures'), 'dir'), mkdir(fullfile(root, 'figures')); end
    print(fig, fullfile(root, 'figures', sprintf('fig_mcrb_rx_%s_%s%s.png', lower(strrep(res.ch.name, '-', '_')), o.tier, o.tag)), '-dpng', '-r120');
end
