function fig = plot_ber_mcrb_panels(res)
% PLOT_BER_MCRB_PANELS  Three stacked panels sharing the SNR axis (old
% plot_ber_mcrb layout): detection BER/SER, delay LB/MCRB0/CRB, Doppler
% LB/MCRB0/CRB -- "BER and MCRB are not the same currency" at a glance.
    root = fileparts(fileparts(mfilename('fullpath')));
    o = res.opts;  snr = res.snr_db(:);  nW = numel(res.waveforms);  cols = lines(nW);
    fig = new_fig([80 80 560 900]);
    subplot(3, 1, 1);  hold on;
    for iw = 1:nW
        b = res.ber(:, iw);  b(b == 0) = NaN;  s = res.ser(:, iw);  s(s == 0) = NaN;
        semilogy(snr, b, '-o', 'Color', cols(iw, :), 'DisplayName', [res.waveforms{iw} ' BER']);
        semilogy(snr, s, '--', 'Color', cols(iw, :), 'DisplayName', [res.waveforms{iw} ' SER']);
    end
    set(gca, 'YScale', 'log');  grid on;  ylabel('error rate');  legend('Location', 'southwest', 'FontSize', 7);
    title(sprintf('%s, %d trials, %s statistic', res.ch.name, res.nT, o.statistic), 'Interpreter', 'none');
    for blk = 1:2
        if blk == 1, tag = 'tau'; sc = 1e9; lab = 'delay RMSE (ns)'; else, tag = 'nu'; sc = 1; lab = 'Doppler RMSE (Hz)'; end
        subplot(3, 1, 1 + blk);  hold on;
        for iw = 1:nW
            c = cols(iw, :);
            semilogy(snr, sc*res.(['rmse_' tag '_lb'])(:, iw),   '-', 'Color', c, 'LineWidth', 2, 'DisplayName', [res.waveforms{iw} ' LB']);
            semilogy(snr, sc*res.(['rmse_' tag '_mcrb'])(:, iw), '-', 'Color', c, 'LineWidth', 1, 'DisplayName', [res.waveforms{iw} ' MCRB0']);
            semilogy(snr, sc*res.(['rmse_' tag '_crb'])(:, iw),  ':', 'Color', c, 'LineWidth', 1.2, 'DisplayName', [res.waveforms{iw} ' CRB']);
        end
        set(gca, 'YScale', 'log');  grid on;  ylabel(lab);
        if blk == 2, xlabel('SNR (dB)'); end
    end
    if ~exist(fullfile(root, 'figures'), 'dir'), mkdir(fullfile(root, 'figures')); end
    print(fig, fullfile(root, 'figures', sprintf('fig_ber_mcrb_panels_%s_%s%s.png', lower(strrep(res.ch.name, '-', '_')), o.tier, o.tag)), '-dpng', '-r120');
end
