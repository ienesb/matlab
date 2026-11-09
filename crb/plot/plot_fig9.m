function fig = plot_fig9(res)
% PLOT_FIG9  BER vs SNR_d in the layout of Raviteja TVT'19 Fig. 9, with the
% impl/otfs_ber reference points overlaid.
    root = fileparts(fileparts(mfilename('fullpath')));
    fig = new_fig([100 100 640 480]);
    sty = {'k--', 'b-o', 'r-s', 'k-^', 'k-v'};
    hold on;
    for c = 1:numel(res.curves)
        semilogy(res.snr_d_db, res.ber(:, c), sty{min(c, numel(sty))}, 'LineWidth', 1.4, 'DisplayName', res.curves{c});
    end
    r = res.ref;
    semilogy(r.snr, r.otfs_ideal, 'k:', 'DisplayName', 'ref OTFS ideal (impl)');
    semilogy(r.snr, r.ofdm, 'k-.', 'DisplayName', 'ref OFDM ideal (impl)');
    set(gca, 'YScale', 'log');  grid on;  xlabel('SNR_d (dB)');  ylabel('BER');
    title(sprintf('%s, N=%d M=%d, l_\\tau=%d, 4-QAM, T=3\\sigma, %d trials', res.ch.name, res.g.N, res.g.M, res.ch.l_tau, res.nT));
    legend('Location', 'southwest', 'FontSize', 8);
    ylim([1e-5 1e-1]);
    if ~exist(fullfile(root, 'figures'), 'dir'), mkdir(fullfile(root, 'figures')); end
    print(fig, fullfile(root, 'figures', sprintf('fig_ber_fig9_%s.png', res.opts.tier)), '-dpng', '-r130');
end
