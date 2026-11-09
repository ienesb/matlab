function fig = plot_ber_vs_snr(res)
% PLOT_BER_VS_SNR  BER vs SNR for every receiver configuration of run_ber_eva,
% with the Gray 4-QAM AWGN and order-1 Rayleigh references and the measurement
% floor 1/(trials * bits).  Zero-BER points are masked.
    root = fileparts(fileparts(mfilename('fullpath')));
    snr = res.snr_db(:);  gb = 10.^(snr/10)/2;
    fig = new_fig([100 100 760 540]);  hold on;
    sty = {'-o', '-s', '-^', '--', '-.', '-d', '--v', '-x'};
    for c = 1:numel(res.configs)
        b = res.ber(:, c);  b(b == 0) = NaN;
        semilogy(snr, b, sty{mod(c-1, numel(sty))+1}, 'LineWidth', 1.3, 'MarkerSize', 4, 'DisplayName', res.configs{c});
    end
    semilogy(snr, 0.5*erfc(sqrt(gb)), 'k:', 'DisplayName', 'AWGN Q(\surd(2\gamma_b))');
    semilogy(snr, 0.5*(1 - sqrt(gb./(1+gb))), 'k-.', 'DisplayName', 'Rayleigh order 1');
    fl = max(res.floor, [], 2);
    semilogy(snr, fl, 'Color', [0.6 0.6 0.6], 'LineStyle', ':', 'DisplayName', 'measurement floor');
    set(gca, 'YScale', 'log');  grid on;  xlabel('SNR (dB), N_0 = 1/SNR');  ylabel('BER');
    ylim([1e-7 1]);
    title(sprintf('%s: %s gains, Doppler x%g, ICI=%d, Tcp=%d Ts, %d trials', res.ch.name, res.ch.gain_model, ...
        res.ch.doppler_scale, res.opts.use_ici, res.opts.Tcp_samples, res.nT), 'Interpreter', 'none');
    legend('Location', 'southwest', 'FontSize', 8);
    if ~exist(fullfile(root, 'figures'), 'dir'), mkdir(fullfile(root, 'figures')); end
    tag = '';  if isfield(res.opts, 'tag'), tag = res.opts.tag; end
    print(fig, fullfile(root, 'figures', sprintf('fig_ber_%s_%s%s.png', strrep(res.opts.channel, '-', '_'), res.opts.tier, tag)), '-dpng', '-r130');
end
