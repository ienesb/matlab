function fig = plot_liu_fig1(res)
% PLOT_LIU_FIG1  Ranging CRB [m^2] vs SNR, one panel per constellation
% (Liu et al. Fig. 1 layout: linear y axis, SC/AFDM/OTFS/OFDM/Jensen).
    if nargin < 1
        root = fileparts(fileparts(mfilename('fullpath')));
        res = load(fullfile(root, 'results', 'crb_mt40_block.mat'));
    end
    root = fileparts(fileparts(mfilename('fullpath')));
    nMod = numel(res.mods);  nW = numel(res.waveforms);
    sty = struct('SC', '-', 'AFDM', '--', 'OTFS', '-.', 'OFDM', ':');
    fig = new_fig([100 100 1100 420]);
    for im = 1:nMod
        subplot(1, nMod, im);  hold on;
        for iw = 1:nW
            plot(res.snr_db, res.mean_m2(:, iw, im), sty.(res.waveforms{iw}), 'LineWidth', 1.6, ...
                'DisplayName', sprintf('%s (drop %.1f%%)', res.waveforms{iw}, 100*res.drop_frac(iw, im)));
        end
        if strcmpi(res.mods{im}, '16psk')
            plot(res.snr_db, res.jensen_m2, 'kx', 'MarkerSize', 7, 'DisplayName', 'Jensen bound');
        else
            plot(res.snr_db, res.jensen_m2, 'k-', 'LineWidth', 1.2, 'DisplayName', 'Jensen bound');
        end
        grid on;  xlabel('SNR (dB)');  ylabel('CRB (m^2), sum over 40 targets');
        title(sprintf('%s  (N=128, L=%d, %d draws, mean over valid, monostatic c/(2B))', upper(res.mods{im}), res.ch.P, res.nT));
        legend('Location', 'northeast');
        ytop = max([res.mean_m2(1, :, im), res.jensen_m2(1)], [], 'omitnan');
        if isfinite(ytop) && ytop > 0, ylim([0, ytop*1.1]); end
    end
    if ~exist(fullfile(root, 'figures'), 'dir'), mkdir(fullfile(root, 'figures')); end
    print(fig, fullfile(root, 'figures', 'fig_liu_fig1.png'), '-dpng', '-r130');
end
