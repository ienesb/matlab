function fig = plot_crb_vs_snr(res, statistic)
% PLOT_CRB_VS_SNR  Delay and Doppler RMSE bounds vs SNR for every waveform /
% ICI case of a run_crb_eva (or run_crb_mt40_frame) result.
%   plot_crb_vs_snr(res, 'mean' | 'median')
% Panels: delay trace, delay path 1, Doppler trace, Doppler path 1 (RMSE = sqrt
% of the CRB), solid = ICI, dashed = no-ICI; black dotted = CT floor, black
% dash-dot = flat-spectrum reference (3x floor).  Legend shows the drop fraction.
    if nargin < 2, statistic = 'mean'; end
    root = fileparts(fileparts(mfilename('fullpath')));
    snr = res.snr_db(:);  snr_lin = 10.^(snr/10);
    nW = numel(res.waveforms);  nI = numel(res.ici_list);
    cols = lines(nW);
    fig = new_fig([80 80 1200 800]);
    panels = {'crb_tau', 'trace', 'delay RMSE bound (ns)', 1e9; 'crb_tau', 'path1', 'delay RMSE bound, path 1 (ns)', 1e9; ...
              'crb_nu',  'trace', 'Doppler RMSE bound (Hz)', 1;  'crb_nu',  'path1', 'Doppler RMSE bound, path 1 (Hz)', 1};
    for k = 1:4
        subplot(2, 2, k);  hold on;
        A = res.(panels{k,1});  fl = res.(strrep(panels{k,1}, 'crb', 'fl'));
        for iw = 1:nW
            for ii = 1:nI
                v = res.valid(:, iw, ii);
                if strcmp(panels{k,2}, 'trace'), x = sum(A(v, :, iw, ii), 2); else, x = A(v, 1, iw, ii); end
                if strcmp(statistic, 'median'), c = median(x); else, c = mean(x); end
                ls = '-';  if ~res.ici_list(ii), ls = '--'; end
                semilogy(snr, panels{k,4}*sqrt(c ./ snr_lin), ls, 'Color', cols(iw,:), 'LineWidth', 1.4, ...
                    'DisplayName', sprintf('%s %s (drop %.0f%%)', res.waveforms{iw}, ternary(res.ici_list(ii), 'ICI', 'no-ICI'), 100*res.drop_frac(iw, ii)));
            end
        end
        if strcmp(panels{k,2}, 'trace'), f0 = mean(sum(fl, 2)); else, f0 = mean(fl(:, 1)); end
        semilogy(snr, panels{k,4}*sqrt(f0 ./ snr_lin), 'k:', 'LineWidth', 1.2, 'DisplayName', 'CT floor');
        semilogy(snr, panels{k,4}*sqrt(3*f0 ./ snr_lin), 'k-.', 'LineWidth', 1.0, 'DisplayName', 'flat (3x floor)');
        set(gca, 'YScale', 'log');  grid on;  xlabel('SNR (dB)');  ylabel(panels{k,3});
        if k == 1, legend('Location', 'southwest', 'FontSize', 7); end
    end
    sgtitle(sprintf('%s: %s gains, Tcp = %d Ts, %d trials, %s over valid trials, %s', res.ch.name, res.opts.gain_model, ...
        res.opts.Tcp_samples, res.nT, statistic, res.range_convention), 'Interpreter', 'none');
    if ~exist(fullfile(root, 'figures'), 'dir'), mkdir(fullfile(root, 'figures')); end
    print(fig, fullfile(root, 'figures', sprintf('fig_crb_%s_%s.png', lower(strrep(res.ch.name, '-', '_')), res.opts.tier)), '-dpng', '-r120');
end

function out = ternary(c, a, b)
    if c, out = a; else, out = b; end
end
