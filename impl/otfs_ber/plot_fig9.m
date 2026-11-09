function plot_fig9(p, BER_ideal, BER_ce, BER_ofdm)
% PLOT_FIG9  Render the reproduction of Fig. 9 (BER vs SNR_d, fractional
% Doppler, full guard symbols) and save it as fig_9.png.

x = p.SNRd_dB;
f = figure('Color','w','Position',[100 100 640 480]);
try f.Theme = 'light'; catch, end              % force light theme (R2025+)
mk = {'-o','-s','-^'};
co = lines(numel(p.SNRp_dB));
ax = axes('Parent',f,'Color','w','XColor','k','YColor','k');
hold(ax,'on'); grid(ax,'on'); box(ax,'on');

% OTFS with channel estimation, one curve per SNR_p
h = gobjects(1, numel(p.SNRp_dB));
for pi = 1:numel(p.SNRp_dB)
    h(pi) = semilogy(x, BER_ce(pi,:), mk{mod(pi-1,3)+1}, ...
        'LineWidth',1.4,'MarkerSize',7,'Color',co(pi,:), ...
        'DisplayName',sprintf('SNR_p = %d dB',p.SNRp_dB(pi)));
end
% references
hO = semilogy(x, BER_ofdm,  'k--', 'LineWidth',1.6,'DisplayName','OFDM, Ideal');
hI = semilogy(x, BER_ideal, 'k:',  'LineWidth',1.6,'DisplayName','OTFS, Ideal');

set(gca,'YScale','log');
xlim([min(x) max(x)]); ylim([1e-5 1e-1]);
xlabel('SNR_d in dB'); ylabel('BER');
title(sprintf('Fig. 9: N=%d, M=%d, l_\\tau=%d, 4-QAM (frac. Doppler, full guard)', ...
      p.N, p.M, p.l_tau));
legend([h hO hI],'Location','southwest');
hold off;

outdir = fileparts(mfilename('fullpath'));
saveas(f, fullfile(outdir,'fig_9.png'));
end
