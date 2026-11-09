% replot.m - regenerate Fig. 1 from saved results (crb_new_results.mat)
clear; clc; close all;
load('crb_new_results.mat');   % SNR_dB, crb_jensen, crb_sc/afdm/otfs/ofdm

make_fig(SNR_dB, crb_sc(:,1), crb_afdm(:,1), crb_otfs(:,1), crb_ofdm(:,1), ...
         crb_jensen, '-',  '(a) Ranging CRB under 16-QAM constellation', ...
         'fig_crb_16qam.png');
make_fig(SNR_dB, crb_sc(:,2), crb_afdm(:,2), crb_otfs(:,2), crb_ofdm(:,2), ...
         crb_jensen, 'kx', '(b) Ranging CRB under 16-PSK constellation', ...
         'fig_crb_16psk.png');

function make_fig(snr, sc, afdm, otfs, ofdm, jensen, jstyle, ttl, fname)
    fig = figure('Color','w','Position',[100 100 700 500]);
    ax  = axes('Parent',fig,'Color','w','XColor','k','YColor','k', ...
               'GridColor',[0.3 0.3 0.3]); hold(ax,'on');
    plot(ax, snr, sc,   '-',  'LineWidth',1.8, 'Color',[0.00 0.45 0.74]);
    plot(ax, snr, afdm, '-.', 'LineWidth',1.8, 'Color',[0.49 0.18 0.56]);
    plot(ax, snr, otfs, '--', 'LineWidth',1.8, 'Color',[0.93 0.69 0.13]);
    plot(ax, snr, ofdm, '-',  'LineWidth',1.8, 'Color',[0.85 0.33 0.10]);
    if strcmp(jstyle,'kx')
        plot(ax, snr, jensen, 'kx', 'MarkerSize',9, 'LineWidth',1.5);
    else
        plot(ax, snr, jensen, '-',  'LineWidth',1.8, 'Color','k');
    end
    grid(ax,'on'); box(ax,'on');
    xlabel(ax,'SNR (dB)','Color','k');
    ylabel(ax,'CRB (m^2)','Color','k');
    legend(ax,{'SC','AFDM','OTFS','OFDM','Jensen Bound'}, ...
           'Location','northeast','TextColor','k','Color','w','EdgeColor','k');
    title(ax,ttl,'Color','k');
    xlim(ax,[snr(1) snr(end)]);
    exportgraphics(fig,fname,'Resolution',150,'BackgroundColor','white');
end
