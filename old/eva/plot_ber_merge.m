function plot_ber_merge()
% PLOT_BER_MERGE  BER vs SNR merging the ICI (full: Doppler+ICI) and the
% no-Doppler cases on one axis.  OFDM = orange, CP-OTFS = blue; case shown by
% line style (solid = ICI+Doppler, dashed = no Doppler).  Saves .fig to figures/.
    thisdir=fileparts(mfilename('fullpath'));
    A=load(fullfile(thisdir,'results','mcrb_ce_ici_results_2000.mat'));      % ICI + Doppler
    B=load(fullfile(thisdir,'results','mcrb_ce_ici_nodopp_2000_results.mat')); % no Doppler
    orange=[0.85 0.33 0.10];  blue=[0 0.45 0.74];

    hf=figure('Name','BER: ICI vs no Doppler','Position',[200 160 860 570],'Color','w');
    hold on; grid on; box on;
    h=gobjects(1,4);
    h(1)=plot_ber(A.SNR_dbs, A.ber_data(:,1), '-o',  orange);   % OFDM  ICI+Dopp
    h(2)=plot_ber(A.SNR_dbs, A.ber_data(:,2), '-o',  blue);     % OTFS  ICI+Dopp
    h(3)=plot_ber(B.SNR_dbs, B.ber_data(:,1), '--s', orange);   % OFDM  no Doppler
    h(4)=plot_ber(B.SNR_dbs, B.ber_data(:,2), '--s', blue);     % OTFS  no Doppler
    set(gca,'YScale','log'); xlabel('SNR [dB]'); ylabel('BER');
    title(sprintf('BER vs SNR, reduced EVA (4-QAM): ICI+Doppler vs no Doppler'));
    legend(h, {'OFDM (ICI+Doppler)','CP-OTFS (ICI+Doppler)', ...
               'OFDM (no Doppler)','CP-OTFS (no Doppler)'}, ...
           'Location','southwest','FontSize',9);
    set(gca,'FontSize',11); if exist('theme','file'), theme(hf,"light"); end

    figdir=fullfile(thisdir,'figures'); if ~exist(figdir,'dir'), mkdir(figdir); end
    figpath=fullfile(figdir,'mcrb_ce_ici_ber_ici_vs_nodopp.fig');
    pngpath=fullfile(thisdir,'results','mcrb_ce_ici_ber_ici_vs_nodopp.png');
    savefig(hf, figpath);
    exportgraphics(hf, pngpath,'Resolution',150,'BackgroundColor','white');
    fprintf('Saved %s\nSaved %s\n', figpath, pngpath);
end

function h=plot_ber(SNR, ber, sty, col)
    y=ber(:); y(y==0)=NaN;                     % mask zeros/NaN for the log axis
    h=semilogy(SNR(:), y, sty, 'Color',col,'LineWidth',2,'MarkerSize',6,'MarkerFaceColor',col);
end
