function plot_mcrb_confirm()
% PLOT_MCRB_CONFIRM  Before/after the OTFS guard shrink at B=15.4 MHz: delay & Doppler
% LB ratio (OTFS/OFDM) vs SNR for full guard (G=13) vs shrunk guard (G=9).  Ratio>1
% => OFDM wins.  Shows OFDM wins BOTH delay and Doppler at every low SNR with G=9.
    thisdir=fileparts(mfilename('fullpath'));
    S=load(fullfile(thisdir,'results','mcrb_sweep_confirm_results.mat'));
    res=S.res; nSNR=numel(S.SNR_dbs); SNR=S.SNR_dbs(:);
    % config 1 = g13 (full), config 2 = g9 (shrunk)
    dr=@(c) [res((c-1)*nSNR+(1:nSNR)).dr]; nr=@(c) [res((c-1)*nSNR+(1:nSNR)).nr];
    orange=[0.85 0.33 0.10]; blue=[0 0.45 0.74];

    hf=figure('Position',[160 130 980 460],'Color','w');

    subplot(1,2,1); hold on; grid on; box on;             % delay
    plot(SNR,dr(1),'--s','Color',orange,'LineWidth',1.8,'MarkerSize',6);
    plot(SNR,dr(2),'-o' ,'Color',orange,'LineWidth',2.2,'MarkerSize',7,'MarkerFaceColor',orange);
    yline(1,'k--','OFDM wins above','LineWidth',1.3,'FontSize',9,'LabelHorizontalAlignment','left');
    xlabel('SNR [dB]'); ylabel('delay LB ratio  LB_{OTFS}/LB_{OFDM}');
    legend({'G=13 (full guard)','G=9 (shrunk guard)'},'Location','northwest','FontSize',9);
    title('Delay'); ylim([0.9 1.2]); set(gca,'FontSize',11);

    subplot(1,2,2); hold on; grid on; box on;             % Doppler
    plot(SNR,nr(1),'--d','Color',blue,'LineWidth',1.8,'MarkerSize',6);
    plot(SNR,nr(2),'-d' ,'Color',blue,'LineWidth',2.2,'MarkerSize',7,'MarkerFaceColor',blue);
    yline(1,'k--','OFDM wins above','LineWidth',1.3,'FontSize',9,'LabelHorizontalAlignment','left');
    xlabel('SNR [dB]'); ylabel('Doppler LB ratio  LB_{OTFS}/LB_{OFDM}');
    legend({'G=13 (full guard)','G=9 (shrunk guard)'},'Location','southwest','FontSize',9);
    title('Doppler'); ylim([0 6]); set(gca,'FontSize',11);

    sgtitle(sprintf(['CP-OFDM vs CP-OTFS sensing LB, reduced EVA, B=15.4 MHz, %d trials: ' ...
                     'shrinking the OTFS guard -> OFDM wins BOTH'], S.nTrials));
    if exist('theme','file'), theme(hf,"light"); end
    figdir=fullfile(thisdir,'figures'); if ~exist(figdir,'dir'), mkdir(figdir); end
    savefig(hf,fullfile(figdir,'mcrb_confirm_guard.fig'));
    exportgraphics(hf,fullfile(thisdir,'results','mcrb_confirm_guard.png'),'Resolution',150,'BackgroundColor','white');
    fprintf('Saved %s and .png\n',fullfile(figdir,'mcrb_confirm_guard.fig'));
end
