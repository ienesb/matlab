function ber_confirm_salvage()
% BER_CONFIRM_SALVAGE  Partial BER from the stopped ber_confirm_config run (the
% script saved only at the end, so these are the points printed to the log).
% Config: reduced EVA, N=512 M=128 deltaf=30kHz (B=15.4 MHz), Doppler+ICI, OTFS
% shrunk guard G=9, 200 trials.  KEY finding: with the shrunk guard, CP-OTFS does
% NOT cliff -- it FLOORS at ~2e-2 (worse than OFDM past ~12 dB), because pulling
% data close to the pilot degrades the OTFS threshold CE.  So the guard shrink
% that gave OFDM the lower sensing LB also erases OTFS's BER advantage.
    thisdir=fileparts(mfilename('fullpath'));
    SNR_of=[0 5 10 15 20 25];  ber_of=[2.109e-1 9.997e-2 3.451e-2 1.027e-2 3.357e-3 1.309e-3];
    SNR_ot=[0 5 10 15];        ber_ot=[2.099e-1 9.257e-2 2.060e-2 1.927e-2];
    orange=[0.85 0.33 0.10]; blue=[0 0.45 0.74];

    hf=figure('Position',[200 160 840 560],'Color','w'); hold on; grid on; box on;
    semilogy(SNR_of,ber_of,'-o','Color',orange,'LineWidth',2,'MarkerSize',7,'MarkerFaceColor',orange);
    semilogy(SNR_ot,ber_ot,'-d','Color',blue,'LineWidth',2,'MarkerSize',7,'MarkerFaceColor',blue);
    set(gca,'YScale','log'); xlabel('SNR [dB]'); ylabel('BER'); xlim([0 25]);
    title(sprintf(['BER (partial), shrunk-guard config: reduced EVA, B=15.4 MHz, Doppler+ICI, OTFS G=9\n' ...
                   'CP-OTFS FLOORS at ~2e-2 (shrunk guard degrades its CE) -- no BER advantage here']));
    legend({'OFDM (dftdd CE)','CP-OTFS (MP), run stopped @15 dB'},'Location','southwest','FontSize',10);
    set(gca,'FontSize',11); if exist('theme','file'), theme(hf,"light"); end
    figdir=fullfile(thisdir,'figures'); if ~exist(figdir,'dir'), mkdir(figdir); end
    save(fullfile(thisdir,'results','ber_confirm_config_partial.mat'),'SNR_of','ber_of','SNR_ot','ber_ot');
    savefig(hf,fullfile(figdir,'ber_confirm_config_partial.fig'));
    exportgraphics(hf,fullfile(thisdir,'results','ber_confirm_config_partial.png'),'Resolution',150,'BackgroundColor','white');
    fprintf('Saved partial results (.mat/.fig/.png)\n');
end
