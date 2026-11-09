function plot_ber_lb_combined(lbfile, berfile)
% PLOT_BER_LB_COMBINED  The whole "BER != sensing" thesis in one figure, for the
% G=13 config (reduced EVA, N=512 M=128 deltaf=30kHz, B=15.4 MHz, Doppler+ICI).
% Left : BER vs SNR -- CP-OTFS cliffs (better COMMS).
% Right: sensing LB ratio LB_OTFS/LB_OFDM vs SNR -- >1 means OFDM has the lower
%        bound; OFDM wins Doppler and ties/edges delay at low SNR (better SENSING),
%        both converging to 1 as detection cleans up.  Same channel & receivers.
    thisdir=fileparts(mfilename('fullpath'));
    if nargin<1||isempty(lbfile),  lbfile =fullfile(thisdir,'results','mcrb_confirm_g13_ext_results.mat'); end
    if nargin<2||isempty(berfile), berfile=fullfile(thisdir,'results','ber_confirm_config_g13_results.mat'); end
    B=load(berfile); L=load(lbfile);
    res=L.res; nS=numel(L.SNR_dbs); SNRl=L.SNR_dbs(:);
    dr=[res(1:nS).dr].'; nr=[res(1:nS).nr].';           % single config
    orange=[0.85 0.33 0.10]; blue=[0 0.45 0.74];

    hf=figure('Position',[150 130 1000 460],'Color','w');

    % --- (a) BER ---
    subplot(1,2,1); hold on; grid on; box on;
    y=B.ber(:,1); y(y==0)=NaN; semilogy(B.SNR_dbs,y,'-o','Color',orange,'LineWidth',2,'MarkerSize',5,'MarkerFaceColor',orange);
    y=B.ber(:,2); y(y==0)=NaN; semilogy(B.SNR_dbs,y,'-d','Color',blue,'LineWidth',2,'MarkerSize',5,'MarkerFaceColor',blue);
    set(gca,'YScale','log'); xlabel('SNR [dB]'); ylabel('BER'); xlim([0 25]);
    legend({'CP-OFDM','CP-OTFS'},'Location','southwest','FontSize',9);
    title('(a) Communications: BER   (CP-OTFS wins)'); set(gca,'FontSize',11);

    % --- (b) sensing LB ratio ---
    subplot(1,2,2); hold on; grid on; box on;
    plot(SNRl,dr,'-o','Color',orange,'LineWidth',2,'MarkerSize',6,'MarkerFaceColor',orange);
    plot(SNRl,nr,'-d','Color',blue,'LineWidth',2,'MarkerSize',6,'MarkerFaceColor',blue);
    yline(1,'k--','OFDM has lower LB above','LineWidth',1.3,'FontSize',9,'LabelHorizontalAlignment','left');
    set(gca,'YScale','log'); xlabel('SNR [dB]'); ylabel('LB ratio  LB_{OTFS}/LB_{OFDM}');
    legend({'delay','Doppler'},'Location','northeast','FontSize',9);
    title('(b) Sensing: LB ratio   (CP-OFDM wins at low SNR)'); set(gca,'FontSize',11);

    sgtitle(sprintf(['CP-OFDM vs CP-OTFS, reduced EVA (B=15.4 MHz, Doppler+ICI): ' ...
                     'OTFS better COMMS (BER), OFDM better low-SNR SENSING (LB)']));
    if exist('theme','file'), theme(hf,"light"); end
    figdir=fullfile(thisdir,'figures'); if ~exist(figdir,'dir'), mkdir(figdir); end
    savefig(hf,fullfile(figdir,'ber_lb_combined_g13.fig'));
    exportgraphics(hf,fullfile(thisdir,'results','ber_lb_combined_g13.png'),'Resolution',150,'BackgroundColor','white');
    fprintf('Saved %s and .png\n',fullfile(figdir,'ber_lb_combined_g13.fig'));
end
