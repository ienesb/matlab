function plot_mcrb_sweep_guard()
% PLOT_MCRB_SWEEP_GUARD  Delay & Doppler LB ratio (OTFS/OFDM) vs OTFS delay-guard
% size G, at the B=15.4 MHz sweet spot.  Ratio>1 => OFDM wins.  Shows the guard is
% OTFS's delay-LB advantage: shrinking it raises the delay ratio above 1 while
% Doppler stays >1 -> OFDM wins BOTH.  (G=8 is a pilot-lattice outlier -- 4 time
% pilots break OFDM's CE -- and is flagged/excluded from the trend.)
    thisdir=fileparts(mfilename('fullpath'));
    S=load(fullfile(thisdir,'results','mcrb_sweep_guard_results.mat'));
    res=S.res; nSNR=numel(S.SNR_dbs); nG=numel(res)/nSNR; SNRs=S.SNR_dbs;
    G=zeros(1,nG); DR=zeros(nG,nSNR); NR=zeros(nG,nSNR);
    for gi=1:nG
        idx=(gi-1)*nSNR+(1:nSNR);
        tok=regexp(res(idx(1)).label,'g(\d+)','tokens','once'); G(gi)=str2double(tok{1});
        DR(gi,:)=[res(idx).dr]; NR(gi,:)=[res(idx).nr];
    end
    [G,ord]=sort(G); DR=DR(ord,:); NR=NR(ord,:);
    bad=(G==8);                                              % pilot-lattice outlier
    orange=[0.85 0.33 0.10]; blue=[0 0.45 0.74]; sty={'-','--',':'};

    hf=figure('Position',[180 140 880 560],'Color','w'); hold on; grid on; box on; h=gobjects(0); lg={};
    for si=1:nSNR
        h(end+1)=plot(G(~bad),DR(~bad,si),[sty{si} 'o'],'Color',orange,'LineWidth',1.8,'MarkerSize',6,'MarkerFaceColor',orange); %#ok<AGROW>
        plot(G(~bad),NR(~bad,si),[sty{si} 'd'],'Color',blue,'LineWidth',1.8,'MarkerSize',6,'MarkerFaceColor',blue);
        lg=[lg,{sprintf('SNR=%d dB',SNRs(si))}]; %#ok<AGROW>
    end
    plot(G(bad),DR(bad,:),'x','Color',[.6 .6 .6],'MarkerSize',9,'LineWidth',1.5);   % G=8 outlier
    plot(G(bad),NR(bad,:),'x','Color',[.6 .6 .6],'MarkerSize',9,'LineWidth',1.5);
    yline(1,'k--','OFDM wins above','LineWidth',1.3,'FontSize',9,'LabelHorizontalAlignment','left');
    set(gca,'XDir','reverse');                               % shrinking guard -> rightward
    xlabel('OTFS delay-guard size  G  [bins]   (\leftarrow larger    smaller \rightarrow)');
    ylabel('LB ratio  LB_{OTFS}/LB_{OFDM}');
    text(9,3.5,'\bf Doppler ratio','Color',blue,'FontSize',10);
    text(9,1.25,'\bf delay ratio','Color',orange,'FontSize',10);
    legend(h,lg,'Location','east','FontSize',9);
    title(sprintf(['LB ratio vs OTFS guard size, N=512 M=128 \\Deltaf=30kHz (B=15.4 MHz), %d trials\n' ...
                   'shrinking the guard lifts delay ratio >1; Doppler stays >1 -> OFDM wins BOTH  (G=8 = pilot outlier, x)'], S.nTrials));
    set(gca,'FontSize',11); if exist('theme','file'), theme(hf,"light"); end
    figdir=fullfile(thisdir,'figures'); if ~exist(figdir,'dir'), mkdir(figdir); end
    savefig(hf,fullfile(figdir,'mcrb_sweep_guard.fig'));
    exportgraphics(hf,fullfile(thisdir,'results','mcrb_sweep_guard.png'),'Resolution',150,'BackgroundColor','white');
    fprintf('Saved %s and .png\n',fullfile(figdir,'mcrb_sweep_guard.fig'));

    fprintf('\n%3s | %-22s | %-22s\n','G','delay_r (SNR 0/2/4)','Doppler_r (SNR 0/2/4)');
    for gi=1:nG
        fl=tern(G(gi)==8,'  <- pilot outlier','');
        fprintf('%3d | %5.2f %5.2f %5.2f | %5.2f %5.2f %5.2f%s\n',G(gi),DR(gi,:),NR(gi,:),fl);
    end
end
function s=tern(c,a,b), if c, s=a; else, s=b; end, end
