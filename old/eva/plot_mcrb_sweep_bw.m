function plot_mcrb_sweep_bw(matfile)
% PLOT_MCRB_SWEEP_BW  Delay & Doppler LB ratio (OTFS/OFDM) vs bandwidth B=N*deltaf
% from mcrb_sweep_grid.  Ratio>1 => OFDM wins.  Shows the bandwidth trend: delay
% ratio climbs toward/past 1 as B grows while Doppler stays >1 (OFDM wins both).
    thisdir=fileparts(mfilename('fullpath'));
    if nargin<1||isempty(matfile), matfile=fullfile(thisdir,'results','mcrb_sweep_grid_results.mat'); end
    S=load(matfile); res=S.res; nSNR=numel(S.SNR_dbs); nG=numel(res)/nSNR;
    B=zeros(1,nG); drm=B; drl=B; drh=B; nrm=B; nrl=B; nrh=B; lab=cell(1,nG);
    for gi=1:nG
        idx=(gi-1)*nSNR+(1:nSNR);
        B(gi)=res(idx(1)).N*res(idx(1)).deltaf/1e6;                 % MHz
        dr=[res(idx).dr]; nr=[res(idx).nr];
        drm(gi)=median(dr); drl(gi)=min(dr); drh(gi)=max(dr);
        nrm(gi)=median(nr); nrl(gi)=min(nr); nrh(gi)=max(nr);
        lab{gi}=res(idx(1)).label;
    end
    [B,ord]=sort(B); drm=drm(ord);drl=drl(ord);drh=drh(ord); nrm=nrm(ord);nrl=nrl(ord);nrh=nrh(ord); lab=lab(ord);
    orange=[0.85 0.33 0.10]; blue=[0 0.45 0.74];

    hf=figure('Position',[180 140 880 560],'Color','w'); hold on; grid on; box on;
    jit=1+0.012*(1:nG);                                             % tiny x-jitter so overlapping-B points show
    errorbar(B.*jit,drm,drm-drl,drh-drm,'o','Color',orange,'MarkerFaceColor',orange,'LineWidth',1.6,'MarkerSize',7,'CapSize',4);
    errorbar(B./jit,nrm,nrm-nrl,nrh-nrm,'d','Color',blue,'MarkerFaceColor',blue,'LineWidth',1.6,'MarkerSize',7,'CapSize',4);
    yline(1,'k--','OFDM wins above','LineWidth',1.3,'FontSize',9,'LabelHorizontalAlignment','left');
    set(gca,'YScale','log'); xlabel('bandwidth  B = N\cdot\Deltaf  [MHz]'); ylabel('LB ratio  (LB_{OTFS}/LB_{OFDM})');
    legend({'delay ratio','Doppler ratio'},'Location','east','FontSize',10);
    title(sprintf(['LB ratio vs bandwidth, reduced EVA, low SNR (%s dB, %d trials)\n' ...
                   'delay ratio crosses 1 as B grows; Doppler stays >1 -> OFDM wins both at high B'], ...
                   num2str(S.SNR_dbs), S.nTrials));
    set(gca,'FontSize',11); if exist('theme','file'), theme(hf,"light"); end
    figdir=fullfile(thisdir,'figures'); if ~exist(figdir,'dir'), mkdir(figdir); end
    savefig(hf,fullfile(figdir,'mcrb_sweep_bw.fig'));
    exportgraphics(hf,fullfile(thisdir,'results','mcrb_sweep_bw.png'),'Resolution',150,'BackgroundColor','white');
    fprintf('Saved %s and .png\n',fullfile(figdir,'mcrb_sweep_bw.fig'));

    fprintf('\n%-28s %8s | %14s | %14s\n','config','B[MHz]','delay_r med[min,max]','dopp_r med[min,max]');
    for gi=1:nG
        fprintf('%-28s %8.2f | %.2f [%.2f,%.2f] | %.2f [%.2f,%.2f]\n', ...
            lab{gi},B(gi),drm(gi),drl(gi),drh(gi),nrm(gi),nrl(gi),nrh(gi));
    end
end
