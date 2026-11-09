function plot_mcrb_nodopp_compare(iciFile, noiciFile)
% PLOT_MCRB_NODOPP_COMPARE  Delay MCRB and LB vs SNR, WITH ICI vs WITHOUT ICI
% (no Doppler), for OFDM-ICI and CP-OTFS-ICI on the reduced-EVA channel.
% Convention: MCRB = A^-1 B A^-1 = N0*G + F ; BB = (eta_bar-eta_0)(.)^T ;
%             LB = MCRB + BB.  Delay = sum over the tau block, reported as RMSE [ns].
% Recomputed from the raw G_bar/F_bar/BB_bar so both files use one definition.
    thisdir=fileparts(mfilename('fullpath'));
    if nargin<1||isempty(iciFile),   iciFile  =fullfile(thisdir,'results','mcrb_ce_ici_results_2000.mat'); end
    if nargin<2||isempty(noiciFile), noiciFile=fullfile(thisdir,'results','mcrb_ce_ici_nodopp_2000_results.mat'); end
    I=load(iciFile); Z=load(noiciFile);

    wf={'OFDM ICI','CP-OTFS ICI'}; ttl={'OFDM','CP-OTFS'};
    for w=1:numel(wf)
        [Is,Im,Il,Ic]=curves(I,wf{w});     % ICI:  SNR, MCRB, LB, CRB  (delay RMSE [s])
        [Zs,Zm,Zl,Zc]=curves(Z,wf{w});     % noICI
        hf=figure('Position',[180 140 840 560],'Color','w'); hold on; grid on; box on;
        cIC=[0.85 0.33 0.10]; cNO=[0 0.45 0.74]; s=1e9;   % s -> ns
        p1=semilogy(Is,Il*s,'-o' ,'Color',cIC,'LineWidth',2.1,'MarkerSize',5);
        p2=semilogy(Is,Im*s,'--s','Color',cIC,'LineWidth',1.7,'MarkerSize',5);
        p3=semilogy(Zs,Zl*s,'-o' ,'Color',cNO,'LineWidth',2.1,'MarkerSize',5);
        p4=semilogy(Zs,Zm*s,'--s','Color',cNO,'LineWidth',1.7,'MarkerSize',5);
        p5=semilogy(Is,Ic*s,':' ,'Color',[.35 .35 .35],'LineWidth',1.6);
        set(gca,'YScale','log'); xlabel('SNR [dB]'); ylabel('delay RMSE bound [ns]');
        legend([p1 p2 p3 p4 p5], ...
            {'LB, with ICI','MCRB, with ICI','LB, no ICI','MCRB, no ICI','CRB (genie)'}, ...
            'Location','northeast','FontSize',9);
        title(sprintf(['%s: delay MCRB & LB vs SNR, reduced EVA, %d trials\n' ...
                       'ICI (nominal Doppler) vs no ICI (nu=0);  LB = MCRB + BB'], ttl{w}, I.nTrials));
        set(gca,'FontSize',11); if exist('theme','file'), theme(hf,"light"); end
        png=fullfile(fileparts(noiciFile),sprintf('mcrb_nodopp_compare_%s.png',lower(ttl{w})));
        exportgraphics(hf,png,'Resolution',150,'BackgroundColor','white');
        fprintf('Saved %s\n',png);
    end

    % ---- console: OFDM delay LB & MCRB, ICI vs no-ICI, and the BB share -------
    [Is,Im,Il,~]=curves(I,'OFDM ICI'); [~,Zm,Zl,~]=curves(Z,'OFDM ICI');
    fprintf('\nOFDM delay RMSE [ns]:  LB(ICI) MCRB(ICI) | LB(noICI) MCRB(noICI)\n');
    for x=[0 3 6 9 12 15 20 25]
        i=find(Is==x); if isempty(i), continue; end
        fprintf('  %2d dB | %7.3f %7.3f | %7.3f %7.3f\n', x, Il(i)*1e9,Im(i)*1e9, Zl(i)*1e9,Zm(i)*1e9);
    end
end

function [SNR,mcrb,lb,crb]=curves(S,wfname)
    ci=find(strcmp(S.configs(:,1),wfname),1); it=S.idx_tau; SNR=S.SNR_dbs(:);
    n=numel(SNR); mcrb=nan(n,1); lb=nan(n,1); crb=nan(n,1);
    for si=1:n
        N0=1/db2pow(SNR(si));
        MC=N0*S.G_bar(:,:,si,ci)+S.F_bar(:,:,si,ci);
        LBm=MC+S.BB_bar(:,:,si,ci); CR=N0*S.Gc_bar(:,:,si,ci);
        mcrb(si)=sqrt(sum(real(diag(MC(it,it)))));
        lb(si)  =sqrt(sum(real(diag(LBm(it,it)))));
        crb(si) =sqrt(sum(real(diag(CR(it,it)))));
    end
end
