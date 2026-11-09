function plot_mcrb_nodopp_wf(matfile)
% PLOT_MCRB_NODOPP_WF  OFDM vs CP-OTFS delay AND Doppler bounds for a reduced-EVA
% mcrb_ce_ici run.  Convention: MCRB = N0*G + F ; BB = (eta_bar-eta_0)(.)^T ;
% LB = MCRB + BB.  Two figures: delay RMSE [ns] and Doppler RMSE [Hz], each with
% LB (solid), MCRB (dashed), CRB (dotted) for both waveforms.  Regime (nominal /
% no-ICI / no-Doppler) is read from the file.  CP-OTFS may stop early (error-free
% beyond its SNR cap -> LB = CRB there).
    thisdir=fileparts(mfilename('fullpath'));
    if nargin<1||isempty(matfile)
        matfile=fullfile(thisdir,'results','mcrb_ce_ici_nodopp_2000_results.mat');
    end
    S=load(matfile);
    wf={'OFDM ICI','CP-OTFS ICI'}; lab={'OFDM','CP-OTFS'};
    cols=[0.85 0.33 0.10; 0 0.45 0.74];

    regime = 'nominal Doppler + ICI';
    if isfield(S,'doppler_scale') && S.doppler_scale==0
        regime = 'NO Doppler (nu=0)';
    elseif isfield(S,'use_ici') && ~S.use_ici
        regime = 'Doppler ~444 Hz, NO ICI';
    end
    [~,base] = fileparts(matfile);  stem = regexprep(base,'_results$','');

    % param, index block, unit scale, ylabel, filename tag
    pars = { 'delay',   S.idx_tau, 1e9, 'delay RMSE bound [ns]',    'delay'
             'Doppler', S.idx_nu,  1,   'Doppler RMSE bound [Hz]',  'doppler' };

    for pp = 1:size(pars,1)
        idx = pars{pp,2};  sc = pars{pp,3};
        hf=figure('Position',[180 140 860 570],'Color','w'); hold on; grid on; box on;
        h=gobjects(0); lg={};
        for w=1:2
            [SNR,mc,lb,cr]=curves(S,wf{w},idx);
            h(end+1)=semilogy(SNR,lb*sc,'-o' ,'Color',cols(w,:),'LineWidth',2.1,'MarkerSize',5); %#ok<AGROW>
            h(end+1)=semilogy(SNR,mc*sc,'--s','Color',cols(w,:),'LineWidth',1.7,'MarkerSize',5); %#ok<AGROW>
            h(end+1)=semilogy(SNR,cr*sc,':' ,'Color',cols(w,:),'LineWidth',1.5); %#ok<AGROW>
            lg=[lg,{[lab{w} ': LB'],[lab{w} ': MCRB'],[lab{w} ': CRB']}]; %#ok<AGROW>
        end
        set(gca,'YScale','log'); xlabel('SNR [dB]'); ylabel(pars{pp,4});
        legend(h,lg,'Location','northeast','FontSize',9);
        title(sprintf(['OFDM vs CP-OTFS: %s bounds, reduced EVA, %s, %d trials\n' ...
                       'LB = MCRB + BB ;  MCRB = N0 G + F ;  CRB = genie'], pars{pp,1}, regime, S.nTrials));
        set(gca,'FontSize',11); if exist('theme','file'), theme(hf,"light"); end
        png=fullfile(fileparts(matfile),sprintf('%s_ofdm_vs_otfs_%s.png',stem,pars{pp,5}));
        exportgraphics(hf,png,'Resolution',150,'BackgroundColor','white'); fprintf('Saved %s\n',png);
        figdir=fullfile(thisdir,'figures'); if ~exist(figdir,'dir'), mkdir(figdir); end
        figp=fullfile(figdir,sprintf('%s_ofdm_vs_otfs_%s.fig',stem,pars{pp,5}));  % .fig -> eva/figures/
        savefig(hf,figp); fprintf('Saved %s\n',figp);
    end

    % console: Doppler LB/MCRB/CRB
    fprintf('\nDoppler RMSE [Hz]:  LB / MCRB / CRB\n');
    fprintf('%4s | %-24s | %-24s\n','SNR','OFDM (LB/MCRB/CRB)','CP-OTFS (LB/MCRB/CRB)');
    [So,mo,lo,co]=curves(S,'OFDM ICI',S.idx_nu); [~,mt,lt,ct]=curves(S,'CP-OTFS ICI',S.idx_nu);
    for x=[0 3 6 9 12 15]
        i=find(So==x); if isempty(i), continue; end
        fprintf('  %2d | %7.2f %7.2f %7.2f | %7.2f %7.2f %7.2f\n', x, ...
            lo(i),mo(i),co(i), lt(i),mt(i),ct(i));
    end
end

function [SNR,mcrb,lb,crb]=curves(S,wfname,idx)
    ci=find(strcmp(S.configs(:,1),wfname),1); SNR=S.SNR_dbs(:);
    n=numel(SNR); mcrb=nan(n,1); lb=nan(n,1); crb=nan(n,1);
    for si=1:n
        N0=1/db2pow(SNR(si));
        MC=N0*S.G_bar(:,:,si,ci)+S.F_bar(:,:,si,ci);
        LBm=MC+S.BB_bar(:,:,si,ci); CR=N0*S.Gc_bar(:,:,si,ci);
        mcrb(si)=sqrt(sum(real(diag(MC(idx,idx)))));
        lb(si)  =sqrt(sum(real(diag(LBm(idx,idx)))));
        crb(si) =sqrt(sum(real(diag(CR(idx,idx)))));
    end
end
