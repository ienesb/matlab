function plot_ber_rayleigh_vs_constmag()
% PLOT_BER_RAYLEIGH_VS_CONSTMAG  Overlay the reduced-EVA 3-tap BER for the two
% tap-fading models: RAYLEIGH (alpha_p~CN(0,power_p)) vs CONST-MAGNITUDE
% (|alpha_p| fixed, random phase -- the usual EVA model).  No Doppler, no ICI.
% OFDM = orange, CP-OTFS = blue; Rayleigh = solid, const-mag = dashed.
    thisdir=fileparts(mfilename('fullpath'));
    R=load(fullfile(thisdir,'results','ber_rayleigh_eva_results.mat'));   % Rayleigh
    C=load(fullfile(thisdir,'results','ber_constmag_eva_results.mat'));   % const-mag
    orange=[0.85 0.33 0.10]; blue=[0 0.45 0.74];

    hf=figure('Name','BER: Rayleigh vs const-mag','Position',[200 150 880 580],'Color','w');
    hold on; grid on; box on; h=gobjects(1,5);
    h(5)=semilogy(R.SNR_dbs,R.ber_theory,':','Color',[.4 .4 .4],'LineWidth',1.5);
    h(1)=pb(R.SNR_dbs,R.ber(:,1),'-o', orange);   % OFDM  Rayleigh
    h(2)=pb(R.SNR_dbs,R.ber(:,2),'-d', blue);     % OTFS  Rayleigh
    h(3)=pb(C.SNR_dbs,C.ber(:,1),'--o',orange);   % OFDM  const-mag
    h(4)=pb(C.SNR_dbs,C.ber(:,2),'--d',blue);     % OTFS  const-mag
    set(gca,'YScale','log'); xlabel('SNR = E_s/N_0 [dB]'); ylabel('BER');
    title(sprintf(['Reduced-EVA 3-tap: Rayleigh vs const-magnitude taps ' ...
                   '(no Doppler, no ICI, %d trials)'], R.nTrials));
    legend(h,{'OFDM, Rayleigh','CP-OTFS, Rayleigh','OFDM, const-|amp|', ...
              'CP-OTFS, const-|amp|','order-1 Rayleigh ref'}, ...
           'Location','southwest','FontSize',9);
    set(gca,'FontSize',11); if exist('theme','file'), theme(hf,"light"); end

    figdir=fullfile(thisdir,'figures'); if ~exist(figdir,'dir'), mkdir(figdir); end
    savefig(hf,fullfile(figdir,'ber_rayleigh_vs_constmag.fig'));
    exportgraphics(hf,fullfile(thisdir,'results','ber_rayleigh_vs_constmag.png'),'Resolution',150,'BackgroundColor','white');
    fprintf('Saved %s and .png\n',fullfile(figdir,'ber_rayleigh_vs_constmag.fig'));

    fprintf('\n%4s | %-21s | %-21s\n','SNR','OFDM (Rayl / const)','CP-OTFS (Rayl / const)');
    for si=1:numel(R.SNR_dbs)
        if mod(R.SNR_dbs(si),4)~=0, continue; end
        fprintf('%4d | %9.2e %9.2e | %9.2e %9.2e\n', R.SNR_dbs(si), ...
            R.ber(si,1),C.ber(si,1), R.ber(si,2),C.ber(si,2));
    end
end

function h=pb(SNR,ber,sty,col)
    y=ber(:); y(y==0)=NaN;
    h=semilogy(SNR(:),y,sty,'Color',col,'LineWidth',2,'MarkerSize',6,'MarkerFaceColor',col);
end
