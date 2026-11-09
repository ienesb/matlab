function plot_lb_g13_absolute(matfile)
% PLOT_LB_G13_ABSOLUTE  Absolute delay & Doppler LB (RMSE bound) vs SNR for the
% G=13 config (reduced EVA, N=512 M=128 deltaf=30kHz, B=15.4 MHz, Doppler+ICI),
% CP-OFDM vs CP-OTFS.  LB = MCRB + BB (the realistic-detector bound).  Two panels:
% delay [ns] and Doppler [Hz].
    thisdir=fileparts(mfilename('fullpath'));
    if nargin<1||isempty(matfile), matfile=fullfile(thisdir,'results','mcrb_confirm_g13_full_results.mat'); end
    S=load(matfile); res=S.res; SNR=S.SNR_dbs(:); nS=numel(SNR);
    dOF=[res(1:nS).dOF].'; dOT=[res(1:nS).dOT].';        % delay LB [s]
    nOF=[res(1:nS).nOF].'; nOT=[res(1:nS).nOT].';        % Doppler LB [Hz]
    orange=[0.85 0.33 0.10]; blue=[0 0.45 0.74];

    hf=figure('Position',[150 130 1000 460],'Color','w');

    subplot(1,2,1); hold on; grid on; box on;             % delay
    semilogy(SNR,dOF*1e9,'-o','Color',orange,'LineWidth',2,'MarkerSize',5,'MarkerFaceColor',orange);
    semilogy(SNR,dOT*1e9,'-d','Color',blue,'LineWidth',2,'MarkerSize',5,'MarkerFaceColor',blue);
    set(gca,'YScale','log'); xlabel('SNR [dB]'); ylabel('delay LB (RMSE) [ns]');
    legend({'CP-OFDM','CP-OTFS'},'Location','northeast','FontSize',9);
    title('(a) Delay LB vs SNR'); set(gca,'FontSize',11);

    subplot(1,2,2); hold on; grid on; box on;             % Doppler
    semilogy(SNR,nOF,'-o','Color',orange,'LineWidth',2,'MarkerSize',5,'MarkerFaceColor',orange);
    semilogy(SNR,nOT,'-d','Color',blue,'LineWidth',2,'MarkerSize',5,'MarkerFaceColor',blue);
    set(gca,'YScale','log'); xlabel('SNR [dB]'); ylabel('Doppler LB (RMSE) [Hz]');
    legend({'CP-OFDM','CP-OTFS'},'Location','northeast','FontSize',9);
    title('(b) Doppler LB vs SNR'); set(gca,'FontSize',11);

    sgtitle(sprintf(['Sensing LB (=MCRB+BB) vs SNR, reduced EVA (B=15.4 MHz, Doppler+ICI, G=13, %d trials)\n' ...
                     'OFDM lower only for Doppler at low SNR; OTFS lower elsewhere (its cliff -> smaller BB)'], S.nTrials));
    if exist('theme','file'), theme(hf,"light"); end
    figdir=fullfile(thisdir,'figures'); if ~exist(figdir,'dir'), mkdir(figdir); end
    savefig(hf,fullfile(figdir,'lb_g13_absolute.fig'));
    exportgraphics(hf,fullfile(thisdir,'results','lb_g13_absolute.png'),'Resolution',150,'BackgroundColor','white');
    fprintf('Saved %s and .png\n',fullfile(figdir,'lb_g13_absolute.fig'));

    fprintf('\n%4s | %-22s | %-22s\n','SNR','delay LB [ns] OFDM/OTFS','Doppler LB [Hz] OFDM/OTFS');
    for si=1:nS
        if mod(SNR(si),3)==0
            fprintf('  %2d | %9.3f %9.3f | %9.3f %9.3f\n',SNR(si),dOF(si)*1e9,dOT(si)*1e9,nOF(si),nOT(si));
        end
    end
end
