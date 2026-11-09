function plot_ber_noici(matfile)
% PLOT_BER_NOICI  BER vs SNR for OFDM and CP-OTFS from a mcrb_ce_ici run (uses the
% stored true BER, ber_data).  Saves a MATLAB .fig (and a .png for preview).
    thisdir=fileparts(mfilename('fullpath'));
    if nargin<1||isempty(matfile)
        matfile=fullfile(thisdir,'results','mcrb_ce_ici_noici_1000t_results.mat');
    end
    S=load(matfile);
    SNR=S.SNR_dbs(:); names=S.configs(:,1);
    cols=[0.85 0.33 0.10; 0 0.45 0.74]; mrk={'-o','-d'};

    regime='nominal Doppler + ICI';
    if isfield(S,'doppler_scale') && S.doppler_scale==0, regime='NO Doppler (nu=0)';
    elseif isfield(S,'use_ici') && ~S.use_ici, regime='Doppler ~444 Hz, NO ICI'; end

    hf=figure('Name','BER vs SNR (no ICI)','Position',[200 160 840 560],'Color','w');
    hold on; grid on; box on; lg={};
    for ci=1:size(names,1)
        y=S.ber_data(:,ci); y(y==0)=NaN;                 % mask zeros/NaN for log axis
        semilogy(SNR,y,mrk{ci},'Color',cols(ci,:),'LineWidth',2,'MarkerSize',6, ...
                 'MarkerFaceColor',cols(ci,:));
        lg=[lg,names(ci)]; %#ok<AGROW>
    end
    set(gca,'YScale','log'); xlabel('SNR [dB]'); ylabel('BER');
    title(sprintf('BER vs SNR, reduced EVA, %s (%d-QAM, %d trials)', regime, 4, S.nTrials));
    legend(lg,'Location','southwest','FontSize',10); set(gca,'FontSize',11);
    if exist('theme','file'), theme(hf,"light"); end

    [~,base]=fileparts(matfile); stem=regexprep(base,'_results$','');
    figdir=fullfile(thisdir,'figures'); if ~exist(figdir,'dir'), mkdir(figdir); end
    figpath=fullfile(figdir,[stem '_ber.fig']);          % .fig files -> eva/figures/
    pngpath=fullfile(fileparts(matfile),[stem '_ber.png']);
    savefig(hf, figpath);
    exportgraphics(hf, pngpath,'Resolution',150,'BackgroundColor','white');
    fprintf('Saved %s\nSaved %s\n', figpath, pngpath);

    fprintf('\nBER vs SNR:\n%4s | %10s %10s\n','SNR',names{1},names{2});
    for si=1:numel(SNR)
        if mod(SNR(si),3)~=0, continue; end
        fprintf('%4d | %10.2e %10.2e\n',SNR(si),S.ber_data(si,1),S.ber_data(si,2));
    end
end
