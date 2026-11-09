function plot_crb_paper_delay(matfile)
% PLOT_CRB_PAPER_DELAY  Pure delay-trace CRB (p=0, no detection errors -> F=0,
% LB = sigma2*G = CRB) vs SNR for the paper channel (16-PSK, N=128, M=1,
% delay-only targets): precoded OFDM / SC / OTFS + pilot-only reference.
% CRB is linear in sigma2 = alpha_ref2/10^(SNR/10), so curves are straight lines;
% the waveform gap is the constant Jensen offset (OFDM = Jensen exactly).
    thisdir=fileparts(mfilename('fullpath'));
    if nargin<1||isempty(matfile), matfile=fullfile(thisdir,'results','compute_mcrb_mc_true_paper_results.mat'); end
    S=load(matfile);
    [~,ip0]=min(abs(S.p_list-0));
    c=3e8; u=(c/2)^2;                                  % tau [s^2] -> range [m^2]
    nC=size(S.configs,1);
    tr0=zeros(1,nC);
    for ci=1:nC, tr0(ci)=median(S.trG_all(:,ip0,ci)); end   % CRB @ sigma2=1 [s^2]

    snr=-10:1:30; s2=S.alpha_ref2./10.^(snr/10);
    cols=[0.85 0.33 0.10; 0 0.45 0.74; 0.93 0.69 0.13];     % OFDM / SC / OTFS
    hf=figure('Position',[190 150 840 560],'Color','w'); hold on; grid on; box on;
    h=gobjects(1,nC); lg=cell(1,nC);
    for ci=1:nC
        h(ci)=semilogy(snr,u*tr0(ci)*s2,'-','Color',cols(ci,:),'LineWidth',2);
        lg{ci}=sprintf('%s  (%+.2f dB)',S.configs{ci,1},10*log10(tr0(ci)/tr0(1)));
    end
    if isfield(S,'trGpil') && isfinite(S.trGpil)
        hp=semilogy(snr,u*S.trGpil*s2,'k-.','LineWidth',1.8);
        h=[h hp]; lg=[lg,{sprintf('OFDM pilot-only (%.0f%%)',100*S.pilot_frac)}];
    end
    set(gca,'YScale','log');
    xlabel('SNR = |\alpha_1|^2/\sigma^2  [dB]');
    ylabel('tr\{CRB_R\} = (c/2)^2 \Sigma_p CRB(\tau_p)  [m^2]');
    title(sprintf(['Paper channel: delay-trace CRB vs SNR (p=0, no mismatch; %s, P=%d, %d trials)\n' ...
                   'OFDM sits exactly on the Jensen bound; spread waveforms +0.4-0.6 dB above'], ...
                   S.mod_type, S.P, S.nTrials));
    legend(h,lg,'Location','southwest','FontSize',10); set(gca,'FontSize',11);
    if exist('theme','file'), theme(hf,"light"); end
    figdir=fullfile(thisdir,'figures'); if ~exist(figdir,'dir'), mkdir(figdir); end
    savefig(hf,fullfile(figdir,'crb_paper_delay_vs_snr.fig'));
    exportgraphics(hf,fullfile(thisdir,'results','crb_paper_delay_vs_snr.png'),'Resolution',150,'BackgroundColor','white');
    fprintf('Saved %s and .png\n',fullfile(figdir,'crb_paper_delay_vs_snr.fig'));

    % ---- linear-axis version (pilot-only omitted: 26x higher, would squash the
    % waveform curves; the Jensen gap shows at the low-SNR end) ----------------
    hf2=figure('Position',[220 170 840 560],'Color','w'); hold on; grid on; box on;
    h2=gobjects(1,nC);
    for ci=1:nC
        h2(ci)=plot(snr,u*tr0(ci)*s2,'-','Color',cols(ci,:),'LineWidth',2);
    end
    xlabel('SNR = |\alpha_1|^2/\sigma^2  [dB]');
    ylabel('tr\{CRB_R\} = (c/2)^2 \Sigma_p CRB(\tau_p)  [m^2]');
    title(sprintf(['Paper channel: delay-trace CRB vs SNR, LINEAR axis (p=0; %s, P=%d, %d trials)\n' ...
                   'pilot-only line omitted (26x higher would flatten the waveform curves)'], ...
                   S.mod_type, S.P, S.nTrials));
    legend(h2,lg(1:nC),'Location','northeast','FontSize',10); set(gca,'FontSize',11);
    if exist('theme','file'), theme(hf2,"light"); end
    savefig(hf2,fullfile(figdir,'crb_paper_delay_vs_snr_linear.fig'));
    exportgraphics(hf2,fullfile(thisdir,'results','crb_paper_delay_vs_snr_linear.png'),'Resolution',150,'BackgroundColor','white');
    fprintf('Saved %s and .png\n',fullfile(figdir,'crb_paper_delay_vs_snr_linear.fig'));
end
