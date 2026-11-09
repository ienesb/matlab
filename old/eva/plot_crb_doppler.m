function plot_crb_doppler(matfile)
% PLOT_CRB_DOPPLER  Trace Doppler CRB vs SNR for the Doppler-only ('paperdopp')
% run of run_fim_avg.  The FIM is linear in 1/sigma2, so CRB(SNR) = CRB(sigma2=1)
% * sigma2 with sigma2 = |alpha_1|^2 / 10^(SNR/10); here |alpha_1|^2 = 1.
% CP-OFDM vs CP-OTFS, ICI and no-ICI.  Saves crb_doppler_vs_snr.png.
    thisdir=fileparts(mfilename('fullpath'));
    if nargin<1||isempty(matfile), matfile=fullfile(thisdir,'results','eva_fim_avg_results.mat'); end
    S=load(matfile);
    assert(isfield(S,'estimate_doppler') && S.estimate_doppler, 'run has no Doppler CRB');

    % trace Doppler CRB at sigma2 = 1 (sum over paths of the trial-mean per-path CRB)
    trI_OF  = sum(mean(S.CRBn_OF_I,1));    trI_OT  = sum(mean(S.CRBn_CPOT_I,1));
    trN_OF  = sum(mean(S.CRBn_OF_N,1));    trN_OT  = sum(mean(S.CRBn_CPOT_N,1));

    alpha_ref2 = abs(S.alpha_vec(1))^2;
    snr = -10:1:30;  s2 = alpha_ref2 ./ 10.^(snr/10);      % sigma2(SNR)

    hf=figure('Position',[180 140 860 570],'Color','w'); hold on; grid on; box on;
    cOF=[0.85 0.33 0.10]; cOT=[0 0.45 0.74];
    semilogy(snr, trI_OF*s2, '-o','Color',cOF,'LineWidth',2.0,'MarkerSize',4);
    semilogy(snr, trI_OT*s2, '-o','Color',cOT,'LineWidth',2.0,'MarkerSize',4);
    semilogy(snr, trN_OF*s2, '--s','Color',cOF,'LineWidth',1.7,'MarkerSize',4);
    semilogy(snr, trN_OT*s2, '--s','Color',cOT,'LineWidth',1.7,'MarkerSize',4);
    set(gca,'YScale','log'); xlabel('SNR = |\alpha_1|^2/\sigma^2  [dB]');
    ylabel('tr\{CRB_\nu\} = \Sigma_p CRB(\nu_p)  [Hz^2]');
    legend({sprintf('CP-OFDM (ICI)   %.0f dB gap',10*log10(trI_OT/trI_OF)), ...
            'CP-OTFS (ICI)', ...
            sprintf('CP-OFDM (no-ICI)  %.0f dB gap',10*log10(trN_OT/trN_OF)), ...
            'CP-OTFS (no-ICI)'}, 'Location','northeast','FontSize',9);
    title(sprintf(['Doppler-only channel (%s): trace Doppler CRB vs SNR\n' ...
                   'N=1, M=%d, %d targets, %s; OFDM sits ~%.1f dB below OTFS (Jensen)'], ...
                   S.channel_type, S.M_o, S.P, S.modulation, 10*log10(trN_OT/trN_OF)));
    set(gca,'FontSize',11); if exist('theme','file'), theme(hf,"light"); end
    png=fullfile(fileparts(matfile),'crb_doppler_vs_snr.png');
    exportgraphics(hf,png,'Resolution',150,'BackgroundColor','white');
    fprintf('Saved %s\n',png);
    fprintf('tr CRB(nu) @ sigma2=1 [Hz^2]:  OFDM(noICI)=%.3e  OTFS(noICI)=%.3e  (OTFS +%.1f dB)\n', ...
            trN_OF, trN_OT, 10*log10(trN_OT/trN_OF));
    fprintf('                               OFDM(ICI)=%.3e  OTFS(ICI)=%.3e  (OTFS +%.1f dB)\n', ...
            trI_OF, trI_OT, 10*log10(trI_OT/trI_OF));
end
