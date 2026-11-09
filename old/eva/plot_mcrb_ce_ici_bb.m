function [hf_tau, hf_nu] = plot_mcrb_ce_ici_bb(matfile)
% PLOT_MCRB_CE_ICI_BB  MCRB decomposition on ONE figure per parameter.
%
% For each waveform, three curves are overlaid so the pseudo-true bias term is
% visible directly as the gap between the first two:
%   solid  : FULL MCRB              sqrt(tr{N0*G + F + BB})
%   dashed : MCRB WITHOUT BB        sqrt(tr{N0*G + F})  = sqrt(tr{A^-1 B A^-1})
%   dotted : the BB term alone      sqrt(tr{BB})        = the pseudo-true bias
%
% Reading it: where solid and dashed coincide, BB is negligible and the bound is
% noise-limited (N0*G); where they separate, the detected-symbol errors are
% biasing the pseudo-true parameter.  Note F is ~15 orders below N0*G in this
% setup, so "without BB" is effectively the noise term.
%
%   plot_mcrb_ce_ici_bb                 % default results file
%   plot_mcrb_ce_ici_bb(matfile)        % a specific results file
    thisdir = fileparts(mfilename('fullpath'));
    if nargin < 1 || isempty(matfile)
        matfile = fullfile(thisdir, 'results', 'mcrb_ce_ici_results.mat');
    end
    assert(isfile(matfile), 'Results file not found: %s\nRun mcrb_ce_ici first.', matfile);
    S = load(matfile);

    SNR = S.SNR_dbs;  nSNR = numel(SNR);  nC = size(S.configs,1);
    it = S.idx_tau;   in = S.idx_nu;
    names = S.configs(:,1);
    col = lines(nC);

    % ---- decompose ---------------------------------------------------------
    [tF_full,tF_noBB,tF_BB,tF_crb] = deal(zeros(nSNR,nC));   % delay  (tau block)
    [nF_full,nF_noBB,nF_BB,nF_crb] = deal(zeros(nSNR,nC));   % Doppler (nu block)
    for ci = 1:nC
        for si = 1:nSNR
            N0 = 1/db2pow(SNR(si));
            gT = N0*sum(real(diag(S.G_bar (it,it,si,ci))));
            fT =    sum(real(diag(S.F_bar (it,it,si,ci))));
            bT =    sum(real(diag(S.BB_bar(it,it,si,ci))));
            cT = N0*sum(real(diag(S.Gc_bar(it,it,si,ci))));   % genie CRB (Xhat = X)
            tF_full(si,ci)=sqrt(gT+fT+bT);  tF_noBB(si,ci)=sqrt(gT+fT);
            tF_BB(si,ci)=sqrt(bT);          tF_crb(si,ci)=sqrt(cT);
            gN = N0*sum(real(diag(S.G_bar (in,in,si,ci))));
            fN =    sum(real(diag(S.F_bar (in,in,si,ci))));
            bN =    sum(real(diag(S.BB_bar(in,in,si,ci))));
            cN = N0*sum(real(diag(S.Gc_bar(in,in,si,ci))));
            nF_full(si,ci)=sqrt(gN+fN+bN);  nF_noBB(si,ci)=sqrt(gN+fN);
            nF_BB(si,ci)=sqrt(bN);          nF_crb(si,ci)=sqrt(cN);
        end
    end

    % ---- delay figure ------------------------------------------------------
    hf_tau = figure('Position',[100 100 760 540],'Color','w'); hold on; grid on; box on;
    lg = {};
    for ci = 1:nC
        plot(SNR, nz(tF_full(:,ci))*1e9, '-o', 'Color',col(ci,:),'LineWidth',2,'MarkerSize',6);
        plot(SNR, nz(tF_noBB(:,ci))*1e9, '--', 'Color',col(ci,:),'LineWidth',1.6);
        plot(SNR, nz(tF_crb(:,ci)) *1e9, ':',  'Color',col(ci,:),'LineWidth',1.8);
        lg = [lg, {sprintf('%s: MCRB (with BB)',names{ci}), ...
                   sprintf('%s: without BB (N_0G+F)',names{ci}), ...
                   sprintf('%s: CRB (no mismatch)',names{ci})}]; %#ok<AGROW>
    end
    set(gca,'YScale','log'); xlabel('SNR [dB]'); ylabel('delay RMSE bound [ns]');
    title(sprintf('%s ICI — delay: MCRB / without-BB / CRB  (P=%d, %d trials)', ...
                  S.evaname, S.P, S.nTrials));
    legend(lg,'Location','southwest','FontSize',8); set(gca,'FontSize',11);
    if exist('theme','file'), theme(hf_tau,"light"); end

    % ---- Doppler figure ----------------------------------------------------
    hf_nu = figure('Position',[100 100 760 540],'Color','w'); hold on; grid on; box on;
    for ci = 1:nC
        plot(SNR, nz(nF_full(:,ci)), '-o', 'Color',col(ci,:),'LineWidth',2,'MarkerSize',6);
        plot(SNR, nz(nF_noBB(:,ci)), '--', 'Color',col(ci,:),'LineWidth',1.6);
        plot(SNR, nz(nF_crb(:,ci)),  ':',  'Color',col(ci,:),'LineWidth',1.8);
    end
    set(gca,'YScale','log'); xlabel('SNR [dB]'); ylabel('Doppler RMSE bound [Hz]');
    title(sprintf('%s ICI — Doppler: MCRB / without-BB / CRB  (P=%d, %d trials)', ...
                  S.evaname, S.P, S.nTrials));
    legend(lg,'Location','southwest','FontSize',8); set(gca,'FontSize',11);
    if exist('theme','file'), theme(hf_nu,"light"); end

    % ---- save + recap ------------------------------------------------------
    resdir = fileparts(matfile);
    p1 = fullfile(resdir,'mcrb_ce_ici_bbsplit_delay.png');
    p2 = fullfile(resdir,'mcrb_ce_ici_bbsplit_doppler.png');
    exportgraphics(hf_tau, p1, 'Resolution',150, 'BackgroundColor','white');
    exportgraphics(hf_nu,  p2, 'Resolution',150, 'BackgroundColor','white');

    fprintf('\ndelay RMSE [ns]:  full | withoutBB | BB alone   (BB share of the bound)\n');
    for ci = 1:nC
        fprintf(' %s\n', names{ci});
        for si = 1:nSNR
            sh = 100*tF_BB(si,ci)^2 / max(tF_full(si,ci)^2, realmin);
            fprintf('   SNR=%2d | %9.3e | %9.3e | %9.3e   (%5.1f%%)\n', SNR(si), ...
                    tF_full(si,ci)*1e9, tF_noBB(si,ci)*1e9, tF_BB(si,ci)*1e9, sh);
        end
    end
    fprintf('\nSaved %s\n      %s\n', p1, p2);
end

function y = nz(y)
% Zeros cannot be drawn on a log axis -> mask them out (BB is exactly 0 once the
% detector makes no symbol errors).
    y(y == 0) = NaN;
end
