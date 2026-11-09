function [hf_tau, hf_nu] = plot_mcrb_ce_ici(matfile, include_BB)
% PLOT_MCRB_CE_ICI  Plot delay and Doppler RMSE bounds (MCRB at pseudo-true from
% realistic detection, and the genie CRB) vs SNR, for OFDM-ICI and CP-OTFS-ICI,
% in TWO SEPARATE figures (delay, Doppler).  Curves are recomputed from the stored
% affine pieces  LB = N0*G + F + BB  so the pseudo-true bias term BB can be toggled.
%
%   plot_mcrb_ce_ici                          % full MCRB (LB = N0*G + F + BB)
%   plot_mcrb_ce_ici(matfile)                 % from a specific results file
%   plot_mcrb_ce_ici(matfile, false)          % covariance part only (drop BB)
%   [hTau,hNu] = plot_mcrb_ce_ici(...)        % handles to the delay/Doppler figures
    thisdir = fileparts(mfilename('fullpath'));
    if nargin < 1 || isempty(matfile)
        matfile = fullfile(thisdir, 'results', 'mcrb_ce_ici_results_2000.mat');
    end
    if nargin < 2 || isempty(include_BB), include_BB = true; end
    assert(isfile(matfile), 'Results file not found: %s\nRun mcrb_ce_ici first.', matfile);
    S = load(matfile);

    SNR = S.SNR_dbs;  nC = size(S.configs,1);  nSNR = numel(SNR);
    it = S.idx_tau;   in = S.idx_nu;
    names = S.configs(:,1);
    col = lines(nC);

    % ---- recompute bounds from stored matrices (BB optionally excluded) ----
    rt_mcrb = zeros(nSNR,nC);  rn_mcrb = zeros(nSNR,nC);
    rt_crb  = zeros(nSNR,nC);  rn_crb  = zeros(nSNR,nC);
    for ci = 1:nC
        for si = 1:nSNR
            N0 = 1/db2pow(SNR(si));
            LB  = N0*S.G_bar(:,:,si,ci) + S.F_bar(:,:,si,ci) + include_BB*S.BB_bar(:,:,si,ci);
            CRB = N0*S.Gc_bar(:,:,si,ci);
            rt_mcrb(si,ci)=sqrt(sum(real(diag(LB(it,it)))));
            rn_mcrb(si,ci)=sqrt(sum(real(diag(LB(in,in)))));
            rt_crb(si,ci) =sqrt(sum(real(diag(CRB(it,it)))));
            rn_crb(si,ci) =sqrt(sum(real(diag(CRB(in,in)))));
        end
    end
    boundlabel = 'LB';  bbsuffix = '';
    if ~include_BB, boundlabel = 'MCRB\\{BB}'; bbsuffix = ' (BB excluded)'; end
    lg = cell(1,2*nC);
    for ci=1:nC
        lg{2*ci-1}=sprintf('%s (%s)',names{ci},boundlabel);
        lg{2*ci}  =sprintf('%s (CRB)',names{ci});
    end

    % ---- delay figure ----
    hf_tau = figure('Position',[100 100 640 480]); hold on; grid on;
    for ci = 1:nC
        semilogy(SNR, rt_mcrb(:,ci)*1e9, '-o', 'Color', col(ci,:), 'LineWidth', 2, 'MarkerSize',7);
        semilogy(SNR, rt_crb(:,ci) *1e9, '--', 'Color', col(ci,:), 'LineWidth', 1.5);
    end
    set(gca,'YScale','log'); xlabel('SNR [dB]'); ylabel('delay RMSE bound [ns]');
    legend(lg,'Location','southwest');
    title(sprintf('%s ICI — delay estimation%s  (P=%d, %d trials)', S.evaname, bbsuffix, S.P, S.nTrials));

    % ---- Doppler figure ----
    hf_nu = figure('Position',[100 100 640 480]); hold on; grid on;
    for ci = 1:nC
        semilogy(SNR, rn_mcrb(:,ci), '-o', 'Color', col(ci,:), 'LineWidth', 2, 'MarkerSize',7);
        semilogy(SNR, rn_crb(:,ci),  '--', 'Color', col(ci,:), 'LineWidth', 1.5);
    end
    set(gca,'YScale','log'); xlabel('SNR [dB]'); ylabel('Doppler RMSE bound [Hz]');
    legend(lg,'Location','southwest');
    title(sprintf('%s ICI — Doppler estimation%s  (P=%d, %d trials)', S.evaname, bbsuffix, S.P, S.nTrials));

    % ---- save (so headless -batch runs actually persist the figures) --------
    resdir = fileparts(matfile);
    tag = '';  if ~include_BB, tag = 'noBB_'; end
    for h = [hf_tau, hf_nu]
        if exist('theme','file'), theme(h,"light"); end
    end
    png_tau = fullfile(resdir, sprintf('mcrb_ce_ici_%sdelay.png',   tag));
    png_nu  = fullfile(resdir, sprintf('mcrb_ce_ici_%sdoppler.png', tag));
    exportgraphics(hf_tau, png_tau, 'Resolution',150, 'BackgroundColor','white');
    exportgraphics(hf_nu,  png_nu,  'Resolution',150, 'BackgroundColor','white');
    fprintf('Saved %s\n      %s\n', png_tau, png_nu);

    % ---- console recap ----
    fprintf('Loaded %s  [%s]\n', matfile, boundlabel);
    fprintf('%6s', 'SNR');
    for ci=1:nC, fprintf(' | %-11s SER  RMSEtau[ns]', names{ci}); end
    fprintf('\n');
    for si = 1:nSNR
        fprintf('%5d ', SNR(si));
        for ci=1:nC, fprintf(' | %.2e  %10.3e', S.ber_data(si,ci), rt_mcrb(si,ci)*1e9); end
        fprintf('\n');
    end
end
