function [hf_tau, hf_nu, hf_ser] = plot_bsc_vs_ce(bscfile, cefile)
% PLOT_BSC_VS_CE  Compare the ABSTRACT-BSC MCRB (compute_mcrb_mc, scenario
% 'evace') with the REALISTIC-DETECTOR MCRB (mcrb_ce_ici) on the SAME reduced-EVA
% 512x128 channel.  The two scripts sweep different axes (BSC over p, CE over
% SNR), so the BSC bound is matched to the CE detector's symbol-error rate at
% each SNR and evaluated at the SAME noise sigma^2 = N0(SNR):
%
%   real (solid) : mcrb_ce_ici MCRB = N0*G + F + BB     (clustered detector errors)
%   BSC  (dashed): evace MCRB at N0, SER-matched p       (i.i.d. BSC errors)
%   CRB  (dotted): genie bound N0*Gc                     (no mismatch)
%
% Where solid > dashed, the abstract BSC UNDERSTATES the sensing penalty -- a
% real detector's structured errors bias the pseudo-true parameter more than
% coin-flip errors of the same rate.  QPSK BSC: SER = 1-(1-p)^2.
%
%   plot_bsc_vs_ce                       % default result files
%   plot_bsc_vs_ce(bscfile, cefile)     % specific files
    thisdir = fileparts(mfilename('fullpath'));
    if nargin < 1 || isempty(bscfile), bscfile = fullfile(thisdir,'results','compute_mcrb_mc_evace_results.mat'); end
    if nargin < 2 || isempty(cefile),  cefile  = fullfile(thisdir,'results','mcrb_ce_ici_results.mat'); end
    B = load(bscfile);   R = load(cefile);

    SNR = R.SNR_dbs(:);  nS = numel(SNR);
    itR = R.idx_tau;  inR = R.idx_nu;  itB = B.idx_tau;  inB = B.idx_nu;
    ser_bsc = 1 - (1 - B.p_list).^2;                 % QPSK BSC symbol-error rate
    ser_min = 0.5 * ser_bsc(2);                      % below this, p snaps to 0 (skip)

    names = {'OFDM ICI','CP-OTFS ICI'};  cols = [0.85 0.33 0.10; 0.00 0.45 0.74];
    nC = numel(names);
    [rt_r,rt_b,rt_c, rn_r,rn_b,rn_c] = deal(nan(nS,nC));   % delay / Doppler RMSE
    for ci = 1:nC
        cR = find(strcmp(R.configs(:,1), names{ci}), 1);
        cB = find(strcmp(B.configs(:,1), names{ci}), 1);
        for si = 1:nS
            N0 = 1/db2pow(SNR(si));  s = R.ser_data(si,cR);
            LB  = N0*R.G_bar(:,:,si,cR) + R.F_bar(:,:,si,cR) + R.BB_bar(:,:,si,cR);
            CRB = N0*R.Gc_bar(:,:,si,cR);
            rt_r(si,ci)=sqrt(sum(real(diag(LB (itR,itR)))));  rt_c(si,ci)=sqrt(sum(real(diag(CRB(itR,itR)))));
            rn_r(si,ci)=sqrt(sum(real(diag(LB (inR,inR)))));  rn_c(si,ci)=sqrt(sum(real(diag(CRB(inR,inR)))));
            if s >= ser_min                          % genuine SER match available
                [~,ip] = min(abs(ser_bsc - s));
                rt_b(si,ci) = sqrt(median(N0*B.trG_all(:,ip,cB)  + B.trF_all(:,ip,cB)  + B.trBB_all(:,ip,cB)));
                rn_b(si,ci) = sqrt(median(N0*B.trGn_all(:,ip,cB) + B.trFn_all(:,ip,cB) + B.trBBn_all(:,ip,cB)));
            end
        end
    end

    hf_tau = mkfig(SNR, rt_r, rt_b, rt_c, cols, names, 1e9, ...
        'delay RMSE bound [ns]', 'delay: realistic-detector vs BSC (SER-matched)', ...
        fullfile(fileparts(cefile),'mcrb_ce_ici_bsc_vs_ce_delay.png'), R);
    hf_nu  = mkfig(SNR, rn_r, rn_b, rn_c, cols, names, 1, ...
        'Doppler RMSE bound [Hz]', 'Doppler: realistic-detector vs BSC (SER-matched)', ...
        fullfile(fileparts(cefile),'mcrb_ce_ici_bsc_vs_ce_doppler.png'), R);

    % ---- detection SER / BER vs SNR: the error rate that DRIVES the mismatch ---
    % Also shows where the BSC comparison is valid: the dotted line is the
    % smallest SER the BSC p-grid can match; below it the dashed BSC curves above
    % simply stop (SER too small for any p on the grid).
    hf_ser = figure('Position',[160 100 820 540],'Color','w'); hold on; grid on; box on;
    lg = {};
    for ci = 1:nC
        cR = find(strcmp(R.configs(:,1), names{ci}), 1);
        semilogy(SNR, nz(R.ser_data(:,cR)), '-o', 'Color',cols(ci,:),'LineWidth',2,'MarkerSize',6);
        lg = [lg, {[names{ci} ': SER']}]; %#ok<AGROW>
        if isfield(R,'ber_data')
            semilogy(SNR, nz(R.ber_data(:,cR)), '--^', 'Color',cols(ci,:),'LineWidth',1.6,'MarkerSize',5);
            lg = [lg, {[names{ci} ': BER']}]; %#ok<AGROW>
        end
    end
    % yline(ser_min, ':', 'BSC min matchable SER', 'Color',[0.4 0.4 0.4], ...
          % 'LineWidth',1.3, 'FontSize',8, 'LabelHorizontalAlignment','left');
    set(gca,'YScale','log'); xlabel('SNR [dB]'); ylabel('detection error rate');
    title(sprintf(['reduced EVA 512x128 ICI — detection SER / BER (drives the MCRB mismatch)\n' ...
                   '(%d trials; real detector = mcrb\\_ce\\_ici; BSC matches this SER at each SNR)'], R.nTrials));
    legend(lg,'Location','southwest','FontSize',8); set(gca,'FontSize',11);
    if exist('theme','file'), theme(hf_ser,"light"); end
    exportgraphics(hf_ser, fullfile(fileparts(cefile),'mcrb_ce_ici_bsc_vs_ce_ser.png'), ...
                   'Resolution',150, 'BackgroundColor','white');
    fprintf('Saved %s\n', fullfile(fileparts(cefile),'mcrb_ce_ici_bsc_vs_ce_ser.png'));

    fprintf('\nOFDM delay RMSE [ns]  real | BSC(SER-matched) | ratio  (BSC blank where SER<%.3f)\n', ser_min);
    for si = 1:nS
        fprintf('  SNR=%2d | SER=%.2e | %9.3e | %9.3e | %s\n', SNR(si), R.ser_data(si,1), ...
                rt_r(si,1)*1e9, rt_b(si,1)*1e9, tern(isnan(rt_b(si,1)),'  --', sprintf('%.2f',rt_r(si,1)/rt_b(si,1))));
    end
end

function hf = mkfig(SNR, rr, rb, rc, cols, names, scl, ylab, ttl, png, R)
    hf = figure('Position',[120 120 820 560],'Color','w'); hold on; grid on; box on;
    lg = {};
    for ci = 1:size(rr,2)
        plot(SNR, rr(:,ci)*scl, '-o',  'Color',cols(ci,:),'LineWidth',2,'MarkerSize',6);
        plot(SNR, rb(:,ci)*scl, '--s', 'Color',cols(ci,:),'LineWidth',1.8,'MarkerSize',6);
        plot(SNR, rc(:,ci)*scl, ':',   'Color',cols(ci,:),'LineWidth',1.6);
        lg = [lg, {[names{ci} ': real (CE)'], [names{ci} ': BSC (SER-matched)'], [names{ci} ': CRB']}]; %#ok<AGROW>
    end
    set(gca,'YScale','log'); xlabel('SNR [dB]'); ylabel(ylab);
    title(sprintf('reduced EVA 512x128 ICI — %s\n(%d trials; BSC = compute\\_mcrb\\_mc evace)', ttl, R.nTrials));
    legend(lg,'Location','southwest','FontSize',8); set(gca,'FontSize',11);
    if exist('theme','file'), theme(hf,"light"); end
    exportgraphics(hf, png, 'Resolution',150, 'BackgroundColor','white');
    fprintf('Saved %s\n', png);
end

function s = tern(c,a,b), if c, s=a; else, s=b; end, end
function y = nz(y), y(y==0) = NaN; end     % mask zeros for the log axis
