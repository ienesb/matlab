function hf = plot_mcrb_real_vs_bsc(matfile)
% PLOT_MCRB_REAL_VS_BSC  Does the MCRB depend only on (SNR, error rate), or on
% the error PATTERN?  Overlays, per waveform, the sensing bound computed from:
%   solid  : the REALISTIC detector's error pattern (clustered in deep fades /
%            MP failure modes)   -- mcrb_ce_ici real arm
%   dashed : a BSC-MATCHED pattern with the SAME per-trial error COUNT but
%            i.i.d. uniform placement -- the compute_mcrb_mc abstraction
% Same channel, data, SNR and symbol-error rate; only the placement differs.
% If solid == dashed, the abstract BSC model is justified and "same SER -> same
% MCRB"; where they separate, error CLUSTERING matters.
%
%   plot_mcrb_real_vs_bsc                 % default results file
%   plot_mcrb_real_vs_bsc(matfile)
    thisdir = fileparts(mfilename('fullpath'));
    if nargin < 1 || isempty(matfile)
        matfile = fullfile(thisdir, 'results', 'mcrb_ce_ici_results.mat');
    end
    assert(isfile(matfile), 'Results file not found: %s', matfile);
    S = load(matfile);
    assert(isfield(S,'G_bar_b'), 'This results file has no BSC-matched arm (run mcrb_ce_ici with do_bsc=true).');

    SNR = S.SNR_dbs;  nSNR = numel(SNR);  nC = size(S.configs,1);
    it = S.idx_tau;   in = S.idx_nu;   names = S.configs(:,1);
    col = lines(nC);

    [rt_r,rt_b,rn_r,rn_b] = deal(zeros(nSNR,nC));
    for ci = 1:nC
        for si = 1:nSNR
            N0 = 1/db2pow(SNR(si));
            LBr = N0*S.G_bar  (:,:,si,ci) + S.F_bar  (:,:,si,ci) + S.BB_bar  (:,:,si,ci);
            LBb = N0*S.G_bar_b(:,:,si,ci) + S.F_bar_b(:,:,si,ci) + S.BB_bar_b(:,:,si,ci);
            rt_r(si,ci)=sqrt(sum(real(diag(LBr(it,it)))));  rt_b(si,ci)=sqrt(sum(real(diag(LBb(it,it)))));
            rn_r(si,ci)=sqrt(sum(real(diag(LBr(in,in)))));  rn_b(si,ci)=sqrt(sum(real(diag(LBb(in,in)))));
        end
    end

    hf = figure('Position',[100 90 820 780],'Color','w');
    lg = {};
    for ci = 1:nC
        lg = [lg, {sprintf('%s: real detector',names{ci}), sprintf('%s: BSC-matched',names{ci})}]; %#ok<AGROW>
    end
    % delay
    subplot(2,1,1); hold on; grid on; box on;
    for ci = 1:nC
        plot(SNR, nz(rt_r(:,ci))*1e9, '-o','Color',col(ci,:),'LineWidth',2,'MarkerSize',5);
        plot(SNR, nz(rt_b(:,ci))*1e9, '--s','Color',col(ci,:),'LineWidth',1.6,'MarkerSize',5);
    end
    set(gca,'YScale','log'); ylabel('delay RMSE bound [ns]');
    title(sprintf('%s ICI — MCRB: realistic vs BSC-matched error pattern  (same SER, P=%d, %d trials)', ...
                  S.evaname, S.P, S.nTrials));
    legend(lg,'Location','southwest','FontSize',8); set(gca,'FontSize',10);
    % Doppler
    subplot(2,1,2); hold on; grid on; box on;
    for ci = 1:nC
        plot(SNR, nz(rn_r(:,ci)), '-o','Color',col(ci,:),'LineWidth',2,'MarkerSize',5);
        plot(SNR, nz(rn_b(:,ci)), '--s','Color',col(ci,:),'LineWidth',1.6,'MarkerSize',5);
    end
    set(gca,'YScale','log'); xlabel('SNR [dB]'); ylabel('Doppler RMSE bound [Hz]');
    legend(lg,'Location','southwest','FontSize',8); set(gca,'FontSize',10);

    if exist('theme','file'), theme(hf,"light"); end
    png = fullfile(fileparts(matfile),'mcrb_ce_ici_real_vs_bsc.png');
    exportgraphics(hf, png, 'Resolution',150, 'BackgroundColor','white');

    fprintf('\ndelay RMSE [ns]  real | BSC-matched | ratio (real/BSC)   [SER matched per trial]\n');
    for ci = 1:nC
        fprintf(' %s\n', names{ci});
        for si = 1:nSNR
            r = rt_r(si,ci)/max(rt_b(si,ci),realmin);
            fprintf('   SNR=%2d | SER=%.2e | %9.3e | %9.3e | %5.2f%s\n', SNR(si), ...
                    S.ser_data(si,ci), rt_r(si,ci)*1e9, rt_b(si,ci)*1e9, r, tick(abs(r-1)>0.15));
        end
    end
    fprintf('\nSaved %s\n', png);
end

function y = nz(y),  y(y==0) = NaN;  end
function s = tick(b),  if b, s = '  <-- differs'; else, s = ''; end,  end
