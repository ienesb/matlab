% plot_mcrb_paper_pt.m
% --------------------------------------------------------------------------
% Paper-scenario PSEUDO-TRUE bound from compute_mcrb_mc.m run with
% scenario='paper': trace RANGE bound over the delay-only targets for the
% precoded OFDM / SC / OTFS waveforms, vs SNR and vs the 16-PSK detection-
% error rate p.  Same as plot_mcrb_paper.m but the bound now carries the
% pseudo-true bias term BB (here eta_0 != eta_bar, so BB > 0).
%
% Naming convention:
%   MCRB = A^-1 B A^-1        = sigma^2 * G + F   (noise term + score-bias floor)
%   BB   = (eta_bar - eta_0)(eta_bar - eta_0)^T   (pseudo-true parameter bias)
%   LB   = MCRB + BB          = sigma^2 * G + F + BB
%
%   Figs 1-2 : LB   (the total, sigma^2*G + F + BB) vs SNR and vs p
%   Fig 3    : MCRB (the A^-1 B A^-1 term, sigma^2*G + F, WITHOUT BB)
%   Fig 4    : BB   (the pseudo-true bias term alone, SNR-independent)
%   tr{LB_tau} = sum_p LB(tau_p,tau_p) ,  range = (c/2)^2 * tr{LB_tau}
%   SNR := |alpha_1|^2 / sigma^2 .
%
% The per-trial tau-traces (trG_all/trF_all/trBB_all) are averaged with a robust
% MEDIAN.  The OFDM pilot-only CRB (known comb pilots, no mismatch) is overlaid.
% --------------------------------------------------------------------------
clear; clc;
close all;
thisdir = fileparts(mfilename('fullpath'));
resdir  = fullfile(thisdir, 'results');

% Which run to plot (<-- switch here):
%   'compute_mcrb_mc_paper_results.mat'    delay-only  (M=1, no ICI)
%   'compute_mcrb_mc_paperv2_results.mat'  delay-Doppler (M=16, ICI) -> also
%                                          emits the Doppler (nu) figures
resultfile = 'compute_mcrb_mc_paperv2_results.mat';
S = load(fullfile(resdir, resultfile));

c_light = 3e8;  u_scale = (c_light/2)^2;      % tau [s^2] -> range [m^2]
it = S.idx_tau;
% Doppler figures only when nu was actually estimated (paperv2).
do_nu = isfield(S,'trGn_all') && isfield(S,'estimate_doppler') && S.estimate_doppler;
nC = size(S.configs,1);
names = S.configs(:,1);
cols  = [0.85 0.33 0.10;    % OFDM  (red-orange)
         0.00 0.45 0.74;    % SC    (blue)
         0.93 0.69 0.13];   % OTFS  (yellow)

% Range trace pseudo-true LB at (sigma^2, config ci, p-index ip): robust
% MEDIAN over trials of sigma^2*trG + trF + trBB  (= MCRB + BB).
statname = 'median';
trR = @(sig2, ci, ip) u_scale * median( ...
        sig2*S.trG_all(:,ip,ci) + S.trF_all(:,ip,ci) + S.trBB_all(:,ip,ci));

% Pilot-only OFDM range CRB: known comb pilots -> no mismatch -> pure CRB.
has_pilot = isfield(S,'trGpil') && isfinite(S.trGpil);
pilotR    = @(sig2) u_scale * sig2 * S.trGpil;
if has_pilot, pilotlab = sprintf('OFDM pilot-only (%.0f%%)', 100*S.pilot_frac); end

% ===== Console: ranking at p = 0 (pure CRB) =============================
fprintf('\nPaper pseudo-true LB (%s, %s, P=%d, %d trials, %s over trials)\n', ...
        S.mod_type, S.chan_desc, S.P, S.nTrials, statname);
fprintf('p=0 trace range CRB [m^2] at SNR=0 dB (sigma^2=1):\n');
[~,ip0] = min(abs(S.p_list - 0));
for ci = 1:nC
    fprintf('  %-5s  %.4f  (+%.2f dB vs OFDM)\n', names{ci}, ...
            trR(1,ci,ip0), 10*log10(trR(1,ci,ip0)/trR(1,1,ip0)));
end

% ===== Figure 1: trace range LB vs SNR =================================
snr_db  = -10:1:30;
sigma2  = S.alpha_ref2 ./ 10.^(snr_db/10);
[~,ipm] = min(abs(S.p_list - 0.10));            % mid-p for the dashed curves
figure('Name','Paper pseudo-true LB vs SNR','Position',[200 200 780 520],'Color','w');
h = gobjects(nC,1);
for ci = 1:nC
    yA = arrayfun(@(s) trR(s,ci,ip0), sigma2);   % p = 0
    yB = arrayfun(@(s) trR(s,ci,ipm), sigma2);   % p = p_mid
    h(ci) = semilogy(snr_db, yA, '-',  'Color',cols(ci,:), 'LineWidth',1.9); hold on;
             semilogy(snr_db, yB, '--', 'Color',cols(ci,:), 'LineWidth',1.4);
end
if has_pilot
    hp = semilogy(snr_db, arrayfun(pilotR, sigma2), 'k-.', 'LineWidth',2.2);
    yd = arrayfun(@(s) trR(s,1,ipm), sigma2);
    xi = find(arrayfun(pilotR,sigma2) < yd, 1);
    if ~isempty(xi)
        fprintf('vs-SNR: pilot-only (%.0f%%) beats data-aided OFDM (p=%.2f) above SNR ~%d dB\n', ...
                100*S.pilot_frac, S.p_list(ipm), snr_db(xi));
    end
end
hold off; grid on; box on;
xlabel('SNR = |\alpha_1|^2 / \sigma^2  [dB]');
ylabel('tr\{LB_R\} = (c/2)^2 \Sigma_p LB(\tau_p)  [m^2]');
title(sprintf('Paper pseudo-true LB vs SNR (%s, %s): solid p=0 (CRB), dashed p=%.2f', ...
              S.mod_type, statname, S.p_list(ipm)));
if has_pilot, lh=[h;hp]; ln=[names; {pilotlab}]; else, lh=h; ln=names; end
legend(lh, ln, 'Location','southwest'); theme(gcf,"light"); set(gca,'FontSize',11);
exportgraphics(gcf, fullfile(resdir,'mcrb_paper_pt_vs_snr.png'), ...
               'Resolution',150, 'BackgroundColor','white');

% ===== Figure 2: trace range LB vs p at fixed SNR ======================
for snr_fix = [0 20]
    sig2f = S.alpha_ref2 / 10^(snr_fix/10);
    figure('Name',sprintf('Paper pseudo-true LB vs p (SNR=%d dB)',snr_fix), ...
           'Position',[240 220 780 520],'Color','w');
    for ci = 1:nC
        yp = arrayfun(@(ip) trR(sig2f,ci,ip), 1:numel(S.p_list));
        semilogy(S.p_list, yp, 'o-', 'Color',cols(ci,:), ...
                 'LineWidth',1.9, 'MarkerSize',5); hold on;
    end
    if has_pilot
        semilogy(S.p_list, pilotR(sig2f)*ones(size(S.p_list)), 'k-.', 'LineWidth',2.2);
        yd = arrayfun(@(ip) trR(sig2f,1,ip), 1:numel(S.p_list));   % OFDM data-aided
        xc = find(yd > pilotR(sig2f), 1);
        if ~isempty(xc)
            fprintf('  SNR=%2d dB: data-aided OFDM exceeds pilot-only at p>=%.3f -> pilots win above this BER\n', ...
                    snr_fix, S.p_list(xc));
        else
            fprintf('  SNR=%2d dB: data-aided OFDM below pilot-only for all p<=%.2f (data always better)\n', ...
                    snr_fix, S.p_list(end));
        end
    end
    hold off; grid on; box on;
    xlabel('16-PSK detection-error rate  p');
    ylabel('tr\{LB_R\} = (c/2)^2 \Sigma_p LB(\tau_p)  [m^2]');
    title(sprintf('Paper pseudo-true LB (%s) vs p at SNR=%d dB  (%s, P=%d)', ...
                  statname, snr_fix, S.mod_type, S.P));
    if has_pilot, ln=[names; {pilotlab}]; else, ln=names; end
    legend(ln, 'Location','northwest'); theme(gcf,"light"); set(gca,'FontSize',11);
    exportgraphics(gcf, fullfile(resdir,sprintf('mcrb_paper_pt_vs_p_snr%02d.png',snr_fix)), ...
                   'Resolution',150, 'BackgroundColor','white');
end

% ===== Figure 3: MCRB = A^-1 B A^-1 (= sigma^2*G + F, WITHOUT the BB bias) vs p =
% The misspecified-estimator MSE from its own covariance + score bias, without
% the parameter-offset (BB) term.  Shown at a fixed (high) SNR where sigma^2*G is
% small, so the MCRB is dominated by the score-bias floor F.
snr_ab = 20;  sig2ab = S.alpha_ref2 / 10^(snr_ab/10);
trAB = @(ci,ip) u_scale * median(sig2ab*S.trG_all(:,ip,ci) + S.trF_all(:,ip,ci));
figure('Name','Paper pseudo-true MCRB vs p','Position',[260 200 780 520],'Color','w');
for ci = 1:nC
    yp = arrayfun(@(ip) trAB(ci,ip), 1:numel(S.p_list));
    semilogy(S.p_list, yp, 'o-', 'Color',cols(ci,:), 'LineWidth',1.9, 'MarkerSize',5); hold on;
end
hold off; grid on; box on;
xlabel('16-PSK detection-error rate  p');
ylabel('tr\{MCRB\}_\tau = (c/2)^2 \Sigma_p (\sigma^2 G + F)(\tau_p)  [m^2]');
title(sprintf('Paper pseudo-true: MCRB = A^{-1}BA^{-1} (\\sigma^2 G + F) vs p at SNR=%d dB  (%s, P=%d)', ...
              snr_ab, S.mod_type, S.P));
legend(names, 'Location','northwest'); theme(gcf,"light"); set(gca,'FontSize',11);
exportgraphics(gcf, fullfile(resdir,'mcrb_paper_pt_ABA_vs_p.png'), ...
               'Resolution',150, 'BackgroundColor','white');

% ===== Figure 4: pseudo-true bias term BB = (eta_bar - eta_0)(.)^T vs p =======
% Parameter-offset floor from the pseudo-true drift; SNR-independent.  Zero at
% p=0 (eta_0 = eta_bar) and grows as detection errors bias the pseudo-true delay.
trBBf = @(ci,ip) u_scale * median(S.trBB_all(:,ip,ci));
figure('Name','Paper pseudo-true BB vs p','Position',[300 180 780 520],'Color','w');
for ci = 1:nC
    yp = arrayfun(@(ip) trBBf(ci,ip), 1:numel(S.p_list));
    semilogy(S.p_list, yp, 'o-', 'Color',cols(ci,:), 'LineWidth',1.9, 'MarkerSize',5); hold on;
end
hold off; grid on; box on;
xlabel('16-PSK detection-error rate  p');
ylabel('tr\{BB\}_\tau = (c/2)^2 \Sigma_p (\tau_{bar}-\tau_0)^2(\tau_p)  [m^2]');
title(sprintf('Paper pseudo-true: bias term BB = (\\eta_{bar}-\\eta_0)(\\eta_{bar}-\\eta_0)^T vs p (SNR-indep.)  (%s, P=%d)', ...
              S.mod_type, S.P));
legend(names, 'Location','northwest'); theme(gcf,"light"); set(gca,'FontSize',11);
exportgraphics(gcf, fullfile(resdir,'mcrb_paper_pt_BB_vs_p.png'), ...
               'Resolution',150, 'BackgroundColor','white');

% ===== Doppler (nu) figures -- paperv2 only ============================
% Trace of the pseudo-true Doppler LB [Hz^2]: median of sig2*trGn+trFn+trBBn.
if do_nu
    trNU    = @(sig2, ci, ip) median(sig2*S.trGn_all(:,ip,ci) + S.trFn_all(:,ip,ci) ...
                                     + S.trBBn_all(:,ip,ci));
    has_pn  = isfield(S,'trGpiln') && isfinite(S.trGpiln);
    pilotNU = @(sig2) sig2 * S.trGpiln;

    fprintf('\np=0 trace Doppler CRB [Hz^2] at SNR=0 dB (sigma^2=1):\n');
    for ci = 1:nC
        fprintf('  %-5s  %.4e  (+%.2f dB vs OFDM)\n', names{ci}, ...
                trNU(1,ci,ip0), 10*log10(trNU(1,ci,ip0)/trNU(1,1,ip0)));
    end

    % ---- vs SNR ----
    figure('Name','Paper pseudo-true LB (nu) vs SNR','Position',[280 180 780 520],'Color','w');
    hn = gobjects(nC,1);
    for ci = 1:nC
        hn(ci) = semilogy(snr_db, arrayfun(@(s) trNU(s,ci,ip0), sigma2), '-', ...
                          'Color',cols(ci,:),'LineWidth',1.9); hold on;
                 semilogy(snr_db, arrayfun(@(s) trNU(s,ci,ipm), sigma2), '--', ...
                          'Color',cols(ci,:),'LineWidth',1.4);
    end
    if has_pn, hpn = semilogy(snr_db, arrayfun(pilotNU, sigma2), 'k-.','LineWidth',2.2); end
    hold off; grid on; box on;
    xlabel('SNR = |\alpha_1|^2 / \sigma^2  [dB]');
    ylabel('tr\{LB_\nu\} = \Sigma_p LB(\nu_p)  [Hz^2]');
    title(sprintf('Paper pseudo-true LB (Doppler) vs SNR (%s, %s): solid p=0 (CRB), dashed p=%.2f', ...
                  S.mod_type, statname, S.p_list(ipm)));
    if has_pn, lh=[hn;hpn]; ln=[names; {pilotlab}]; else, lh=hn; ln=names; end
    legend(lh, ln, 'Location','southwest'); theme(gcf,"light"); set(gca,'FontSize',11);
    exportgraphics(gcf, fullfile(resdir,'mcrb_paper_pt_nu_vs_snr.png'), ...
                   'Resolution',150, 'BackgroundColor','white');

    % ---- vs p at fixed SNR ----
    for snr_fix = [0 20]
        sig2f = S.alpha_ref2 / 10^(snr_fix/10);
        figure('Name',sprintf('Paper pseudo-true LB (nu) vs p (SNR=%d dB)',snr_fix), ...
               'Position',[320 160 780 520],'Color','w');
        for ci = 1:nC
            semilogy(S.p_list, arrayfun(@(ip) trNU(sig2f,ci,ip), 1:numel(S.p_list)), ...
                     'o-','Color',cols(ci,:),'LineWidth',1.9,'MarkerSize',5); hold on;
        end
        if has_pn
            semilogy(S.p_list, pilotNU(sig2f)*ones(size(S.p_list)), 'k-.','LineWidth',2.2);
        end
        hold off; grid on; box on;
        xlabel('16-PSK detection-error rate  p');
        ylabel('tr\{LB_\nu\} = \Sigma_p LB(\nu_p)  [Hz^2]');
        title(sprintf('Paper pseudo-true LB (Doppler, %s) vs p at SNR=%d dB  (%s, P=%d)', ...
                      statname, snr_fix, S.mod_type, S.P));
        if has_pn, ln=[names; {pilotlab}]; else, ln=names; end
        legend(ln, 'Location','northwest'); theme(gcf,"light"); set(gca,'FontSize',11);
        exportgraphics(gcf, fullfile(resdir,sprintf('mcrb_paper_pt_nu_vs_p_snr%02d.png',snr_fix)), ...
                       'Resolution',150, 'BackgroundColor','white');
    end
    fprintf('Saved mcrb_paper_pt_nu_vs_snr.png, mcrb_paper_pt_nu_vs_p_snr00/20.png\n');
end

fprintf(['Saved mcrb_paper_pt_vs_snr.png, mcrb_paper_pt_vs_p_snr00.png, ' ...
         'mcrb_paper_pt_vs_p_snr20.png,\n  mcrb_paper_pt_ABA_vs_p.png, ' ...
         'mcrb_paper_pt_BB_vs_p.png in %s\n'], resdir);
