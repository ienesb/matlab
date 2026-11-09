% plot_mcrb_paper.m
% --------------------------------------------------------------------------
% Paper-scenario bound from compute_mcrb_mc_true.m run with scenario='paper':
% trace RANGE bound over the 40 delay-only targets, for the precoded
% OFDM / SC / OTFS waveforms, vs SNR and vs the 16-PSK detection-error rate p.
%
% Naming convention (see also plot_mcrb_paper_pt.m):
%   MCRB = A^-1 B A^-1        = sigma^2 * G + F   (noise term + score-bias floor)
%   BB   = (eta_bar - eta_0)(eta_bar - eta_0)^T   (pseudo-true parameter bias)
%   LB   = MCRB + BB          = sigma^2 * G + F + BB
%
% This script evaluates at the TRUE parameter, so BB = 0 and LB = MCRB.  The
% plotted quantity is therefore the MCRB, affine in noise power:
%   tr{MCRB_tau} = sum_p MCRB(tau_p, tau_p)  [s^2] ,  range = (c/2)^2 * tr{MCRB_tau}
%   SNR := |alpha_1|^2 / sigma^2 .
%
%   p = 0  -> F = 0, so tr{MCRB_tau} = sigma^2 * G = the pure CRB (Fig. 1(b)).
%   p > 0  -> F is an SNR-independent floor: detection errors cap the ranging
%             accuracy no matter how high the SNR.
% --------------------------------------------------------------------------
clear; clc; close all;
thisdir = fileparts(mfilename('fullpath'));
resdir  = fullfile(thisdir, 'results');

% Which run to plot (<-- switch here):
%   'compute_mcrb_mc_true_paper_results.mat'    delay-only  (M=1, no ICI)
%   'compute_mcrb_mc_true_paperv2_results.mat'  delay-Doppler (M=16, ICI) -> also
%                                               emits the Doppler (nu) figures
resultfile = 'compute_mcrb_mc_true_paperv2_results.mat';
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

% Range trace bound at (sigma^2, config ci, p-index ip).  Prefer the per-trial
% MEDIAN (robust to the heavy-tailed spread-waveform mismatch floor); fall back
% to the averaged pieces (mean) if only the _bar matrices were saved.
if isfield(S,'trG_all') && isfield(S,'trF_all')
    statname = 'median';
    trR = @(sig2, ci, ip) u_scale * median(sig2*S.trG_all(:,ip,ci) + S.trF_all(:,ip,ci));
else
    statname = 'mean';
    gT = zeros(nC, numel(S.p_list));  fT = zeros(nC, numel(S.p_list));
    for ci = 1:nC
        for ip = 1:numel(S.p_list)
            gT(ci,ip) = sum(real(diag(S.G_bar(it,it,ip,ci))));
            fT(ci,ip) = sum(real(diag(S.F_bar(it,it,ip,ci))));
        end
    end
    trR = @(sig2, ci, ip) u_scale * (sig2.*gT(ci,ip) + fT(ci,ip));
end

% Pilot-only OFDM range CRB: known comb pilots -> no mismatch -> pure CRB (no
% floor), so it scales as u*sigma^2*trGpil and is independent of p.
has_pilot = isfield(S,'trGpil') && isfinite(S.trGpil);
pilotR    = @(sig2) u_scale * sig2 * S.trGpil;
if has_pilot, pilotlab = sprintf('OFDM pilot-only (%.0f%%)', 100*S.pilot_frac); end

% ===== Console: ranking at p = 0 (pure CRB) =============================
fprintf('\nPaper MCRB (%s, %s, P=%d, %d trials, %s over trials)\n', ...
        S.mod_type, S.chan_desc, S.P, S.nTrials, statname);
fprintf('p=0 trace range CRB [m^2] at SNR=0 dB (sigma^2=1):\n');
[~,ip0] = min(abs(S.p_list - 0));
for ci = 1:nC
    fprintf('  %-5s  %.4f  (+%.2f dB vs OFDM)\n', names{ci}, ...
            trR(1,ci,ip0), 10*log10(trR(1,ci,ip0)/trR(1,1,ip0)));
end

% ===== Figure 1: trace range MCRB vs SNR ===============================
% Solid = p 0 (pure CRB), dashed = a representative p > 0 (mismatch floor).
snr_db  = -10:1:30;
sigma2  = S.alpha_ref2 ./ 10.^(snr_db/10);
[~,ipm] = min(abs(S.p_list - 0.10));            % mid-p for the dashed curves
figure('Name','Paper MCRB vs SNR','Position',[200 200 780 520],'Color','w');
h = gobjects(nC,1);
for ci = 1:nC
    yA = arrayfun(@(s) trR(s,ci,ip0), sigma2);   % p = 0
    yB = arrayfun(@(s) trR(s,ci,ipm), sigma2);   % p = p_mid
    h(ci) = semilogy(snr_db, yA, '-',  'Color',cols(ci,:), 'LineWidth',1.9); hold on;
             semilogy(snr_db, yB, '--', 'Color',cols(ci,:), 'LineWidth',1.4);
end
if has_pilot
    hp = semilogy(snr_db, arrayfun(pilotR, sigma2), 'k-.', 'LineWidth',2.2);
    yd = arrayfun(@(s) trR(s,1,ipm), sigma2);            % OFDM data-aided at p_mid
    xi = find(arrayfun(pilotR,sigma2) < yd, 1);
    if ~isempty(xi)
        fprintf('vs-SNR: pilot-only (%.0f%%) beats data-aided OFDM (p=%.2f) above SNR ~%d dB\n', ...
                100*S.pilot_frac, S.p_list(ipm), snr_db(xi));
    end
end
hold off; grid on; box on;
xlabel('SNR = |\alpha_1|^2 / \sigma^2  [dB]');
ylabel('tr\{MCRB_R\} = (c/2)^2 \Sigma_p MCRB(\tau_p)  [m^2]');
title(sprintf('Paper MCRB vs SNR  (%s, %s): solid p=0 (CRB), dashed p=%.2f', ...
              S.mod_type, statname, S.p_list(ipm)));
if has_pilot, lh=[h;hp]; ln=[names; {pilotlab}]; else, lh=h; ln=names; end
legend(lh, ln, 'Location','southwest'); theme(gcf,"light"); set(gca,'FontSize',11);
exportgraphics(gcf, fullfile(resdir,'mcrb_paper_vs_snr.png'), ...
               'Resolution',150, 'BackgroundColor','white');

% ===== Figure 2: trace range MCRB vs p at fixed SNR ====================
for snr_fix = [0 20]
    sig2f = S.alpha_ref2 / 10^(snr_fix/10);
    figure('Name',sprintf('Paper MCRB vs p (SNR=%d dB)',snr_fix), ...
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
    ylabel('tr\{MCRB_R\} = (c/2)^2 \Sigma_p MCRB(\tau_p)  [m^2]');
    title(sprintf('Paper MCRB (%s) vs p at SNR=%d dB  (%s, P=%d)', ...
                  statname, snr_fix, S.mod_type, S.P));
    if has_pilot, ln=[names; {pilotlab}]; else, ln=names; end
    legend(ln, 'Location','northwest'); theme(gcf,"light"); set(gca,'FontSize',11);
    exportgraphics(gcf, fullfile(resdir,sprintf('mcrb_paper_vs_p_snr%02d.png',snr_fix)), ...
                   'Resolution',150, 'BackgroundColor','white');
end

% ===== Doppler (nu) figures -- paperv2 only ============================
% Trace of the Doppler MCRB [Hz^2]: same median construction on the nu block.
if do_nu
    trNU    = @(sig2, ci, ip) median(sig2*S.trGn_all(:,ip,ci) + S.trFn_all(:,ip,ci));
    has_pn  = isfield(S,'trGpiln') && isfinite(S.trGpiln);
    pilotNU = @(sig2) sig2 * S.trGpiln;

    fprintf('\np=0 trace Doppler CRB [Hz^2] at SNR=0 dB (sigma^2=1):\n');
    for ci = 1:nC
        fprintf('  %-5s  %.4e  (+%.2f dB vs OFDM)\n', names{ci}, ...
                trNU(1,ci,ip0), 10*log10(trNU(1,ci,ip0)/trNU(1,1,ip0)));
    end

    % ---- vs SNR ----
    figure('Name','Paper MCRB (nu) vs SNR','Position',[280 180 780 520],'Color','w');
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
    ylabel('tr\{MCRB_\nu\} = \Sigma_p MCRB(\nu_p)  [Hz^2]');
    title(sprintf('Paper MCRB (Doppler) vs SNR  (%s, %s): solid p=0 (CRB), dashed p=%.2f', ...
                  S.mod_type, statname, S.p_list(ipm)));
    if has_pn, lh=[hn;hpn]; ln=[names; {pilotlab}]; else, lh=hn; ln=names; end
    legend(lh, ln, 'Location','southwest'); theme(gcf,"light"); set(gca,'FontSize',11);
    exportgraphics(gcf, fullfile(resdir,'mcrb_paper_nu_vs_snr.png'), ...
                   'Resolution',150, 'BackgroundColor','white');

    % ---- vs p at fixed SNR ----
    for snr_fix = [0 20]
        sig2f = S.alpha_ref2 / 10^(snr_fix/10);
        figure('Name',sprintf('Paper MCRB (nu) vs p (SNR=%d dB)',snr_fix), ...
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
        ylabel('tr\{MCRB_\nu\} = \Sigma_p MCRB(\nu_p)  [Hz^2]');
        title(sprintf('Paper MCRB (Doppler, %s) vs p at SNR=%d dB  (%s, P=%d)', ...
                      statname, snr_fix, S.mod_type, S.P));
        if has_pn, ln=[names; {pilotlab}]; else, ln=names; end
        legend(ln, 'Location','northwest'); theme(gcf,"light"); set(gca,'FontSize',11);
        exportgraphics(gcf, fullfile(resdir,sprintf('mcrb_paper_nu_vs_p_snr%02d.png',snr_fix)), ...
                       'Resolution',150, 'BackgroundColor','white');
    end
    fprintf('Saved mcrb_paper_nu_vs_snr.png, mcrb_paper_nu_vs_p_snr00.png, mcrb_paper_nu_vs_p_snr20.png\n');
end

fprintf(['Saved mcrb_paper_vs_snr.png, mcrb_paper_vs_p_snr00.png, ' ...
         'mcrb_paper_vs_p_snr20.png in %s\n'], resdir);
