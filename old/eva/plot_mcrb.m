% plot_mcrb.m
% --------------------------------------------------------------------------
% Average pseudo-true lower bound on the reduced EVA channel (multi-target),
% trace over the delay (tau) block, for {OFDM, CP-OTFS} x {no-ICI, ICI}.
% Reads compute_mcrb_mc.m output (averaged affine pieces G_bar, F_bar, BB_bar):
%
%   LB(sigma^2) = sigma^2 * G + F + BB,   BB = (eta_bar - eta_0)(eta_bar - eta_0)^T
%   tr{LB_tau}  = sum_p LB(tau_p, tau_p)    [s^2],
%   SNR := |alpha_1|^2 / sigma^2.
%
% Two figures, four curves each:
%   Fig 1: fixed p = 0.1,    x-axis = SNR [dB]
%   Fig 2: fixed SNR = 0 dB, x-axis = p
% --------------------------------------------------------------------------
clear; clc;
close all;
thisdir = fileparts(mfilename('fullpath'));
resdir  = fullfile(thisdir, 'results');

S = load(fullfile(resdir, 'compute_mcrb_mc_true_results.mat'), ...
         'G_bar','F_bar','BB_bar','eta0_bar','eta_bar','configs','p_list', ...
         'alpha_ref2','idx_tau','P','nTrials');
% S = load(fullfile(resdir, 'compute_mcrb_mc_true_results.mat'), ...
%          'G_bar','F_bar','BB_bar','eta0_bar','eta_bar','configs','p_list', ...
%          'alpha_ref2','idx_tau','P','nTrials');
nC = size(S.configs,1);  

it = S.idx_tau;
% it = 3;

% Scalar tau-traces of the three averaged pieces, per (config, p).
gT = zeros(nC, numel(S.p_list));   % CRB-like (scales with sigma^2)
fT = zeros(nC, numel(S.p_list));   % score-bias floor
bT = zeros(nC, numel(S.p_list));   % pseudo-true bias floor
for ci = 1:nC
    for ip = 1:numel(S.p_list)
        gT(ci,ip) = sum(real(diag(S.G_bar (it,it,ip,ci))));
        fT(ci,ip) = sum(real(diag(S.F_bar (it,it,ip,ci))));
        bT(ci,ip) = sum(real(diag(S.BB_bar(it,it,ip,ci))));
    end
end
trLB    = @(sig2, ci, ip) sig2.*gT(ci,ip) + fT(ci,ip) + bT(ci,ip);   % tr{LB_tau}
trLB_AB = @(sig2, ci, ip) sig2.*gT(ci,ip) + fT(ci,ip);               % tr{A^-1 B A^-1} (no BB bias)
trLB_BB = @(sig2, ci, ip) bT(ci,ip) + 0.*sig2;                       % tr{BB}_tau only (SNR-indep.)

styles = {'-o','-s','--^','--d'};
colors = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.47 0.67 0.19; 0.49 0.18 0.56];

% ===== Console: compare the four configs at SNR = -40 dB ================
snr_cmp  = -40;                                  % SNR [dB] for the comparison
sig2_cmp = S.alpha_ref2 / 10^(snr_cmp/10);       % matching sigma^2
fprintf('\ntr{LB_tau} at SNR = %d dB  [s^2]   (reduced EVA, P=%d, %d trials)\n', ...
        snr_cmp, S.P, S.nTrials);
fprintf('%-8s', 'p');  fprintf('%16s', S.configs{:,1});  fprintf('\n');
for ip = 1:numel(S.p_list)
    fprintf('%-8.1f', S.p_list(ip));
    for ci = 1:nC
        fprintf('%16.6e', trLB(sig2_cmp, ci, ip));
    end
    fprintf('\n');
end
fprintf('\n');

% Ranking (lowest -> highest tr{LB_tau}) per p.
fprintf('Ranking lowest -> highest tr{LB_tau} at SNR = %d dB\n', snr_cmp);
for ip = 1:numel(S.p_list)
    vals = arrayfun(@(ci) trLB(sig2_cmp, ci, ip), 1:nC);
    [vs, ord] = sort(vals, 'ascend');
    fprintf('p = %.1f :  ', S.p_list(ip));
    for r = 1:nC
        fprintf('%s (%.3e)', S.configs{ord(r),1}, vs(r));
        if r < nC, fprintf('  <  '); end
    end
    fprintf('\n');
end
fprintf('\n');

% ===== Figure 1: tr{LB_tau} vs SNR at fixed p = 0.1 ======================
for p_fix = 0:0.1:0.5
    ip = find(abs(S.p_list - p_fix) < 1e-9, 1);
    assert(~isempty(ip), 'p = %.2f not in p_list', p_fix);
    snr_db  = -40:2:30;
    sigma2  = S.alpha_ref2 ./ 10.^(snr_db/10);
    figure('Name','Avg LB vs SNR (fixed p)','Position',[200 220 780 520]);
    for ci = 1:nC
        semilogy(snr_db, trLB(sigma2, ci, ip), styles{ci}, 'Color', colors(ci,:), ...
                 'LineWidth', 1.8, 'MarkerSize', 6); hold on;
    end
    hold off; grid on;
    xlabel('SNR = |\alpha_1|^2 / \sigma^2  [dB]');
    ylabel('tr\{LB_\tau\} = \Sigma_p LB(\tau_p)  [s^2]');
    title(sprintf('Average pseudo-true bound vs SNR  (reduced EVA, P=%d, p=%.1f, %d trials)', ...
                  S.P, p_fix, S.nTrials));
    legend(S.configs(:,1), 'Location','southwest');  set(gca,'FontSize',11);
    png1 = fullfile(resdir, sprintf('mcrb_lb_vs_snr_p%02d.png', round(100*p_fix)));
    exportgraphics(gcf, png1, 'Resolution',150);
    fprintf('Saved %s\n', png1);
end


% ===== Figure 2: tr{LB_tau} vs p at fixed SNR = 0 dB ====================
snr_fix = -40;  sig2_fix = S.alpha_ref2 / 10^(snr_fix/10);
figure('Name','Avg LB vs p (fixed SNR)','Position',[280 180 780 520]);
for ci = 1:nC
    semilogy(S.p_list, trLB(sig2_fix, ci, 1:numel(S.p_list)), styles{ci}, ...
             'Color', colors(ci,:), 'LineWidth', 1.8, 'MarkerSize', 6); hold on;
end
hold off; grid on;
xlabel('BSC error rate  p');
ylabel('tr\{LB_\tau\} = \Sigma_p LB(\tau_p)  [s^2]');
title(sprintf('Average pseudo-true bound vs p  (reduced EVA, P=%d, SNR=%d dB, %d trials)', ...
              S.P, snr_fix, S.nTrials));
legend(S.configs(:,1), 'Location','northwest');  set(gca,'FontSize',11);
png2 = fullfile(resdir,'mcrb_lb_vs_p.png');  exportgraphics(gcf, png2, 'Resolution',150);

% ===== Figure 2: tr{LB_tau} vs p at fixed SNR = 0 dB ====================
snr_fix = 30;  sig2_fix = S.alpha_ref2 / 10^(snr_fix/10);
figure('Name','Avg LB vs p (fixed SNR)','Position',[280 180 780 520]);
for ci = 1:nC
    semilogy(S.p_list, trLB(sig2_fix, ci, 1:numel(S.p_list)), styles{ci}, ...
             'Color', colors(ci,:), 'LineWidth', 1.8, 'MarkerSize', 6); hold on;
end
hold off; grid on;
xlabel('BSC error rate  p');
ylabel('tr\{LB_\tau\} = \Sigma_p LB(\tau_p)  [s^2]');
title(sprintf('Average pseudo-true bound vs p  (reduced EVA, P=%d, SNR=%d dB, %d trials)', ...
              S.P, snr_fix, S.nTrials));
legend(S.configs(:,1), 'Location','northwest');  set(gca,'FontSize',11);
png2 = fullfile(resdir,'mcrb_lb_vs_p.png');  exportgraphics(gcf, png2, 'Resolution',150);

% ===== Figure 2: tr{LB_tau} vs p at fixed SNR = 0 dB ====================
snr_fix = 30;  sig2_fix = S.alpha_ref2 / 10^(snr_fix/10);
figure('Name','Avg LB vs p (fixed SNR)','Position',[280 180 780 520]);
for ci = 1:nC
    semilogy(S.p_list, trLB_BB(sig2_fix, ci, 1:numel(S.p_list)), styles{ci}, ...
             'Color', colors(ci,:), 'LineWidth', 1.8, 'MarkerSize', 6); hold on;
end
hold off; grid on;
xlabel('BSC error rate  p');
ylabel('tr\{LB_\tau\} = \Sigma_p LB(\tau_p)  [s^2]');
title(sprintf('Average pseudo-true bound vs p  (reduced EVA, P=%d, SNR=%d dB, %d trials)', ...
              S.P, snr_fix, S.nTrials));
legend(S.configs(:,1), 'Location','northwest');  set(gca,'FontSize',11);
png2 = fullfile(resdir,'mcrb_lb_vs_p.png');  exportgraphics(gcf, png2, 'Resolution',150);

return;

% ===== Figures 3 & 4: tr{LB_tau} vs SNR, one curve per p ================
% Fig 3: OFDM ICI (ci = 2);  Fig 4: CP-OTFS ICI (ci = 4).  Six p-curves each.
snr_db = -10:2:30;
sigma2 = S.alpha_ref2 ./ 10.^(snr_db/10);
np     = numel(S.p_list);
pcol   = parula(np);                                   % one colour per p
ici_ci = [2 4];                                        % {OFDM ICI, OTFS ICI}
png34  = cell(1, numel(ici_ci));

for k = 1:numel(ici_ci)
    ci = ici_ci(k);
    figure('Name', sprintf('%s: LB vs SNR per p', S.configs{ci,1}), ...
           'Position', [300+40*k 160 780 520]);
    leg = cell(1, np);
    for ip = 1:np
        semilogy(snr_db, trLB(sigma2, ci, ip), '-o', 'Color', pcol(ip,:), ...
                 'LineWidth', 1.8, 'MarkerSize', 5); hold on;
        leg{ip} = sprintf('p = %.1f', S.p_list(ip));
    end
    hold off; grid on;
    xlabel('SNR = |\alpha_1|^2 / \sigma^2  [dB]');
    ylabel('tr\{LB_\tau\} = \Sigma_p LB(\tau_p)  [s^2]');
    title(sprintf('%s: pseudo-true bound vs SNR per p  (reduced EVA, P=%d, %d trials)', ...
                  S.configs{ci,1}, S.P, S.nTrials));
    legend(leg, 'Location','southwest');  set(gca,'FontSize',11);
    fname = sprintf('mcrb_lb_vs_snr_perp_%s.png', lower(strrep(S.configs{ci,2},' ','_')));
    png34{k} = fullfile(resdir, fname);
    exportgraphics(gcf, png34{k}, 'Resolution',150);
end

% ===== Figures 5 & 6: tr{A^-1 B A^-1} vs SNR, one curve per p ===========
% Same as Figs 3 & 4 but the bound DROPS the pseudo-true bias BB, i.e.
% tr{A^-1 B A^-1} = sigma^2*G + F (no (eta_bar-eta_0)(.)^T term).
png56 = cell(1, numel(ici_ci));
for k = 1:numel(ici_ci)
    ci = ici_ci(k);
    figure('Name', sprintf('%s: A^{-1}BA^{-1} vs SNR per p', S.configs{ci,1}), ...
           'Position', [340+40*k 120 780 520]);
    leg = cell(1, np);
    for ip = 1:np
        semilogy(snr_db, trLB_AB(sigma2, ci, ip), '-o', 'Color', pcol(ip,:), ...
                 'LineWidth', 1.8, 'MarkerSize', 5); hold on;
        leg{ip} = sprintf('p = %.1f', S.p_list(ip));
    end
    hold off; grid on;
    xlabel('SNR = |\alpha_1|^2 / \sigma^2  [dB]');
    ylabel('tr\{A^{-1} B A^{-1}\}_\tau = \Sigma_p (A^{-1}BA^{-1})(\tau_p)  [s^2]');
    title(sprintf('%s: tr\\{A^{-1}BA^{-1}\\} vs SNR per p  (reduced EVA, P=%d, %d trials)', ...
                  S.configs{ci,1}, S.P, S.nTrials));
    legend(leg, 'Location','southwest');  set(gca,'FontSize',11);
    fname = sprintf('mcrb_ABA_vs_snr_perp_%s.png', lower(strrep(S.configs{ci,2},' ','_')));
    png56{k} = fullfile(resdir, fname);
    exportgraphics(gcf, png56{k}, 'Resolution',150);
end

% ===== Figures 7 & 8: tr{BB}_tau vs SNR, one curve per p ===============
% Only the pseudo-true parameter-difference term, BB = (eta_bar-eta_0)(.)^T.
% It is SNR-independent, so each p-curve is flat in SNR (a bias floor).
png78 = cell(1, numel(ici_ci));
for k = 1:numel(ici_ci)
    ci = ici_ci(k);
    figure('Name', sprintf('%s: BB (pseudo-true bias) vs SNR per p', S.configs{ci,1}), ...
           'Position', [380+40*k 80 780 520]);
    leg = cell(1, np);
    for ip = 1:np
        semilogy(snr_db, trLB_BB(sigma2, ci, ip), '-o', 'Color', pcol(ip,:), ...
                 'LineWidth', 1.8, 'MarkerSize', 5); hold on;
        leg{ip} = sprintf('p = %.1f', S.p_list(ip));
    end
    hold off; grid on;
    xlabel('SNR = |\alpha_1|^2 / \sigma^2  [dB]');
    ylabel('tr\{BB\}_\tau = \Sigma_p (\eta_{bar}-\eta_0)^2(\tau_p)  [s^2]');
    title(sprintf('%s: pseudo-true bias term vs SNR per p  (reduced EVA, P=%d, %d trials)', ...
                  S.configs{ci,1}, S.P, S.nTrials));
    legend(leg, 'Location','southwest');  set(gca,'FontSize',11);
    fname = sprintf('mcrb_BB_vs_snr_perp_%s.png', lower(strrep(S.configs{ci,2},' ','_')));
    png78{k} = fullfile(resdir, fname);
    exportgraphics(gcf, png78{k}, 'Resolution',150);
end

% ===== Figures 9 & 10: tr|eta_bar - eta_0|_tau vs p ====================
% Absolute pseudo-true delay bias (un-squared), summed over the tau block:
%   tr|eta_bar - eta_0|_tau = Sigma_p |eta_bar - eta_0|(tau_p)   [s],
% using the averaged pseudo-true eta_0 = eta0_bar.  SNR-independent, so plotted
% vs p.  Fig 9 = OFDM ICI (ci=2), Fig 10 = CP-OTFS ICI (ci=4).
absbiasT = zeros(nC, numel(S.p_list));            % Sigma_p |eta_bar - eta0_bar|(tau_p)
for ci = 1:nC
    for ip = 1:numel(S.p_list)
        absbiasT(ci,ip) = sum(S.eta_bar(it) - S.eta0_bar(it,ip,ci));
    end
end
png910 = cell(1, numel(ici_ci));
for k = 1:numel(ici_ci)
    ci = ici_ci(k);
    figure('Name', sprintf('%s: |eta_bar - eta_0|_tau vs p', S.configs{ci,1}), ...
           'Position', [420+40*k 60 780 520]);
    semilogy(S.p_list, absbiasT(ci,:), '-o', 'Color', colors(min(ci,size(colors,1)),:), ...
             'LineWidth', 1.8, 'MarkerSize', 6);
    grid on;
    xlabel('BSC error rate  p');
    ylabel('tr|\eta_{bar} - \eta_0|_\tau = \Sigma_p |\eta_{bar} - \eta_0|(\tau_p)  [s]');
    title(sprintf('%s: absolute pseudo-true delay bias vs p  (reduced EVA, P=%d, %d trials)', ...
                  S.configs{ci,1}, S.P, S.nTrials));
    set(gca,'FontSize',11);
    fname = sprintf('mcrb_absbias_tau_vs_p_%s.png', lower(strrep(S.configs{ci,2},' ','_')));
    png910{k} = fullfile(resdir, fname);
    exportgraphics(gcf, png910{k}, 'Resolution',150);
end

fprintf(['Saved %s\n       %s\n       %s\n       %s\n       %s\n       %s\n' ...
         '       %s\n       %s\n       %s\n       %s\n'], ...
        png1, png2, png34{1}, png34{2}, png56{1}, png56{2}, ...
        png78{1}, png78{2}, png910{1}, png910{2});
