% plot_mcrb_mmle.m
% --------------------------------------------------------------------------
% Combined view of the pseudo-true lower bound (MCRB) and the empirical
% misspecified-MLE (MMLE) delay error on the reduced EVA channel (multi-target),
% trace over the delay (tau) block, for {OFDM, CP-OTFS} x {no-ICI, ICI}.
% Overlays plot_mcrb.m and plot_mmle.m on the same axes:
%
%   MCRB (dashed):  tr{LB_tau}  = sum_p LB(tau_p),  LB(sigma^2)=sigma^2 G + F + BB
%                   from compute_mcrb_mc_results.mat (averaged affine pieces).
%   MMLE (solid):   tr{MSE_tau} = sum_p MSE_bar(tau_p,tau_p)
%                   from compute_mmle_mc_results.mat (empirical, about eta_bar).
%   SNR := |alpha_1|^2 / sigma^2.  Both are measured about eta_bar, so the MMLE
%   should sit at or above its MCRB.  Color = config, dashed = MCRB, solid = MMLE.
%
% Two figures, eight curves each (four configs x {MCRB, MMLE}):
%   Fig 1: fixed p = 0.1,    x-axis = SNR [dB]
%   Fig 2: fixed SNR = 0 dB, x-axis = p
% --------------------------------------------------------------------------
clear; clc; close all;
thisdir = fileparts(mfilename('fullpath'));
resdir  = fullfile(thisdir, 'results');

% ----- MCRB (pseudo-true lower bound) ------------------------------------
C = load(fullfile(resdir, 'compute_mcrb_mc_results.mat'), ...
         'G_bar','F_bar','BB_bar','configs','p_list','alpha_ref2','idx_tau','P','nTrials');
% ----- MMLE (empirical MSE) ----------------------------------------------
U = load(fullfile(resdir, 'compute_mmle_mc_results.mat'), ...
         'MSE_bar','configs','p_list','snr_db_list','idx_tau','P','nTrials');

assert(isequal(C.configs(:,1), U.configs(:,1)), 'config mismatch between MCRB and MMLE results');
nC = size(C.configs,1);
itC = C.idx_tau;
% itC = 3;

itU = U.idx_tau;
% itU = 3;

% MCRB tau-traces of the three averaged pieces, per (config, p).
npC = numel(C.p_list);
gT = zeros(nC, npC);  fT = zeros(nC, npC);  bT = zeros(nC, npC);
for ci = 1:nC
    for ip = 1:npC
        gT(ci,ip) = sum(real(diag(C.G_bar (itC,itC,ip,ci))));
        fT(ci,ip) = sum(real(diag(C.F_bar (itC,itC,ip,ci))));
        bT(ci,ip) = sum(real(diag(C.BB_bar(itC,itC,ip,ci))));
    end
end
trLB    = @(sig2, ci, ip) sig2.*gT(ci,ip) + fT(ci,ip) + bT(ci,ip);   % tr{LB_tau}
trLB_AB = @(sig2, ci, ip) sig2.*gT(ci,ip) + fT(ci,ip);               % tr{A^-1 B A^-1} (no BB bias)

% MMLE tau-trace of the empirical MSE, per (config, p, SNR).
npU = numel(U.p_list);  nsU = numel(U.snr_db_list);
mseT = zeros(nC, npU, nsU);
for ci = 1:nC
    for ip = 1:npU
        for is = 1:nsU
            mseT(ci,ip,is) = sum(real(diag(U.MSE_bar(itU,itU,ip,is,ci))));
        end
    end
end

colors  = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.47 0.67 0.19; 0.49 0.18 0.56];
markers = {'o','s','^','d'};

% ===== MCRB (dashed) vs MMLE (solid) vs SNR, fixed p = 0.1 =============
% One figure per ICI config: Fig 1 = OFDM ICI (ci=2), Fig 2 = CP-OTFS ICI (ci=4).
p_fix = 0.1;
ipC = find(abs(C.p_list - p_fix) < 1e-9, 1);  assert(~isempty(ipC), 'p=%.2f not in MCRB p_list', p_fix);
ipU = find(abs(U.p_list - p_fix) < 1e-9, 1);  assert(~isempty(ipU), 'p=%.2f not in MMLE p_list', p_fix);

snr_fine = linspace(min(U.snr_db_list), max(U.snr_db_list), 200);   % MCRB on a fine grid
sig2_fine = C.alpha_ref2 ./ 10.^(snr_fine/10);

ici_ci = [2 4];                       % {OFDM ICI, CP-OTFS ICI}
png    = cell(1, numel(ici_ci));
for k = 1:numel(ici_ci)
    ci = ici_ci(k);
    figure('Name', sprintf('MCRB vs MMLE: %s, delay error vs SNR (p=%.1f)', C.configs{ci,1}, p_fix), ...
           'Position', [200+40*k 220 820 560]);
    hb = semilogy(snr_fine, trLB(sig2_fine, ci, ipC), '--', ...
                  'Color', colors(ci,:), 'LineWidth', 1.6); hold on;
    hm = semilogy(U.snr_db_list, squeeze(mseT(ci,ipU,:)), ['-' markers{ci}], ...
                  'Color', colors(ci,:), 'LineWidth', 1.9, 'MarkerSize', 6);
    hold off; grid on;
    xlabel('SNR = |\alpha_1|^2 / \sigma^2  [dB]');
    ylabel('tr\{\cdot_\tau\} = \Sigma_p (\tau_p)  [s^2]');
    title(sprintf('%s: MCRB (dashed) vs MMLE (solid) vs SNR  (reduced EVA, P=%d, p=%.1f, %d trials)', ...
                  C.configs{ci,1}, U.P, p_fix, U.nTrials));
    legend([hb hm], {[C.configs{ci,1} ' (MCRB)'], [C.configs{ci,1} ' (MMLE)']}, ...
           'Location','southwest');  set(gca,'FontSize',11);
    fname = sprintf('mcrb_mmle_vs_snr_%s.png', lower(strrep(C.configs{ci,2},' ','_')));
    png{k} = fullfile(resdir, fname);  exportgraphics(gcf, png{k}, 'Resolution',150);
end

return

% ===== A^-1 B A^-1 (dashed) vs MMLE (solid) vs SNR, fixed p = 0.1 ======
% Same two ICI configs, but the bound DROPS the pseudo-true bias BB, i.e.
% tr{A^-1 B A^-1} = sigma^2*G + F.  Fig 3 = OFDM ICI, Fig 4 = CP-OTFS ICI.
png_ab = cell(1, numel(ici_ci));
for k = 1:numel(ici_ci)
    ci = ici_ci(k);
    figure('Name', sprintf('A^{-1}BA^{-1} vs MMLE: %s, delay error vs SNR (p=%.1f)', C.configs{ci,1}, p_fix), ...
           'Position', [240+40*k 180 820 560]);
    hb = semilogy(snr_fine, trLB_AB(sig2_fine, ci, ipC), '--', ...
                  'Color', colors(ci,:), 'LineWidth', 1.6); hold on;
    hm = semilogy(U.snr_db_list, squeeze(mseT(ci,ipU,:)), ['-' markers{ci}], ...
                  'Color', colors(ci,:), 'LineWidth', 1.9, 'MarkerSize', 6);
    hold off; grid on;
    xlabel('SNR = |\alpha_1|^2 / \sigma^2  [dB]');
    ylabel('tr\{\cdot_\tau\} = \Sigma_p (\tau_p)  [s^2]');
    title(sprintf('%s: A^{-1}BA^{-1} (dashed) vs MMLE (solid) vs SNR  (reduced EVA, P=%d, p=%.1f, %d trials)', ...
                  C.configs{ci,1}, U.P, p_fix, U.nTrials));
    legend([hb hm], {[C.configs{ci,1} ' (A^{-1}BA^{-1})'], [C.configs{ci,1} ' (MMLE)']}, ...
           'Location','southwest');  set(gca,'FontSize',11);
    fname = sprintf('mcrb_mmle_ABA_vs_snr_%s.png', lower(strrep(C.configs{ci,2},' ','_')));
    png_ab{k} = fullfile(resdir, fname);  exportgraphics(gcf, png_ab{k}, 'Resolution',150);
end

fprintf('Saved %s\n       %s\n       %s\n       %s\n', png{1}, png{2}, png_ab{1}, png_ab{2});
