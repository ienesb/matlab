% plot_mcrb_compare_true_pt.m
% --------------------------------------------------------------------------
% Compare the two paper-scenario MCRBs on the SAME channel / grid / p-sweep:
%   TRUE-param (compute_mcrb_mc_true):  eta_0 = eta_bar, BB = 0,
%       LB = A^-1 B A^-1 = sigma^2 G + F   (A,B at the true parameter)
%   PSEUDO-TRUE (compute_mcrb_mc):       eta_0 = argmin projection,
%       LB = A^-1 B A^-1 + BB             (A,B at eta_0; BB = (eta_bar-eta_0)(.)^T)
% Trace RANGE CRB over the delay-only targets, robust MEDIAN over trials.
% --------------------------------------------------------------------------
clear; clc; close all;
thisdir = fileparts(mfilename('fullpath'));
resdir  = fullfile(thisdir, 'results');
T = load(fullfile(resdir,'compute_mcrb_mc_true_paper_results.mat'));   % true-param
Q = load(fullfile(resdir,'compute_mcrb_mc_paper_results.mat'));        % pseudo-true

u  = (3e8/2)^2;                              % tau [s^2] -> range [m^2]
nC = size(Q.configs,1);  names = Q.configs(:,1);
cols = [0.85 0.33 0.10; 0.00 0.45 0.74; 0.93 0.69 0.13];
assert(numel(T.p_list)==numel(Q.p_list) && max(abs(T.p_list(:)-Q.p_list(:)))<1e-12, ...
       'the two runs used different p grids');
p = Q.p_list;

% Median tau-trace range helpers (per config ci, p-index ip, noise sigma^2 = s).
TRUE_AB = @(s,ci,ip) u*median(s*T.trG_all(:,ip,ci) + T.trF_all(:,ip,ci));   % true total = A^-1 B A^-1
PT_AB   = @(s,ci,ip) u*median(s*Q.trG_all(:,ip,ci) + Q.trF_all(:,ip,ci));   % pt A^-1 B A^-1
PT_BB   = @(ci,ip)   u*median(Q.trBB_all(:,ip,ci));                          % pt bias term
PT_TOT  = @(s,ci,ip) u*median(s*Q.trG_all(:,ip,ci) + Q.trF_all(:,ip,ci) + Q.trBB_all(:,ip,ci)); % pt total

% ===== Console comparison at SNR = 20 dB ================================
s20 = T.alpha_ref2 / 100;                    % sigma^2 at 20 dB
fprintf('\nComparison at SNR = 20 dB  [range m^2, median over trials]\n');
fprintf('%-5s %-5s %12s %12s %12s %12s %9s\n', ...
        'wave','p','true(ABA)','pt(ABA)','pt(BB)','pt(total)','tot/true');
for ci = 1:nC
    for pp = [0.02 0.05 0.10]
        [~,ip] = min(abs(p-pp));
        a=TRUE_AB(s20,ci,ip); b=PT_AB(s20,ci,ip); c=PT_BB(ci,ip); d=PT_TOT(s20,ci,ip);
        fprintf('%-5s %-5.2f %12.4e %12.4e %12.4e %12.4e %9.1f\n', names{ci}, pp, a,b,c,d, d/a);
    end
end

% ===== Figure 1: OFDM, vs p at 20 dB (decomposition + comparison) =======
ci = 1;
figure('Color','w','Position',[200 200 820 540]);
semilogy(p, arrayfun(@(ip)TRUE_AB(s20,ci,ip),1:numel(p)), '-',  'Color',[0 0.45 0.74],'LineWidth',2.0); hold on;
semilogy(p, arrayfun(@(ip)PT_AB(s20,ci,ip),  1:numel(p)), '--', 'Color',[0 0.45 0.74],'LineWidth',1.6);
semilogy(p, arrayfun(@(ip)PT_BB(ci,ip),      1:numel(p)), ':',  'Color',[0.4 0.4 0.4],'LineWidth',1.8);
semilogy(p, arrayfun(@(ip)PT_TOT(s20,ci,ip), 1:numel(p)), '-',  'Color',[0.85 0.33 0.10],'LineWidth',2.4);
hold off; grid on; box on;
xlabel('16-PSK detection-error rate  p');
ylabel('tr\{CRB_R\}_\tau  [m^2]');
title('OFDM @ SNR=20 dB:  true-param  vs  pseudo-true  MCRB');
legend({'true: A^{-1}BA^{-1} (=total, BB=0)','pt: A^{-1}BA^{-1}', ...
        'pt: BB (bias)','pt: total = A^{-1}BA^{-1}+BB'}, 'Location','southeast');
theme(gcf,"light"); set(gca,'FontSize',11);
exportgraphics(gcf, fullfile(resdir,'mcrb_compare_true_pt_vs_p.png'),'Resolution',150,'BackgroundColor','white');

% ===== Figure 2: OFDM, vs SNR at p=0.1 (floor divergence) ==============
[~,ipm] = min(abs(p-0.10));
snr = -10:1:30;  s = T.alpha_ref2 ./ 10.^(snr/10);
figure('Color','w','Position',[240 220 820 540]);
semilogy(snr, arrayfun(@(x)TRUE_AB(x,ci,ipm),s), '--','Color',[0 0.45 0.74],'LineWidth',2.0); hold on;
semilogy(snr, arrayfun(@(x)PT_AB(x,ci,ipm),  s), ':', 'Color',[0 0.45 0.74],'LineWidth',1.6);
semilogy(snr, arrayfun(@(x)PT_TOT(x,ci,ipm), s), '-', 'Color',[0.85 0.33 0.10],'LineWidth',2.4);
hold off; grid on; box on;
xlabel('SNR = |\alpha_1|^2 / \sigma^2  [dB]');
ylabel('tr\{CRB_R\}_\tau  [m^2]');
title(sprintf('OFDM @ p=%.2f:  true-param (A^{-1}BA^{-1})  vs  pseudo-true (total)', p(ipm)));
legend({'true: A^{-1}BA^{-1} (total)','pt: A^{-1}BA^{-1}','pt: total (A^{-1}BA^{-1}+BB)'}, ...
       'Location','southwest');
theme(gcf,"light"); set(gca,'FontSize',11);
exportgraphics(gcf, fullfile(resdir,'mcrb_compare_true_pt_vs_snr.png'),'Resolution',150,'BackgroundColor','white');

fprintf('\nSaved mcrb_compare_true_pt_vs_p.png and mcrb_compare_true_pt_vs_snr.png in %s\n', resdir);
