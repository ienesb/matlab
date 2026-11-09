% pseudotrue_diff.m
% Pseudo-true vs true DELAY difference for OFDM and CP-OTFS on the reduced EVA
% channel (multi-target), as a function of the BSC error rate p at fixed
% SNR = 0 dB.  Four curves: {OFDM, CP-OTFS} x {no-ICI, ICI}.
%
% For each p and data realization the pseudo-true parameter
%   eta_0 = argmin_eta || mu(d, eta_bar) - mu_tilde(dhat, eta) ||   (pseudo_true.tex)
% is computed (pseudoTrue_OFDM for OFDM, pseudoTrue_CP_OTFS for CP-OTFS), and we
% MC-average the trace of the squared delay error
%   tr{(d_tau)^2} = sum_p ( tau0_p - tau_bar_p )^2     [s^2].
% Data: i.i.d. unit-energy QPSK with a BSC of rate p, applied in the TF grid for
% OFDM and the delay-Doppler grid for CP-OTFS.  No A/B matrices are formed.
%
% The pseudo-true is the noise-free LS projection of the true mean, so it does
% NOT depend on sigma^2 / SNR; the SNR = 0 dB label is only for consistency with
% the MCRB plots.  The OFDM traces are reused from pseudotrue_diff_ofdm.m's saved
% results when the settings match (recomputed otherwise).
clear; clc; close all;
thisdir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisdir, 'utils'));
resdir = fullfile(thisdir, 'results');
rng(1);

%% Parameters (must match pseudotrue_diff_ofdm.m for the OFDM cache to be reused)
deltaf = 15e3;  Ts = 1/deltaf;
N = 128;  M = 64;  nTrials = 200;
p_list = 0:0.1:0.5;  snr_db = 0;

%% Reduced EVA channel, multi-target
fc = 2e9;  v = 120/3.6;  c = 3e8;  numax = v*fc/c;
tau   = [0; 370; 150]*1e-9;
nu    = numax*cos([0.5; 2.0; 3.5]);
alpha = 10.^([0; -0.6; -1.4]/20) .* exp(1j*[0.7; 1.1; -0.5]);
P = numel(alpha);

configs = { 'OFDM no-ICI', 'ofdm', false
            'OFDM ICI',    'ofdm', true
            'OTFS no-ICI', 'otfs', false
            'OTFS ICI',    'otfs', true };
nC = size(configs,1);
trace_dtau = nan(nC, numel(p_list));

%% Reuse cached OFDM traces when the settings match
cache = fullfile(resdir, 'pseudotrue_diff_ofdm_results.mat');
have_ofdm = false;
if isfile(cache)
    C = load(cache);
    ok = C.N==N && C.M==M && C.nTrials==nTrials && size(C.trace_dtau,1)==2 ...
         && numel(C.p_list)==numel(p_list) && max(abs(C.p_list(:)-p_list(:)))<1e-12 ...
         && norm(C.tau-tau)<1e-18 && norm(C.nu-nu)<1e-9 && norm(C.alpha-alpha)<1e-12;
    if ok
        trace_dtau(1:2,:) = C.trace_dtau;  have_ofdm = true;
        fprintf('Reusing cached OFDM traces (pseudotrue_diff_ofdm_results.mat)\n');
    end
end

%% Monte-Carlo over the remaining configs
for ci = 1:nC
    if ci <= 2 && have_ofdm, continue; end
    mod = configs{ci,2};  use_ici = configs{ci,3};
    tic;
    for ip = 1:numel(p_list)
        pp = p_list(ip);  acc = 0;
        for t = 1:nTrials
            d  = qpsk(N,M);
            dh = bsc_qpsk(d, pp);
            if strcmp(mod,'otfs')
                [~,~,tau0] = pseudoTrue_CP_OTFS(d, dh, deltaf, Ts, alpha, tau, nu, use_ici);
            else
                [~,~,tau0] = pseudoTrue_OFDM   (d, dh, deltaf, Ts, alpha, tau, nu, use_ici);
            end
            acc = acc + sum((tau0 - tau).^2);
        end
        trace_dtau(ci,ip) = acc/nTrials;
    end
    fprintf('  %-12s done (%.0f s)\n', configs{ci,1}, toc);
end

save(fullfile(resdir,'pseudotrue_diff_results.mat'), ...
     'trace_dtau','p_list','configs','N','M','P','nTrials','alpha','tau','nu','snr_db');

%% Plot: 4 curves (OFDM solid, OTFS dashed)
figure('Position',[240 200 820 540]);
styles = {'-o','-s','--^','--d'};
colors = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.47 0.67 0.19; 0.49 0.18 0.56];
for ci = 1:nC
    semilogy(p_list, trace_dtau(ci,:), styles{ci}, 'Color', colors(ci,:), ...
             'LineWidth', 1.8, 'MarkerSize', 6); hold on;
end
hold off; grid on;
pos = trace_dtau(trace_dtau > 1e-30);    % p=0 is a numerical zero -> off-axis
if ~isempty(pos), ylim([10^floor(log10(min(pos))-0.5) 10^ceil(log10(max(pos))+0.5)]); end
xlabel('BSC error rate  p');
ylabel('tr\{(\Delta\tau)^2\} = \Sigma_p (\tau_{0,p}-\tau_p)^2  [s^2]');
title(sprintf('Pseudo-true vs true delay difference: OFDM vs CP-OTFS  (reduced EVA, P=%d, SNR=%d dB, %d trials)', ...
              P, snr_db, nTrials));
legend(configs(:,1), 'Location','southeast');  set(gca,'FontSize',11);

out_png = fullfile(resdir,'pseudotrue_diff_compare.png');
exportgraphics(gcf, out_png, 'Resolution',150);
fprintf('Saved %s\n', out_png);

%% ---- helpers ----
function X = qpsk(N,M)
    b = randi([0 3], N, M);  X = (1/sqrt(2))*((2*bitand(b,1)-1)+1j*(2*floor(b/2)-1));
end
function Xh = bsc_qpsk(X, p)
    sR = sign(real(X)); sI = sign(imag(X));
    fR = 1 - 2*(rand(size(X)) < p);  fI = 1 - 2*(rand(size(X)) < p);
    Xh = (sR.*fR + 1j*(sI.*fI))/sqrt(2);
end
