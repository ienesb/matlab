% pseudotrue_diff_ofdm.m
% Difference between the PSEUDO-TRUE and TRUE parameters for OFDM on the reduced
% EVA channel (multi-target), versus the BSC error rate p at fixed SNR = 0 dB.
%
% For each p and each i.i.d. QPSK data realization, the pseudo-true parameter
%   eta_0 = argmin_eta || mu(X, eta_bar) - mu_tilde(Xhat, eta) ||   (pseudo_true.tex)
% is computed with pseudoTrue_OFDM (Xhat = BSC decisions of X at rate p).  We
% report the Monte-Carlo-averaged TRACE of the squared delay error
%   tr{(d_tau)^2} = sum_p ( tau0_p - tau_bar_p )^2     [s^2],
% i.e. the parameter-space analogue of tr{MCRB_tau}.  Two curves: OFDM no-ICI
% and OFDM ICI.  No A/B matrices are formed -- only the parameter difference.
%
% NOTE: the pseudo-true parameter is the noise-free least-squares projection of
% the true mean, so it does NOT depend on sigma^2 / SNR.  The SNR = 0 dB setting
% is kept only for consistency with the MCRB plots; it does not affect the curve.
clear; clc; close all;
thisdir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisdir, 'utils'));
resdir = fullfile(thisdir, 'results');
rng(1);

%% Parameters
deltaf = 15e3;  Ts = 1/deltaf;
N = 128;  M = 64;                        % resolvable grid (Jacobi-scaled cond ~ 70)
nTrials = 200;                           % MC trials per p
p_list = 0:0.1:0.5;                      % BSC error-rate sweep
snr_db = 0;                              % fixed (does not enter the pseudo-true)

%% Reduced EVA channel, multi-target (eva_channel.tex, Sec. "Reduced EVA Channel")
fc = 2e9;  v = 120/3.6;  c = 3e8;  numax = v*fc/c;     % max Doppler ~222.2 Hz
tau  = [0;   370;  150]*1e-9;            % TRUE delays   [s]
nu   = numax*cos([0.5; 2.0; 3.5]);      % TRUE Dopplers [Hz]
alpha = 10.^([0; -0.6; -1.4]/20) .* exp(1j*[0.7; 1.1; -0.5]);   % TRUE gains
P = numel(alpha);

configs = { 'OFDM no-ICI', false
            'OFDM ICI',    true };
nC = size(configs,1);

%% Monte-Carlo: trace of the squared delay difference vs p
trace_dtau = zeros(nC, numel(p_list));
for ci = 1:nC
    use_ici = configs{ci,2};
    tic;
    for ip = 1:numel(p_list)
        pp = p_list(ip);  acc = 0;
        for t = 1:nTrials
            X  = qpsk(N,M);
            Xh = bsc_qpsk(X, pp);
            [~, ~, tau0] = pseudoTrue_OFDM(X, Xh, deltaf, Ts, alpha, tau, nu, use_ici);
            acc = acc + sum((tau0 - tau).^2);
        end
        trace_dtau(ci,ip) = acc/nTrials;
    end
    fprintf('  %-12s done (%.0f s)\n', configs{ci,1}, toc);
end

save(fullfile(resdir,'pseudotrue_diff_ofdm_results.mat'), ...
     'trace_dtau','p_list','configs','N','M','P','nTrials','alpha','tau','nu','snr_db');

%% Plot: trace of tau difference vs p
figure('Name','Pseudo-true vs true delay difference (OFDM, reduced EVA)', ...
       'Position',[240 220 760 520]);
styles = {'-o','-s'};
colors = [0.00 0.45 0.74; 0.85 0.33 0.10];
for ci = 1:nC
    semilogy(p_list, trace_dtau(ci,:), styles{ci}, 'Color', colors(ci,:), ...
             'LineWidth', 1.9, 'MarkerSize', 7); hold on;
end
hold off; grid on;
% p=0 gives a numerical zero (eta_0 = eta_bar); clip the axis to the meaningful range.
pos = trace_dtau(trace_dtau > 1e-30);
if ~isempty(pos), ylim([10^floor(log10(min(pos))-0.5) 10^ceil(log10(max(pos))+0.5)]); end
xlabel('BSC error rate  p');
ylabel('tr\{(\Delta\tau)^2\} = \Sigma_p (\tau_{0,p}-\tau_p)^2  [s^2]');
title(sprintf('Pseudo-true vs true delay difference  (reduced EVA, OFDM, P=%d, SNR=%d dB, %d trials)', ...
              P, snr_db, nTrials));
legend(configs(:,1), 'Location','southeast');
set(gca,'FontSize',11);

out_png = fullfile(resdir,'pseudotrue_diff_tau.png');
exportgraphics(gcf, out_png, 'Resolution',150);
fprintf('Saved %s\n', out_png);

%% ---- helpers ----
function X = qpsk(N,M)
    b = randi([0 3], N, M);
    X = (1/sqrt(2)) * ((2*bitand(b,1)-1) + 1j*(2*floor(b/2)-1));
end
function Xh = bsc_qpsk(X, p)
    sR = sign(real(X)); sI = sign(imag(X));
    fR = 1 - 2*(rand(size(X)) < p);
    fI = 1 - 2*(rand(size(X)) < p);
    Xh = (sR.*fR + 1j*(sI.*fI))/sqrt(2);
end
