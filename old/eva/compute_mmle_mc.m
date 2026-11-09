% compute_mmle_mc.m
% Monte-Carlo MISSPECIFIED MLE (MMLE) on the reduced EVA channel (multi-target),
% for the four configurations {OFDM, CP-OTFS} x {no-ICI, ICI}, swept over the BSC
% error rate p and the SNR.  This is the estimator counterpart of
% compute_mcrb_mc.m: where that script forms the MCRB component matrices, this one
% actually RUNS the misspecified ML estimator and records its estimates -- NO
% MCRB / A,B,Jtilde,b is computed here.
%
% MMLE (ofdm_mmle.tex / cp_otfs_mmle.tex): the receiver fits the assumed mean
% mu_tilde(eta) built from the DETECTED symbols Xhat (not the transmitted X),
% which makes the model misspecified.  Concentrating the gain out reduces the
% estimator to the delay-Doppler matched-filter (periodogram) search
%   (tau_hat,nu_hat) = argmax |c_tilde^H y|^2 / ||c_tilde||^2,
%   alpha_hat        = c_tilde^H y / ||c_tilde||^2,        c_tilde from Xhat,
% extended to P targets by successive (CLEAN) interference cancellation.  See
% mmle_OFDM.m (the periodogram + refine + CLEAN core) and mmle_CP_OTFS.m (the
% (I)SFFT-equivalence wrapper that maps the DD grids to TF and delegates).
%
% For every realization:
%   1. transmitted symbols  X (TF, OFDM) or x (DD, CP-OTFS);
%      detected symbols     Xhat / xhat  via a BSC of rate p.
%   2. true noise-free mean mu(X, eta_bar) from the TRANSMITTED symbols and the
%      true parameter eta_bar, plus CN(0,sigma^2) noise -> observation y.
%   3. eta_hat = MMLE(y, Xhat)  (the misspecified estimator on the detected ref).
%   4. associate the CLEAN-ordered estimates to the true targets (nearest in the
%      (tau,nu) plane) and accumulate the error about eta_bar.
%
% At p = 0 (Xhat = X) the MMLE is the ordinary (genie) MLE and converges to the
% true eta_bar; with detection errors it converges to the pseudo-true eta_0.
% The empirical MSE about eta_bar recorded here is directly comparable to the
% pseudo-true lower bound LB(sigma^2) of compute_mcrb_mc.m / plot_mcrb.m.
clear; clc;
thisdir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisdir, 'utils'));
resdir = fullfile(thisdir, 'results');
rng(1);   % NOTE: the trial loop runs under parfor, so workers draw from their own
          % RNG streams -- results are statistically equivalent but not bit-identical
          % across runs / worker counts. Start a pool first to control worker count.

%% Parameters (same reduced EVA channel as compute_mcrb_mc.m)
deltaf = 15e3;  Ts = 1/deltaf;          % subcarrier spacing / symbol duration [s]
N = 128;  M = 64;                        % grid
nTrials = 1000;                            % MC trials per (config, p, SNR)
p_list      = 0:0.1:0.5;                 % BSC crossover sweep
snr_db_list = 0:5:20;                    % SNR sweep [dB], SNR = |alpha_1|^2/sigma^2
Zpad = 4;                                % FFT zero-pad factor for the MMLE search

%% Reduced EVA channel, multi-target (eva_channel.tex, "Reduced EVA Channel")
fc = 2e9;  v = 120/3.6;  c = 3e8;  numax = v*fc/c;     % max Doppler ~222.2 Hz
tau   = [0; 370; 150]*1e-9;             % TRUE delays   [s]
nu    = numax*cos([0.5; 2.0; 3.5]);     % TRUE Dopplers [Hz]
alpha = 10.^([0; -0.6; -1.4]/20) .* exp(1j*[0.7; 1.1; -0.5]);   % TRUE gains
P  = numel(alpha);
nE = 4*P;
eta_bar = packEta(alpha, tau, nu);      % ground-truth parameter (4P x 1)
idx_tau = 4*(0:P-1) + 3;
idx_nu  = 4*(0:P-1) + 4;
alpha_ref2 = abs(alpha(1))^2;           % SNR reference power
sigma2_list = alpha_ref2 ./ 10.^(snr_db_list/10);

s_tau = 1/(N*deltaf);  s_nu = 1/(M*Ts); % (tau,nu) bin scales for target matching

configs = { 'OFDM no-ICI', 'ofdm', false
            'OFDM ICI',    'ofdm', true
            'OTFS no-ICI', 'otfs', false
            'OTFS ICI',    'otfs', true };
nC = size(configs,1);
np = numel(p_list);
ns = numel(snr_db_list);

%% Storage: every MC estimate + averages
eta_hat_all = zeros(nE, nTrials, np, ns, nC);   % aligned to eta_bar target order
MSE_bar     = zeros(nE, nE, np, ns, nC);        % empirical MSE about eta_bar
eta_hat_bar = zeros(nE, np, ns, nC);            % mean estimate

fprintf('MC misspecified MLE (MMLE) on the detected-symbol reference\n');
fprintf('reduced EVA, N=%d M=%d P=%d, %d trials, p in [%s], SNR in [%s] dB\n', ...
        N, M, P, nTrials, num2str(p_list), num2str(snr_db_list));

nWork = nC*ns*np*nTrials;   doneWork = 0;   tAll = tic;   % overall progress / ETA

for ci = 1:nC
    name = configs{ci,1};  modn = configs{ci,2};  use_ici = configs{ci,3};
    tic;
    for is = 1:ns
        sigma2 = sigma2_list(is);
        for ip = 1:np
            pp = p_list(ip);

            % Per-iteration sliced temporary (indexed only by the parfor variable t).
            eta_t = zeros(nE, nTrials);     % estimates, aligned to eta_bar target order

            parfor t = 1:nTrials
                % ---- data, true mean, noisy observation ----
                if strcmp(modn,'ofdm')
                    X  = qpsk(N,M);
                    Xh = bsc_qpsk(X, pp);
                    mu = ofdm_mean(X, deltaf, Ts, alpha, tau, nu, use_ici);  % TF
                    Y  = mu + noise(N, M, sigma2);
                    eta_hat = mmle_OFDM(Y, Xh, deltaf, Ts, P, use_ici, Zpad);
                else
                    x  = qpsk(N,M);                 % DD symbols
                    xh = bsc_qpsk(x, pp);
                    Xtf = isfft_u(x);
                    muTF = ofdm_mean(Xtf, deltaf, Ts, alpha, tau, nu, use_ici);
                    y   = sfft_u(muTF) + noise(N, M, sigma2);   % DD observation
                    eta_hat = mmle_CP_OTFS(y, xh, deltaf, Ts, P, use_ici, Zpad);
                end

                % ---- align CLEAN-ordered estimates to the true targets ----
                [ah, th, nh] = unpackEta(eta_hat);
                perm = match_targets(th, nh, tau, nu, s_tau, s_nu);
                eta_t(:,t) = packEta(ah(perm), th(perm), nh(perm));
            end

            eta_hat_all(:,:,ip,is,ci) = eta_t;
            E = eta_t - eta_bar;                            % nE x nTrials errors
            MSE_bar(:,:,ip,is,ci) = (E*E.')/nTrials;        % empirical MSE about eta_bar
            eta_hat_bar(:,ip,is,ci) = mean(eta_t, 2);
            doneWork = doneWork + nTrials;
            print_progress(doneWork, nWork, tAll);
        end
    end
    is0 = 1;  ip0 = 1;  ipL = np;
    rt0 = sqrt(sum(real(diag(MSE_bar(idx_tau,idx_tau,ip0,is0,ci)))));
    rtL = sqrt(sum(real(diag(MSE_bar(idx_tau,idx_tau,ipL,is0,ci)))));
    fprintf('  %-12s done (%.0f s) : RMSE_tau @SNR=%gdB  p=%.1f -> %.3e ,  p=%.1f -> %.3e  [s]\n', ...
            name, toc, snr_db_list(is0), p_list(ip0), rt0, p_list(ipL), rtL);
end

save(fullfile(resdir,'compute_mmle_mc_results.mat'), ...
     'eta_hat_all','MSE_bar','eta_hat_bar','eta_bar', ...
     'configs','p_list','snr_db_list','sigma2_list','alpha_ref2', ...
     'N','M','P','nTrials','alpha','tau','nu','deltaf','Ts', ...
     'idx_tau','idx_nu', '-v7.3');
fprintf('\nSaved per-iteration estimates + empirical MSE_bar to\n  results/compute_mmle_mc_results.mat\n');

%% ---- helpers ----
function print_progress(done, total, tAll)
% One-line progress with percent-done, elapsed and estimated time remaining.
    el  = toc(tAll);
    eta = el/done*(total-done);
    fprintf('    progress %5.1f%%  (%d/%d)  elapsed %s  ETA %s\n', ...
            100*done/total, done, total, fmt_hms(el), fmt_hms(eta));
end

function s = fmt_hms(sec)
    sec = max(0, round(sec));
    s = sprintf('%02d:%02d:%02d', floor(sec/3600), floor(mod(sec,3600)/60), mod(sec,60));
end

function Z = noise(N, M, sigma2)
    Z = sqrt(sigma2/2) * (randn(N,M) + 1j*randn(N,M));
end

function perm = match_targets(tau_h, nu_h, tau_b, nu_b, s_tau, s_nu)
% Greedy nearest assignment (in (tau,nu) bin units): for each true target j (in
% descending-strength order) pick the closest still-unused estimate.  perm(j) is
% the estimate index assigned to true target j, so eta_hat(perm) aligns with the
% eta_bar target ordering for the MSE.
    P = numel(tau_b);
    d = ((tau_h(:).' - tau_b(:))/s_tau).^2 + ((nu_h(:).' - nu_b(:))/s_nu).^2;  % P x P
    perm = zeros(P,1);  usedE = false(1,P);
    for j = 1:P
        dj = d(j,:);  dj(usedE) = inf;
        [~, i] = min(dj);  perm(j) = i;  usedE(i) = true;
    end
end

function mu = ofdm_mean(S, deltaf, Ts, alpha, tau, nu, use_ici)
% Noise-free OFDM mean mu(S, eta) (ofdm_mcrb.tex eqs (3)-(4)); S = transmitted.
    [N, M] = size(S);  P = numel(alpha);
    ns = (0:N-1).';  ms = (0:M-1);  mu = zeros(N, M);
    if use_ici
        T = 1/deltaf;  [NP, NN] = ndgrid(0:N-1, 0:N-1);  diff_n = NP - NN;
    end
    for p = 1:P
        phase = exp(-1j*2*pi*deltaf*ns*tau(p)) * exp(1j*2*pi*Ts*ms*nu(p));
        A = S .* phase;
        if use_ici
            mu = mu + alpha(p) * (get_ici_coeffs(diff_n, nu(p), T).' * A);
        else
            mu = mu + alpha(p) * A;
        end
    end
end

function eta = packEta(alpha, tau, nu)
    eta = reshape([real(alpha(:)).'; imag(alpha(:)).'; tau(:).'; nu(:).'], [], 1);
end

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
