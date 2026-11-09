% mcrb_avg.m
% --------------------------------------------------------------------------
% Monte-Carlo averaged MCRB and CRB for OFDM, multi-target, ICI channel,
% on the reduced (3-path) EVA channel.  Adapted from run_fim_avg.m.
%
%   CRB  = J^{-1}                  (J  = getFIM_OFDM_multi_ICI at X)
%   MCRB = A^{-1} B A^{-1}         (A,B = getAB_OFDM_multi_ICI at (X, Xhat))
% with Xhat the receiver-side symbols from a BSC of crossover p.  Uses the
% REALIZATION functions (not the expected/closed-form versions); the bounds
% are averaged over nTrials i.i.d. QPSK data draws for each p in p_list.
% At p = 0 we have Xhat = X, so A = -J, B = J and MCRB = CRB (sanity check).
% --------------------------------------------------------------------------

clear; clc;
thisdir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisdir, 'utils'));
resdir = fullfile(thisdir, 'results');
warning('off','MATLAB:nearlySingularMatrix');

% ====== Grid / system parameters (EVA) ==================================
N = 512;  M = 128;  deltaf = 15e3;
T = 1/deltaf;  Ts = T;                 % Tcp = 0
sigma2 = 1;
fc = 4e9;  c = 3e8;  velocity = 120*1e3/3600;
nu_max = velocity * fc / c;

% ====== Reduced EVA channel (3 paths) ===================================
time_delays     = [0, 150, 370] * 1e-9;
relative_powers = [0, -1.4, -0.6];
rp_lin   = 10.^(relative_powers/10);
path_amp = sqrt(rp_lin / sum(rp_lin));
P = numel(time_delays);

rng(42);
theta     = 2*pi*rand(P,1) - pi;
nu_vec    = nu_max * cos(theta);
nu_vec(1) = nu_max;
tau_vec   = time_delays(:);
alpha_vec = path_amp(:);

% ====== Sweep / MC settings =============================================
p_list  = [0, 0.1, 0.5];
nTrials = 50;
modulation = 'QPSK';

% tau-index and nu-index within the stacked 4P parameter vector.
idx_tau = 4*(0:P-1) + 3;
idx_nu  = 4*(0:P-1) + 4;

fprintf('OFDM multi-target ICI, reduced EVA (%d paths), N=%d M=%d, sigma2=%.2g\n', P, N, M, sigma2);
fprintf('nu_max=%.1f Hz, nTrials=%d, p = [%s]\n\n', nu_max, nTrials, num2str(p_list));
for q = 1:P
    fprintf('  path %d: tau=%6.1f ns, nu=%+7.1f Hz, |alpha|^2=%.3f\n', ...
            q, tau_vec(q)*1e9, nu_vec(q), abs(alpha_vec(q))^2);
end
fprintf('\n');

% ====== Monte-Carlo: affine MCRB components and CRB =====================
% At base sigma2 the per-realization bounds are
%   CRB  = sigma2 * inv(Jhat),                 Jhat = sigma2 * FIM(X)
%   MCRB = sigma2 * G + F,  G = A^-1 Jtilde A^-1,  F = A^-1 (b b^T) A^-1,
% all computed once at sigma2 = 1; the SNR sweep is exact (linear in sigma2
% for CRB, affine for MCRB).  We average the sigma2=1 component matrices.
Jinv_acc = zeros(4*P);                          % mean inv(FIM) at sigma2=1 = CRB
G_acc = cell(numel(p_list),1);  F_acc = cell(numel(p_list),1);
for ip = 1:numel(p_list), G_acc{ip} = zeros(4*P); F_acc{ip} = zeros(4*P); end

tic
for t = 1:nTrials
    X = gen_symbols(N, M, modulation);

    J = getFIM_OFDM_multi_ICI(X, deltaf, Ts, 1, alpha_vec, tau_vec, nu_vec);
    Jinv_acc = Jinv_acc + inv(J);

    for ip = 1:numel(p_list)
        Xhat = bsc_qpsk(X, p_list(ip));
        [A, ~, Jt, b] = getAB_OFDM_core(X, Xhat, deltaf, Ts, 1, ...
                                        alpha_vec, tau_vec, nu_vec, true);
        G = A \ Jt / A;                 % sigma2-scaling term
        F = A \ (b*b.') / A;            % bias noise floor (sigma2-independent)
        G_acc{ip} = G_acc{ip} + (G + G.')/2;
        F_acc{ip} = F_acc{ip} + (F + F.')/2;
    end
end
fprintf('MC time: '); toc

CRB_bar = Jinv_acc / nTrials;                   % sigma2 = 1
Gbar = cell(numel(p_list),1);  Fbar = cell(numel(p_list),1);
MCRB_bar = cell(numel(p_list),1);               % sigma2 = 1 (= G + F)
for ip = 1:numel(p_list)
    Gbar{ip} = G_acc{ip} / nTrials;
    Fbar{ip} = F_acc{ip} / nTrials;
    MCRB_bar{ip} = Gbar{ip} + Fbar{ip};
end

% ====== Per-path bounds =================================================
CRB_tau = real(diag(CRB_bar(idx_tau, idx_tau))).';   % s^2
CRB_nu  = real(diag(CRB_bar(idx_nu,  idx_nu ))).';   % Hz^2

MCRB_tau = zeros(numel(p_list), P);
MCRB_nu  = zeros(numel(p_list), P);
for ip = 1:numel(p_list)
    MCRB_tau(ip,:) = real(diag(MCRB_bar{ip}(idx_tau, idx_tau))).';
    MCRB_nu(ip,:)  = real(diag(MCRB_bar{ip}(idx_nu,  idx_nu ))).';
end

% ====== Report ==========================================================
fprintf('\n===== CRB (J^-1), averaged =====\n');
for q = 1:P
    fprintf('  path %d:  CRB(tau)=%.4e s^2   CRB(nu)=%.4e Hz^2\n', q, CRB_tau(q), CRB_nu(q));
end
for ip = 1:numel(p_list)
    fprintf('\n===== MCRB (A^-1 B A^-1), p = %.2f =====\n', p_list(ip));
    for q = 1:P
        fprintf('  path %d:  MCRB(tau)=%.4e s^2  (x%.3f CRB)   MCRB(nu)=%.4e Hz^2  (x%.3f CRB)\n', ...
            q, MCRB_tau(ip,q), MCRB_tau(ip,q)/CRB_tau(q), ...
               MCRB_nu(ip,q),  MCRB_nu(ip,q)/CRB_nu(q));
    end
end

% ====== Save ============================================================
save(fullfile(resdir, 'mcrb_avg_results.mat'), ...
     'p_list','nTrials','N','M','deltaf','Ts','sigma2', ...
     'alpha_vec','tau_vec','nu_vec','P','nu_max','idx_tau','idx_nu', ...
     'CRB_bar','MCRB_bar','Gbar','Fbar','CRB_tau','CRB_nu','MCRB_tau','MCRB_nu');

warning('on','MATLAB:nearlySingularMatrix');

% ------------------------------------------------------------------------
function S = gen_symbols(N, M, mod_type)
    switch upper(mod_type)
        case 'QPSK'
            b = randi([0 3], N, M);
            S = (1/sqrt(2)) * ((2*bitand(b,1) - 1) + 1j*(2*floor(b/2) - 1));
        otherwise
            error('Only QPSK supported here.');
    end
end

function Xh = bsc_qpsk(X, p)
% Gray-coded QPSK: flip each component bit independently with probability p.
    sR = sign(real(X));  sI = sign(imag(X));
    fR = 1 - 2*(rand(size(X)) < p);
    fI = 1 - 2*(rand(size(X)) < p);
    Xh = (sR.*fR + 1j*(sI.*fI)) / sqrt(2);
end
