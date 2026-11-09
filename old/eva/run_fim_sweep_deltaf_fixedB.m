% run_fim_sweep_deltaf_fixedB.m
% --------------------------------------------------------------------------
% Single-waveform FIM/CRB driver that SWEEPS the subcarrier spacing deltaf
% while the total bandwidth B is held FIXED (so the number of subcarriers
% N = B/deltaf VARIES).  Companion to run_fim_sweep_deltaf_fixedN.m, which
% instead fixes N and lets B float.  Single selectable waveform (see `waveform`
% below), single deltaf-independent channel: no EVA modes.
%
% With B fixed, sweeping deltaf trades subcarriers for spacing: N = B/deltaf
% shrinks as deltaf grows, while the occupied bandwidth (and hence the delay
% resolution scale ~1/B) stays constant.  The channel (tau, nu, alpha) is drawn
% ONCE in physical units and reused at every deltaf, so the CRB curves isolate
% the effect of the grid layout at constant bandwidth.
% --------------------------------------------------------------------------

clear; clc;
thisdir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisdir, 'utils'));
resdir = fullfile(thisdir, 'results');
warning('off','MATLAB:nearlySingularMatrix');

% ====== Fixed bandwidth and sweep axis ===================================
B          = 30.72e6;             % total bandwidth [Hz] (FIXED)
                                  % Chosen so the delay bin 1/B ~= 33 ns is finer
                                  % than the 40 ns min path separation, keeping the
                                  % multi-path delay FIM identifiable across N.
M          = 16;                  % OFDM symbols / time slots
modulation = '16-PSK';            % 'QPSK' | '16-QAM' | '16-PSK'
waveform   = 'cp_ofdm';           % 'cp_ofdm' | 'rcp_ofdm' | 'cp_otfs' | 'rcp_otfs'

deltaf_vec = [240e3 480e3 960e3 1.92e6 3.84e6 7.68e6];   % [Hz] -> N = [128 64 32 16 8 4]
N_vec      = round(B ./ deltaf_vec);               % subcarriers per deltaf (VARIES)
nDf        = numel(deltaf_vec);
if any(N_vec < 1)
    error('deltaf too large for B=%.3g Hz: N = B/deltaf < 1.', B);
end

% Map the selected waveform to its (ICI, no-ICI) multi-target FIM builders.
% All four families share the signature f(S, deltaf, Ts, sigma2, alpha, tau, nu)
% with S the data symbols on the N x M grid (TF for OFDM, DD for OTFS).
[fimICI, fimNoICI, wf_label] = select_waveform_fim(waveform);

% ====== System / channel constants (deltaf-independent) ==================
c        = 3e8;
fc       = 4e9;
velocity = 0;                     % no Doppler (set to 120*1e3/3600 for 120 km/h)
% velocity = 120 * 1e3 / 3600;    % 120 km/h
nu_max   = velocity * fc / c;     % peak Doppler [Hz]
sigma2   = 1;

% ----- Channel: drawn once, in physical units (s, Hz), reused for all deltaf
L_paths  = 5;
tau_max  = 300e-9;                % max delay [s]
min_sep  = 40e-9;                 % min delay separation [s] (>= ~1 delay bin, keeps
                                  % paths resolvable / FIM well-conditioned)
rng(20260712);
usable   = tau_max - (L_paths-1)*min_sep;
if usable <= 0
    error('tau_max too small for %d paths with min separation %.1f ns', L_paths, min_sep*1e9);
end
% Stick-breaking draw: random but every consecutive pair is >= min_sep apart.
tau_vec  = sort(usable * rand(L_paths, 1)) + min_sep*(0:L_paths-1).';  % delays [s]
nu_vec   = nu_max * (2*rand(L_paths, 1) - 1);     % Doppler [Hz]
alpha_vec = ones(L_paths, 1);
P        = L_paths;

% nu is a genuine unknown here (channel carries Doppler).
estimate_doppler = any(nu_vec ~= 0);

% ----- Echo the active configuration ------------------------------------
fprintf('%s deltaf sweep (B fixed)\n', wf_label);
fprintf('B=%.3g Hz, M=%d, modulation=%s, P=%d paths\n', B, M, modulation, P);
fprintf('fc=%.3g Hz, v=%.2f m/s, nu_max=%.3g Hz, sigma2=%.2f\n', ...
        fc, velocity, nu_max, sigma2);
fprintf('deltaf sweep : [%s] Hz  (N = B/deltaf = [%s])\n', ...
        strjoin(compose('%.3g', deltaf_vec), ', '), ...
        strjoin(compose('%d', N_vec), ', '));
for p = 1:P
    fprintf('  path %2d:  tau = %7.2f ns,  nu = %+7.2f Hz,  |alpha|^2 = %.3f\n', ...
            p, tau_vec(p)*1e9, nu_vec(p), abs(alpha_vec(p))^2);
end
fprintf('\n');

% ====== Monte Carlo over independent data draws, per deltaf ==============
nTrials = 500;

% Trace-averaged CRBs vs deltaf (mean over trials and paths).
CRBt_OF_I = zeros(nDf, 1);   CRBn_OF_I = zeros(nDf, 1);   % with ICI
CRBt_OF_N = zeros(nDf, 1);   CRBn_OF_N = zeros(nDf, 1);   % no ICI

tic
for id = 1:nDf
    deltaf = deltaf_vec(id);
    N      = N_vec(id);           % subcarriers for this deltaf (B/deltaf)
    Ts     = 1/deltaf;            % Tcp = 0, so Ts = T = 1/deltaf

    crbt_I = zeros(nTrials, P);   crbn_I = zeros(nTrials, P);
    crbt_N = zeros(nTrials, P);   crbn_N = zeros(nTrials, P);
    for t = 1:nTrials
        X = gen_symbols(N, M, modulation);

        J_OF_I = fimICI  (X, deltaf, Ts, sigma2, alpha_vec, tau_vec, nu_vec);
        J_OF_N = fimNoICI(X, deltaf, Ts, sigma2, alpha_vec, tau_vec, nu_vec);

        [crbt_I(t,:), crbn_I(t,:)] = crb_multi(J_OF_I, P, estimate_doppler);
        [crbt_N(t,:), crbn_N(t,:)] = crb_multi(J_OF_N, P, estimate_doppler);
    end

    % Trace-average: sum over paths of the per-path mean CRB.
    CRBt_OF_I(id) = sum(mean(crbt_I, 1));
    CRBn_OF_I(id) = sum(mean(crbn_I, 1));
    CRBt_OF_N(id) = sum(mean(crbt_N, 1));
    CRBn_OF_N(id) = sum(mean(crbn_N, 1));

    fprintf('deltaf=%8.3g Hz (N=%4d, B=%8.3g Hz):  CRB(tau) ICI=%.4e  noICI=%.4e', ...
            deltaf, N, N*deltaf, CRBt_OF_I(id), CRBt_OF_N(id));
    if estimate_doppler
        fprintf('  |  CRB(nu) ICI=%.4e  noICI=%.4e', CRBn_OF_I(id), CRBn_OF_N(id));
    end
    fprintf('\n');
end
fprintf('MC time: '); toc

% ====== Save ============================================================
save(fullfile(resdir, 'fim_sweep_deltaf_fixedB_results.mat'), ...
     'deltaf_vec','B','N_vec','M','modulation','waveform', ...
     'CRBt_OF_I','CRBn_OF_I','CRBt_OF_N','CRBn_OF_N', ...
     'sigma2','fc','velocity','nu_max', ...
     'alpha_vec','tau_vec','nu_vec','P','nTrials','estimate_doppler');

warning('on','MATLAB:nearlySingularMatrix');

% ====== Plot: trace CRLB of the tau parameters vs deltaf =================
%   trace CRLB(tau) = sum_p CRB(tau_p)   [s^2]
% CRBt_OF_I / CRBt_OF_N already hold the per-deltaf trace CRB(tau) for the
% ICI and no-ICI CP-OFDM FIMs.
figure('Name',sprintf('%s trace CRLB(\\tau) vs \\Deltaf (B fixed)', wf_label));
loglog(deltaf_vec, CRBt_OF_I, '-o', 'LineWidth', 2); hold on;
loglog(deltaf_vec, CRBt_OF_N, '-s', 'LineWidth', 2);
grid on;
xlabel('\Delta f  [Hz]');
ylabel('trace CRLB(\tau) = \Sigma_p CRB(\tau_p)   [s^2]');
title(sprintf('%s delay trace-CRLB vs \\Deltaf   (B=%.3g Hz fixed, M=%d, %s, P=%d)', ...
              wf_label, B, M, modulation, P));
legend(sprintf('%s (ICI)', wf_label), sprintf('%s (no-ICI)', wf_label), 'Location', 'best');

% ------------------------------------------------------------------------
function [fimICI, fimNoICI, label] = select_waveform_fim(waveform)
% Map a waveform tag to its (ICI, no-ICI) multi-target FIM builders + label.
    switch lower(waveform)
        case 'cp_ofdm'
            fimICI   = @getFIM_OFDM_multi_ICI;
            fimNoICI = @getFIM_OFDM_multi_noICI;
        case 'rcp_ofdm'
            fimICI   = @getFIM_RCP_OFDM_multi_ICI;
            fimNoICI = @getFIM_RCP_OFDM_multi_noICI;
        case 'cp_otfs'
            fimICI   = @getFIM_CP_OTFS_multi_ICI;
            fimNoICI = @getFIM_CP_OTFS_multi_noICI;
        case 'rcp_otfs'
            fimICI   = @getFIM_RCP_OTFS_multi_ICI;
            fimNoICI = @getFIM_RCP_OTFS_multi_noICI;
        otherwise
            error('Unknown waveform "%s" (use cp_ofdm, rcp_ofdm, cp_otfs, rcp_otfs)', waveform);
    end
    label = upper(strrep(waveform, '_', '-'));
end

function S = gen_symbols(N, M, mod_type)
% i.i.d. zero-mean, unit-power, zero-pseudovariance data symbols on the N x M grid.
    switch upper(mod_type)
        case 'QPSK'
            b = randi([0 3], N, M);
            S = (1/sqrt(2)) * ((2*bitand(b,1) - 1) + 1j*(2*floor(b/2) - 1));
        case '16-QAM'
            levels = [-3 -1 1 3];
            re = levels(randi(4, N, M));
            im = levels(randi(4, N, M));
            S  = (re + 1j*im) / sqrt(10);          % E[|s|^2] = 1
        case '16-PSK'
            k = randi(16, N, M) - 1;
            S = exp(1j*2*pi*k/16);
        otherwise
            error('Unknown modulation "%s" (use QPSK, 16-QAM, 16-PSK)', mod_type);
    end
end

function [crb_tau, crb_nu] = crb_multi(J, P, estimate_doppler)
% Extract per-path CRB(tau_p) (and optionally CRB(nu_p)) from a 4P x 4P
% multi-target FIM via Schur complement w.r.t. the (alphaR, alphaI) block.
    if estimate_doppler
        idx_a = zeros(1, 2*P);
        idx_k = zeros(1, 2*P);
        for p = 1:P
            base = 4*(p-1);
            idx_a(2*p-1:2*p) = [base+1, base+2];
            idx_k(2*p-1:2*p) = [base+3, base+4];
        end
        Jaa = J(idx_a, idx_a);
        Jak = J(idx_a, idx_k);
        Jkk = J(idx_k, idx_k);
        Jschur  = Jkk - Jak' * (Jaa \ Jak);
        C       = pinv_equilibrated((Jschur + Jschur')/2);
        diag_C  = diag(C);
        crb_tau = diag_C(1:2:end).';
        crb_nu  = diag_C(2:2:end).';
    else
        idx_a = zeros(1, 2*P);
        idx_t = zeros(1, P);
        for p = 1:P
            base = 4*(p-1);
            idx_a(2*p-1:2*p) = [base+1, base+2];
            idx_t(p)         = base+3;
        end
        Jaa = J(idx_a, idx_a);
        Jat = J(idx_a, idx_t);
        Jtt = J(idx_t, idx_t);
        Jschur  = Jtt - Jat' * (Jaa \ Jat);
        C       = pinv_equilibrated((Jschur + Jschur')/2);
        crb_tau = diag(C).';
        crb_nu  = nan(1, P);
    end
end

function C = pinv_equilibrated(J)
% pinv of a symmetric FIM after symmetric diagonal equilibration.
    d = sqrt(abs(diag(J)));
    d(d <= 0) = 1;
    Dinv = diag(1 ./ d);
    C = Dinv * pinv(Dinv * J * Dinv) * Dinv;
end
