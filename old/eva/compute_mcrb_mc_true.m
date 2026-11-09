% compute_mcrb_mc_true.m
% Same as compute_mcrb_mc.m, EXCEPT the pseudo-true parameter is assumed to
% coincide with the TRUE parameter:  eta_0 := eta_bar.  A and B are therefore
% evaluated AT the true parameter (getAB_OFDM_core with eta = eta_bar), not at
% the BSC-induced pseudo-true point.
%
% Everything else is unchanged so the existing plot scripts (which read the
% pseudo-true variable names eta0_all, eta0_bar, BB_all, BB_bar, ...) keep
% working.  Because eta_0 = eta_bar:
%     bias = eta_bar - eta_0 = 0  ==>  BB = 0  (kept, just identically zero),
% while A, B, Jtilde and b still feel the data mismatch through
%     Delta = sum_p alpha_p * Phi_{E,p},  E = X - Xhat
% (X transmitted, Xhat the BSC-corrupted symbols at rate p).  The lower bound
% reduces to LB(sigma^2) = sigma^2*G + F  (no pseudo-true bias floor).
%
% For every realization (X, Xhat):
%   1. eta_0 = eta_bar               (true parameter -- the assumption here).
%   2. A(eta_bar), B(eta_bar), Jtilde(eta_bar), b(eta_bar)  at sigma^2 = 1
%      (getAB_OFDM_core evaluated AT the true parameter; eta_bar = eta).
%   3. lower bound (MSE matrix):
%        LB(sigma^2) = A^{-1} B A^{-1} + (eta_bar - eta_0)(eta_bar - eta_0)^T
%                    = sigma^2 * G + F + BB,
%        G  = A^{-1} Jtilde A^{-1}         (CRB-like, scales with noise)
%        F  = A^{-1} (b b^T) A^{-1}        (score-bias floor, sigma^2-independent)
%        BB = (eta_bar - eta_0)(eta_bar - eta_0)^T = 0   (eta_0 = eta_bar)
%      all evaluated at sigma^2 = 1.  plot_mcrb.m sweeps the SNR via sigma^2.
%
% A is ill-conditioned only through unit disparity (tau~1e-7, nu~1e2, alpha~1),
% so the A^{-1}(.)A^{-1} products use a Jacobi (magnitude) rescale (real result).
% At p = 0 (Xhat = X): A = -J, b = 0, BB = 0, F = 0, G = J^{-1},
% so LB = sigma^2 J^{-1} (the CRB).
clear; clc;
thisdir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisdir, 'utils'));
resdir = fullfile(thisdir, 'results');
warning('off','MATLAB:nearlySingularMatrix');
rng(1);   % NOTE: the trial loop runs under parfor, so workers draw from their own
          % RNG streams -- results are statistically equivalent but not bit-identical
          % across runs / worker counts. Start a pool first to control worker count.

%% ===== Scenario selector =================================================
% Toggle between the original EVA setup and the paper Fig.1 reproduction.
%   'eva'   : original 64x32 full/reduced-EVA, QPSK, physical OFDM/OTFS(isfft_u).
%   'paper' : Fig.1 scenario matching run_fim_avg -- 128x1 delay-only 40-target
%             channel, 16-PSK, precoded OFDM/SC/OTFS (Q^H = F_N U) fed to the
%             same OFDM core.  Delay resolves over the full N=128 steering, so
%             all 40 targets are resolvable and OTFS keeps its true F_32^H(x)I_4
%             structure (the physical isfft_u OTFS would instead be SC here).
% Everything scenario-specific is set inside the switch; the rest is shared.
scenario = 'paperv2';    % 'eva' | 'paper' | 'paperv2'   <-- switch here

% Common random numbers across waveforms (paired comparison).  With CRN every
% config sees the SAME data symbols and the SAME bit flips at a given (p,trial),
% so the OFDM/SC/OTFS difference is paired and the MC noise largely cancels --
% essential to resolve the few-percent gaps in paperv2 (independent draws need
% ~100x more trials for the same confidence).  It is also the faithful
% "equal-BER, different precoder" comparison.  false = independent draws.
use_crn  = true;
crn_seed = 20260707;

switch lower(scenario)
    case 'eva'
        %% Parameters
        deltaf = 15e3;  Ts = 1/deltaf;   % subcarrier spacing / symbol duration [s]
        N = 64;   M = 32;                % grid (reduced from 128x64 to leave the
                                         % large-NM regime; coarser resolution ->
                                         % A more singular; watch kappa_max).
        % N = 128;  M = 64;              % original values
        nTrials = 10000;                 % MC trials per (config, p)
        p_list  = 0:0.02:0.5;            % BSC crossover sweep
        use_full_eva = true;             % false = reduced EVA (3 taps), true = 9 paths
        store_all = true;                % keep per-iteration arrays (P=9 is small)
        kappa_thresh = inf;              % no Omega_N filtering (preserve original)

        %% EVA channel, multi-target (eva_channel.tex). Doppler angles (Jakes
        % theta_l) and gain phases are deterministic literals (not EVA-specified).
        fc = 2e9;  v = 120/3.6;  c = 3e8;  numax = v*fc/c;   % max Doppler ~222.2 Hz
        if use_full_eva
            % Full EVA, 9 paths (Table "EVA channel power delay profile")
            tau_ns = [0;   30;  150;  310;  370;  710; 1090; 1730; 2510];   % [ns]
            pow_db = [0; -1.5; -1.4; -3.6; -0.6; -9.1; -7.0; -12.0; -16.9]; % [dB]
            theta  = [0.5; 2.0; 3.5; 0.9; 1.6; 4.2; 2.7; 5.1; 3.0];         % angles
            phase  = [0.7; 1.1; -0.5; 0.3; -1.2; 2.0; -0.8; 1.5; 0.2];      % phases
        else
            % Reduced EVA, 3 strongest taps
            tau_ns = [0;  370;  150];
            pow_db = [0; -0.6; -1.4];
            theta  = [0.5; 2.0; 3.5];
            phase  = [0.7; 1.1; -0.5];
        end
        tau   = tau_ns*1e-9;                    % TRUE delays   [s]
        nu    = numax*cos(theta);               % TRUE Dopplers [Hz]
        alpha = 10.^(pow_db/20) .* exp(1j*phase);   % TRUE gains
        chan_desc = 'full EVA';  if ~use_full_eva, chan_desc = 'reduced EVA'; end
        mod_type  = 'qpsk';
        % {name, key, use_ici, transform d -> TF symbols}
        configs = { 'OFDM no-ICI', 'ofdm', false, @(d) d
                    'OFDM ICI',    'ofdm', true,  @(d) d
                    'OTFS no-ICI', 'otfs', false, @(d) isfft_u(d)
                    'OTFS ICI',    'otfs', true,  @(d) isfft_u(d) };
        outfile = 'compute_mcrb_mc_true_results.mat';

    case 'paper'
        %% Parameters (paper Fig.1: B=160 MHz, N=128 single block, delay-only)
        B = 160e6;  N = 128;  M = 1;  deltaf = B/N;  Ts = 1/deltaf;
        nTrials = 10000;                 % raised for smooth averages (reduction path)
        p_list  = 0:0.002:0.1;           % fine 16-PSK BER sweep over the practical low-BER band (51 pts)
        store_all = false;               % use the O(nE^2) reduction path (below); the
                                         % per-iter _all arrays would be tens of GB here.
        kappa_thresh = 1e4;              % Omega_N restriction: average the mismatch
                                         % bound only over well-conditioned trials.
                                         % F ~ A^-2 is heavy-tailed, so rare near-
                                         % singular draws otherwise dominate the mean
                                         % and make the vs-p curve jump.  Raise/lower
                                         % to trade smoothness vs retained fraction.
        pilot_frac = 0.20;               % pilot-only OFDM CRB: fraction of the N
                                         % subcarriers used as known comb pilots
                                         % (>=~15% needed to resolve L=10 targets).

        %% Paper channel: L delay-only targets spread over [0,127]/B.
        % L=40 (paper density) makes the alpha+tau FIM near-rank-critical
        % (3L=120 params in N=128 obs), so the mismatch bound F ~ A^-2 explodes
        % at p>0.  L=10 keeps the FIM comfortably over-determined (30 in 128) so
        % the MCRB vs p is well-posed and smooth.  The p=0 CRB / Fig.1(b) at the
        % full 40-target density is produced separately by run_fim_avg.
        L_paths = 10;  bin_t = 1/(N*deltaf);            % 1/B  (<-- switch L here)
        tau_max = (N-1)*bin_t;  min_sep = 1.0;
        rng(20260514);
        usable = tau_max - (L_paths-1)*min_sep*bin_t;
        u      = sort(usable*rand(L_paths,1));
        tau    = u + min_sep*bin_t*(0:L_paths-1).';     % TRUE delays [s]
        nu     = zeros(L_paths,1);                      % delay-only
        alpha  = ones(L_paths,1);                       % unit gains
        chan_desc = 'paper 40-target delay-only';
        mod_type  = '16psk';

        %% Precoders Q^H = F_N U (paper eq.8): OFDM U=F_N^H, SC U=I,
        %% OTFS U=F_N1^H (x) I_N2.  Only meaningful through use_ici=false (nu=0).
        N1p = 32;  N2p = 4;
        FN  = fft(eye(N))  / sqrt(N);
        FN1 = fft(eye(N1p)) / sqrt(N1p);
        QH_OFDM = eye(N);
        QH_SC   = FN;
        QH_OTFS = FN * kron(FN1', eye(N2p));
        configs = { 'OFDM', 'ofdm', false, @(d) QH_OFDM*d
                    'SC',   'sc',   false, @(d) QH_SC  *d
                    'OTFS', 'otfs', false, @(d) QH_OTFS*d };
        outfile = 'compute_mcrb_mc_true_paper_results.mat';

    case 'paperv2'
        %% Paper grid extended in TIME so Doppler is observable (cf. run_fim_avg
        %% 'paperv2'): same B and N, but M=16 slots -> frame M*Ts = 12.8 us.
        %% At M=1 the Doppler columns are structurally null and nu is NOT
        %% estimable, hence M>1 here.  nu_max = half a Doppler bin, i.e. a
        %% sensing-scale target Doppler (vehicular ~500 Hz is invisible in 12.8 us).
        B = 160e6;  N = 128;  M = 16;  deltaf = B/N;  Ts = 1/deltaf;
        nTrials = 5000;                  % ICI kernel makes getAB ~5-10x slower;
                                         % with use_crn the few-% gaps resolve here
        p_list  = 0:0.002:0.1;           % fine 16-PSK BER sweep (low-BER band)
        store_all = false;
        kappa_thresh = 1e4;
        pilot_frac   = 0.20;

        %% Channel: SAME delays as 'paper', plus a Doppler shift on each path.
        L_paths = 10;  bin_t = 1/(N*deltaf);            % 1/B  (<-- switch L here)
        tau_max = (N-1)*bin_t;  min_sep = 1.0;
        rng(20260514);                                  % same delay draw as 'paper'
        usable  = tau_max - (L_paths-1)*min_sep*bin_t;
        u       = sort(usable*rand(L_paths,1));
        tau     = u + min_sep*bin_t*(0:L_paths-1).';    % TRUE delays [s]
        nu_max  = 0.5 / (M*Ts);                         % half a Doppler bin [Hz]
        rng(20260707);                                  % independent Doppler draw
        nu      = nu_max * (2*rand(L_paths,1) - 1);     % TRUE Dopplers [Hz]
        alpha   = ones(L_paths,1);
        chan_desc = sprintf('paperv2 delay-Doppler (nu_max=%.1f kHz)', nu_max/1e3);
        mod_type  = '16psk';

        %% Precoders Q^H = F_N U, applied PER OFDM SYMBOL: tfun(x) = QH*x on the
        %% N x M grid multiplies column-wise, i.e. precodes each symbol's N
        %% subcarriers (natural extension of the paper's single-block U).
        N1p = 32;  N2p = 4;
        FN  = fft(eye(N))  / sqrt(N);
        FN1 = fft(eye(N1p)) / sqrt(N1p);
        QH_OFDM = eye(N);
        QH_SC   = FN;
        QH_OTFS = FN * kron(FN1', eye(N2p));
        % Doppler => ICI: subcarriers are no longer orthogonal, so use_ici = TRUE.
        configs = { 'OFDM', 'ofdm', true, @(x) QH_OFDM*x
                    'SC',   'sc',   true, @(x) QH_SC  *x
                    'OTFS', 'otfs', true, @(x) QH_OTFS*x };
        outfile = 'compute_mcrb_mc_true_paperv2_results.mat';

    otherwise
        error('Unknown scenario "%s" (use ''eva'', ''paper'' or ''paperv2'')', scenario);
end

P  = numel(alpha);
nE = 4*P;
eta_bar = packEta(alpha, tau, nu);      % ground-truth parameter (4P x 1)
idx_tau = 4*(0:P-1) + 3;
idx_nu  = 4*(0:P-1) + 4;
alpha_ref2 = abs(alpha(1))^2;           % SNR reference power (for plot_mcrb)
snr_db_ref = 0;                         % SNR at the saved sigma^2_ref = 1
nC = size(configs,1);
np = numel(p_list);

% Delay-only channels (nu == 0, e.g. the paper scenario at M=1) have a
% structurally null Doppler block in the 4P FIM, so cond(full A) = Inf even when
% the delays are perfectly resolvable.  Restrict the conditioning check to the
% observable (alpha, tau) parameters, so kappa flags genuine delay/gain
% unresolvability rather than the intentionally-absent Doppler directions.
estimate_doppler = any(nu ~= 0);
active = true(nE,1);
if ~estimate_doppler, active(idx_nu) = false; end

% ---- Pilot-only OFDM CRB (paper scenario) ------------------------------
% A fraction pilot_frac of the N subcarriers carry KNOWN comb pilots (spread
% across the band to keep delay resolution).  Because the pilots are known at
% the RX there is NO data mismatch: the pilot-only bound is a pure CRB (only
% the sigma^2*G term, no detection-error floor F).  It uses fewer subcarriers
% so it starts higher than the full data-aided CRB, but it keeps falling with
% SNR -- so above a crossover BER/SNR it beats the floor-limited data-aided
% MCRB.  Deterministic (data-independent for unit-modulus pilots): the FIM
% depends on the pilots only through |Xpil|^2 = 1 on the comb.  Computed once.
if any(strcmpi(scenario,{'paper','paperv2'}))
    Npil      = round(pilot_frac*N);
    pilot_idx = unique(round(linspace(1, N, Npil)));    % comb across the band
    % Comb in frequency on ALL M symbols, so the pilots observe delay AND
    % (when M>1) Doppler.  Unit-modulus pilots; data cells zeroed.
    Xpil      = zeros(N, M);  Xpil(pilot_idx,:) = 1;
    if estimate_doppler                                 % Doppler => ICI kernel
        Jpil = getFIM_OFDM_multi_ICI  (Xpil, deltaf, Ts, 1, alpha, tau, nu);
    else
        Jpil = getFIM_OFDM_multi_noICI(Xpil, deltaf, Ts, 1, alpha, tau, nu);
    end
    Cpil      = stab_quad(Jpil, Jpil, active);          % = Jpil^{-1} (robust pinv)
    trGpil    = sum(real(diag(Cpil(idx_tau, idx_tau))));% tau-trace CRB @ sigma^2=1
    trGpiln   = sum(real(diag(Cpil(idx_nu,  idx_nu )))); % nu-trace CRB @ sigma^2=1
    fprintf('Pilot-only OFDM: %d/%d comb pilots x %d symbols (%.0f%%)\n', ...
            numel(pilot_idx), N, M, 100*pilot_frac);
    fprintf('  tau-trace CRB @sig2=1 = %.4e s^2  (range %.4f m^2)\n', trGpil, (3e8/2)^2*trGpil);
    if estimate_doppler
        fprintf('  nu -trace CRB @sig2=1 = %.4e Hz^2\n', trGpiln);
    end
    fprintf('\n');
else
    pilot_idx = [];  trGpil = NaN;  trGpiln = NaN;  pilot_frac = NaN;
end

%% Storage: averages always; per-iteration arrays only when store_all (they are
%% nE^2 * nTrials * np * nC and blow up for the paper P=40 / nE=160 grid).
if store_all
    LB_all = zeros(nE,nE,nTrials,np,nC);  A_all = zeros(nE,nE,nTrials,np,nC);
    B_all  = zeros(nE,nE,nTrials,np,nC);  G_all = zeros(nE,nE,nTrials,np,nC);
    F_all  = zeros(nE,nE,nTrials,np,nC);  BB_all = zeros(nE,nE,nTrials,np,nC);
    eta0_all = zeros(nE,nTrials,np,nC);
end
LB_bar = zeros(nE,nE,np,nC);  A_bar = zeros(nE,nE,np,nC);  B_bar = zeros(nE,nE,np,nC);
G_bar  = zeros(nE,nE,np,nC);  F_bar = zeros(nE,nE,np,nC);  BB_bar = zeros(nE,nE,np,nC);
eta0_bar = zeros(nE,np,nC);
kappa_max = zeros(np,nC);     % worst-over-trials Jacobi-scaled cond(A) per (p,config)
nvalid    = nTrials*ones(np,nC);  % trials kept by the Omega_N (kappa<=thresh) filter
cond_warn = 1e8;              % above this, A^-1 (and the bound) is not trustworthy
if ~store_all
    % Per-trial tau-traces of G (CRB part) and F (mismatch floor), so the plot
    % can take a robust MEDIAN over trials -- the mean of F is heavy-tailed for
    % the spread waveforms (SC/OTFS) even when A is well-conditioned.
    trG_all = zeros(nTrials, np, nC);    % tau-block traces
    trF_all = zeros(nTrials, np, nC);
    trGn_all = zeros(nTrials, np, nC);   % nu-block traces (paperv2: Doppler MCRB)
    trFn_all = zeros(nTrials, np, nC);
end

fprintf('MC lower bound at TRUE parameter (eta_0 = eta_bar)  LB(sigma2) = sigma2*G + F + BB,  BB=0\n');
fprintf('scenario=%s (%s), mod=%s, N=%d M=%d P=%d, %d trials, p in [%s], affine in SNR (sigma2_ref=1)\n', ...
        scenario, chan_desc, mod_type, N, M, P, nTrials, num2str(p_list));

nWork = nC*np*nTrials;   doneWork = 0;   tAll = tic;   % overall progress / ETA

for ci = 1:nC
    name = configs{ci,1};  use_ici = configs{ci,3};  tfun = configs{ci,4};
    tic;
    for ip = 1:np
        pp = p_list(ip);

        if store_all
            % ---- per-trial storage path (small nE, e.g. EVA scenario) --------
            % Keeps every realization (indexed by the parfor variable t) so the
            % _all archives can be saved; memory ~ nE^2 * nTrials.
            LB_t = zeros(nE,nE,nTrials);  A_t = zeros(nE,nE,nTrials);  B_t = zeros(nE,nE,nTrials);
            G_t  = zeros(nE,nE,nTrials);  F_t = zeros(nE,nE,nTrials);  BB_t = zeros(nE,nE,nTrials);
            eta0_t = zeros(nE,nTrials);   kappa_t = zeros(1,nTrials);
            parfor t = 1:nTrials
                [LB,A,B,G,F,BB,eta0,kappa] = mcrb_trial(mod_type, pp, tfun, N, M, ...
                    deltaf, Ts, eta_bar, alpha, tau, nu, use_ici, active);
                LB_t(:,:,t)=LB; A_t(:,:,t)=A; B_t(:,:,t)=B;
                G_t(:,:,t)=G;   F_t(:,:,t)=F; BB_t(:,:,t)=BB;
                eta0_t(:,t)=eta0;  kappa_t(t)=kappa;
            end
            LB_all(:,:,:,ip,ci)=LB_t; A_all(:,:,:,ip,ci)=A_t; B_all(:,:,:,ip,ci)=B_t;
            G_all(:,:,:,ip,ci)=G_t;   F_all(:,:,:,ip,ci)=F_t; BB_all(:,:,:,ip,ci)=BB_t;
            eta0_all(:,:,ip,ci)=eta0_t;
            LB_bar(:,:,ip,ci)=mean(LB_t,3); A_bar(:,:,ip,ci)=mean(A_t,3);
            B_bar(:,:,ip,ci)=mean(B_t,3);   G_bar(:,:,ip,ci)=mean(G_t,3);
            F_bar(:,:,ip,ci)=mean(F_t,3);   BB_bar(:,:,ip,ci)=mean(BB_t,3);
            eta0_bar(:,ip,ci)=mean(eta0_t,2);
            kappa_max(ip,ci)=max(kappa_t);
        else
            % ---- reduction path (large nE, e.g. paper): accumulate sums in the
            % parfor so memory is O(nE^2) irrespective of nTrials -> trials can be
            % raised freely for smoother averages.  Only the _bar means are kept.
            LBs=zeros(nE); As=zeros(nE); Bs=zeros(nE);
            Gs=zeros(nE); Fs=zeros(nE); BBs=zeros(nE);
            eta0s=zeros(nE,1); kap=0; nval=0;
            trG_t = zeros(nTrials,1);  trF_t = zeros(nTrials,1);   % per-trial tau-traces
            trGn_t = zeros(nTrials,1); trFn_t = zeros(nTrials,1);  % per-trial nu-traces
            it_par = idx_tau;  in_par = idx_nu;                    % broadcast indices
            parfor t = 1:nTrials
                if use_crn
                    % Key the stream on (p-index, trial) ONLY -- not on the config --
                    % so every waveform sees identical data + identical bit flips,
                    % making the OFDM/SC/OTFS comparison paired.
                    rs = RandStream('Threefry', 'Seed', crn_seed);
                    rs.Substream = (ip-1)*nTrials + t;
                    RandStream.setGlobalStream(rs);
                end
                [LB,A,B,G,F,BB,eta0,kappa] = mcrb_trial(mod_type, pp, tfun, N, M, ...
                    deltaf, Ts, eta_bar, alpha, tau, nu, use_ici, active);
                w = double(kappa <= kappa_thresh);   % Omega_N indicator (1 keep, 0 drop)
                LBs=LBs+w*LB; As=As+w*A; Bs=Bs+w*B; Gs=Gs+w*G; Fs=Fs+w*F; BBs=BBs+w*BB;
                eta0s=eta0s+w*eta0; kap=max(kap,kappa); nval=nval+w;
                trG_t(t) = sum(real(diag(G(it_par,it_par))));   % all trials (median is robust)
                trF_t(t) = sum(real(diag(F(it_par,it_par))));
                trGn_t(t) = sum(real(diag(G(in_par,in_par))));
                trFn_t(t) = sum(real(diag(F(in_par,in_par))));
            end
            nrm = max(nval,1);              % average over the retained (well-cond) trials
            LB_bar(:,:,ip,ci)=LBs/nrm; A_bar(:,:,ip,ci)=As/nrm;
            B_bar(:,:,ip,ci)=Bs/nrm;   G_bar(:,:,ip,ci)=Gs/nrm;
            F_bar(:,:,ip,ci)=Fs/nrm;   BB_bar(:,:,ip,ci)=BBs/nrm;
            eta0_bar(:,ip,ci)=eta0s/nrm; kappa_max(ip,ci)=kap; nvalid(ip,ci)=nval;
            trG_all(:,ip,ci)=trG_t;   trF_all(:,ip,ci)=trF_t;
            trGn_all(:,ip,ci)=trGn_t; trFn_all(:,ip,ci)=trFn_t;
        end
        doneWork = doneWork + nTrials;
        print_progress(doneWork, nWork, tAll);
    end
    tr0 = sum(real(diag(LB_bar(idx_tau,idx_tau,1,ci))));
    trL = sum(real(diag(LB_bar(idx_tau,idx_tau,np,ci))));
    kmax = max(kappa_max(:,ci));  nbad = nnz(kappa_max(:,ci) > cond_warn);
    fprintf('  %-12s done (%.0f s) : tr{LB_tau}@SNR=0dB  p=%.1f -> %.3e ,  p=%.1f -> %.3e  [s^2]\n', ...
            name, toc, p_list(1), tr0, p_list(end), trL);
    if nbad > 0
        fprintf('    !! cond(A) up to %.1e at %d/%d p-points (> %.0e): bound UNTRUSTWORTHY there (unresolvable targets)\n', ...
                kmax, nbad, np, cond_warn);
    else
        fprintf('    max cond(A) = %.1e (well-posed)\n', kmax);
    end
    if isfinite(kappa_thresh)
        frac = 100 * nvalid(:,ci) / nTrials;
        fprintf('    Omega_N (kappa<=%.0e): retained %.1f%%-%.1f%% of trials across p\n', ...
                kappa_thresh, min(frac), max(frac));
    end
end

common_vars = {'LB_bar','A_bar','B_bar','G_bar','F_bar','BB_bar','eta0_bar','eta_bar', ...
     'configs','p_list','snr_db_ref','alpha_ref2','N','M','P','nTrials', ...
     'alpha','tau','nu','deltaf','Ts','idx_tau','idx_nu', ...
     'kappa_max','cond_warn','nvalid','kappa_thresh','scenario','mod_type','chan_desc', ...
     'trGpil','trGpiln','pilot_idx','pilot_frac','estimate_doppler'};
if store_all
    save(fullfile(resdir, outfile), ...
         'LB_all','A_all','B_all','G_all','F_all','BB_all','eta0_all', ...
         common_vars{:}, '-v7.3');
else
    save(fullfile(resdir, outfile), 'trG_all','trF_all','trGn_all','trFn_all', ...
         common_vars{:}, '-v7.3');
end
fprintf('\nSaved averaged LB, A, B, G, F, BB, eta_0 (and eta_bar) to results/%s  (store_all=%d)\n', ...
        outfile, store_all);
warning('on','MATLAB:nearlySingularMatrix');

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

function [Y, kappa] = stab_quad(A, C, active)
% Y = A^-1 C A^-1 with a Jacobi (diagonal, magnitude-based) preconditioner, made
% robust to a rank-deficient A via the Moore-Penrose pseudoinverse (no backslash,
% so no nearlySingularMatrix warning).  For resolvable targets As is well
% conditioned and pinv == inv, reproducing the previous result.  For the
% unresolvable full-EVA case As has near-null directions: pinv truncates them,
% keeping Y finite; kappa = cond(As) flags that those entries are not meaningful.
%
% active (optional logical, default all-true): parameters to include in the
% conditioning number.  Excluding a structurally null block (e.g. the Doppler
% directions of a delay-only channel) keeps kappa meaningful instead of Inf.
    d  = 1 ./ sqrt(abs(diag(A)));
    d(~isfinite(d)) = 1;                 % guard against a zero FIM diagonal
    D  = diag(d);
    As = D*A*D;  Cs = D*C*D;
    As = (As + As.')/2;                  % A is real-symmetric by construction
    Ai = pinv(As);                       % robust inverse of the scaled matrix
    Y  = real(D*(Ai*Cs*Ai)*D);
    if nargout > 1
        if nargin < 3 || isempty(active), active = true(size(A,1),1); end
        kappa = cond(As(active, active));
    end
end

function eta = packEta(alpha, tau, nu)
    eta = reshape([real(alpha(:)).'; imag(alpha(:)).'; tau(:).'; nu(:).'], [], 1);
end

function [LB,A,B,G,F,BB,eta0,kappa] = mcrb_trial(mod_type, pp, tfun, N, M, ...
        deltaf, Ts, eta_bar, alpha, tau, nu, use_ici, active)
% One Monte-Carlo realization of the true-parameter MCRB pieces at sigma^2 = 1.
% Draws data, applies the waveform transform, corrupts it (detection errors),
% and returns A, B, and the affine bound components G = A^-1 Jtilde A^-1,
% F = A^-1 (b b^T) A^-1, BB = 0 (eta_0 = eta_bar), LB = G + F + BB.
    % The client-side warning('off',...) does NOT propagate to parfor workers.
    warning('off','MATLAB:nearlySingularMatrix');
    warning('off','MATLAB:singularMatrix');
    d  = gen_data(N, M, mod_type);
    dh = corrupt_data(d, pp, mod_type);

    eta0 = eta_bar;                     % pseudo-true assumed equal to true
    Xtf  = tfun(d);                     % waveform transform / precoder (scenario-set)
    Xhtf = tfun(dh);
    [a0, t0, n0] = unpackEta(eta0);
    [A, B, Jt, b] = getAB_OFDM_core(Xtf, Xhtf, deltaf, Ts, 1, ...
                                    a0, t0, n0, use_ici, alpha, tau, nu);

    [G, kappa] = stab_quad(A, Jt, active);   % A^-1 Jtilde A^-1 (+ cond over active block)
    F  = stab_quad(A, b*b.');                % A^-1 (b b^T) A^-1
    bias = eta_bar - eta0;                   % = 0 (eta_0 = eta_bar)
    BB = bias*bias.';                        % pseudo-true bias outer product (=0)
    LB = G + F + BB;                         % = A^-1 B A^-1 + BB  (at sigma2 = 1)
    LB = (LB + LB.')/2;
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

% ----- modulation dispatch (scenario-selected via mod_type) -------------
function X = gen_data(N, M, mod_type)
% i.i.d. unit-power data symbols on the N x M grid.
    switch lower(mod_type)
        case 'qpsk'
            X = qpsk(N, M);
        case '16psk'
            k = randi(16, N, M) - 1;
            X = exp(1j*2*pi*k/16);
        otherwise
            error('Unknown mod_type "%s" (use ''qpsk'' or ''16psk'')', mod_type);
    end
end

function Xh = corrupt_data(X, p, mod_type)
% Detection-error / mismatch model: crossover probability p.
    switch lower(mod_type)
        case 'qpsk'
            Xh = bsc_qpsk(X, p);
        case '16psk'
            Xh = bsc_16psk(X, p);
        otherwise
            error('Unknown mod_type "%s" (use ''qpsk'' or ''16psk'')', mod_type);
    end
end

function Xh = bsc_16psk(X, p)
% 16-PSK corruption: independent BSC at rate p on the 4 Gray-coded bits of each
% symbol's phase index -- the direct 16-ary analogue of bsc_qpsk's 2-bit BSC.
% A single bit flip on the Gray label moves to an adjacent (or nearby) phase.
    k = mod(round(angle(X)*16/(2*pi)), 16);   % phase index 0..15
    g = bitxor(k, floor(k/2));                 % Gray encode
    for bit = 0:3
        flip = double(rand(size(X)) < p) * 2^bit;
        g = bitxor(g, flip);
    end
    kd = g;                                     % Gray decode (4-bit)
    kd = bitxor(kd, floor(kd/2));
    kd = bitxor(kd, floor(kd/4));
    Xh = exp(1j*2*pi*mod(kd,16)/16);
end
