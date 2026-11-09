% compute_mcrb_mc.m
% Monte-Carlo misspecified lower bound at the PSEUDO-TRUE parameter.
%
% Two scenarios via the `scenario` switch:
%   'eva'   : original 48x24 reduced/full-EVA, QPSK, physical OFDM/CP-OTFS,
%             {OFDM,OTFS}x{no-ICI,ICI}, swept over the BSC rate p.  Untouched.
%   'paper' : Fig.1 scenario matching compute_mcrb_mc_true -- 128x1 delay-only
%             L-target channel, 16-PSK, precoded OFDM/SC/OTFS (Q^H = F_N U) fed
%             to the SAME pseudo-true solver / OFDM core.  Adds the pilot-only
%             OFDM CRB curve.  Stores per-trial tau-traces for a robust MEDIAN.
%
% For every realization (X, Xhat):
%   1. eta_0 = pseudo-true parameter (pseudoTrue_OFDM / pseudoTrue_CP_OTFS),
%      the projection argmin_eta || mu(X,eta_bar) - mu_tilde(Xhat,eta) ||.
%   2. A(eta_0), B(eta_0), Jtilde(eta_0), b(eta_0) at sigma^2 = 1.
%   3. LB(sigma^2) = A^{-1} B A^{-1} + (eta_bar-eta_0)(eta_bar-eta_0)^T
%                  = sigma^2*G + F + BB,   evaluated at sigma^2 = 1, with
%        G  = A^{-1} Jtilde A^{-1}        (CRB-like, scales with noise)
%        F  = A^{-1} (b b^T) A^{-1}       (score-bias floor, sigma^2-indep.)
%        BB = (eta_bar-eta_0)(eta_bar-eta_0)^T  (pseudo-true bias, sigma^2-indep.)
%
% Loop order: trials are the parfor dimension; p is swept INSIDE each trial with
% d fixed across p, so the pseudo-true lsqnonlin is WARM-STARTED from the
% previous p's solution.  At p = 0: eta_0 = eta_bar, A = -J, b = 0, BB = 0,
% F = 0, G = J^{-1}, so LB = sigma^2 J^{-1} (the CRB).
clear; clc;
thisdir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisdir, 'utils'));
resdir = fullfile(thisdir, 'results');
warning('off','MATLAB:nearlySingularMatrix');
rng(1);   % parfor workers draw from their own streams (statistically equivalent).

%% ===== Scenario selector =================================================
scenario   = 'paper';
use_varpro = true;       % pseudo-true solver: varpro (fast) vs full 4P lsqnonlin

% Common random numbers across waveforms (paired comparison).  Keyed on the
% TRIAL only (d is drawn once per trial and held fixed across p), so every
% config sees the SAME data and the SAME per-p bit flips -> the OFDM/SC/OTFS
% difference is paired and the MC noise largely cancels.  Essential: the
% paperv2 gaps are only ~1-5%, which independent draws cannot resolve.
use_crn  = true;
crn_seed = 20260707;

switch lower(scenario)
    case 'eva'
        %% Parameters
        deltaf = 15e3;  Ts = 1/deltaf;
        N = 48;   M = 24;                % reduced grid (see original notes)
        nTrials = 1000;
        p_list  = 0:0.02:0.5;
        store_all    = true;             % keep full per-iteration arrays (small nE)
        kappa_thresh = inf;              % no Omega_N filtering (preserve original)
        pilot_frac   = NaN;              % no pilot curve in eva
        use_full_eva = false;            % false = reduced EVA (3 taps), true = 9 paths

        fc = 2e9;  v = 120/3.6;  c = 3e8;  numax = v*fc/c;   % max Doppler ~222.2 Hz
        if use_full_eva
            tau_ns = [0;   30;  150;  310;  370;  710; 1090; 1730; 2510];   % [ns]
            pow_db = [0; -1.5; -1.4; -3.6; -0.6; -9.1; -7.0; -12.0; -16.9]; % [dB]
            theta  = [0.5; 2.0; 3.5; 0.9; 1.6; 4.2; 2.7; 5.1; 3.0];
            phase  = [0.7; 1.1; -0.5; 0.3; -1.2; 2.0; -0.8; 1.5; 0.2];
        else
            tau_ns = [0;  370;  150];
            pow_db = [0; -0.6; -1.4];
            theta  = [0.5; 2.0; 3.5];
            phase  = [0.7; 1.1; -0.5];
        end
        tau   = tau_ns*1e-9;
        nu    = numax*cos(theta);
        alpha = 10.^(pow_db/20) .* exp(1j*phase);
        chan_desc = 'full EVA';  if ~use_full_eva, chan_desc = 'reduced EVA'; end
        mod_type  = 'qpsk';
        % {name, pseudo-true solver key, use_ici, transform d -> TF symbols}
        configs = { 'OFDM no-ICI', 'ofdm_direct', false, @(x) x
                    'OFDM ICI',    'ofdm_direct', true,  @(x) x
                    'OTFS no-ICI', 'otfs_dd',     false, @(x) isfft_u(x)
                    'OTFS ICI',    'otfs_dd',     true,  @(x) isfft_u(x) };
        outfile = 'compute_mcrb_mc_results.mat';

    case 'paper'
        %% Parameters (paper Fig.1: B=160 MHz, N=128 single block, delay-only)
        B = 160e6;  N = 128;  M = 1;  deltaf = B/N;  Ts = 1/deltaf;
        nTrials = 2000;                  % pseudo-true is slow but varpro is warm-started
        p_list  = [0 0.1];           % fine 16-PSK BER sweep (low-BER band)
        store_all    = false;            % per-trial tau-traces only (median plot)
        kappa_thresh = 1e4;              % Omega_N well-conditioned threshold (reporting)
        pilot_frac   = 0.20;             % pilot-only OFDM CRB uses this subcarrier frac

        %% Paper channel: L delay-only targets, SAME rng draw as the _true script.
        L_paths = 40;  bin_t = 1/(N*deltaf);            % 1/B  (<-- switch L here)
        tau_max = (N-1)*bin_t;  min_sep = 1.0;
        rng(20260514);
        usable  = tau_max - (L_paths-1)*min_sep*bin_t;
        u       = sort(usable*rand(L_paths,1));
        tau     = u + min_sep*bin_t*(0:L_paths-1).';    % TRUE delays [s]
        nu      = zeros(L_paths,1);                      % delay-only
        alpha   = ones(L_paths,1);                       % unit gains
        chan_desc = 'paper delay-only';
        mod_type  = '16psk';

        %% Precoders Q^H = F_N U: OFDM U=F_N^H, SC U=I, OTFS U=F_N1^H (x) I_N2.
        N1p = 32;  N2p = 4;
        FN  = fft(eye(N))  / sqrt(N);
        FN1 = fft(eye(N1p)) / sqrt(N1p);
        QH_OFDM = eye(N);
        QH_SC   = FN;
        QH_OTFS = FN * kron(FN1', eye(N2p));
        % All three use the OFDM pseudo-true solver on the PRECODED TF symbols.
        configs = { 'OFDM', 'ofdm_direct', false, @(x) QH_OFDM*x
                    'SC',   'ofdm_direct', false, @(x) QH_SC  *x
                    'OTFS', 'ofdm_direct', false, @(x) QH_OTFS*x };
        outfile = 'compute_mcrb_mc_paper_results.mat';

    case 'paperv2'
        %% Paper grid extended in TIME so Doppler is observable (same definition
        %% as compute_mcrb_mc_true 'paperv2' / run_fim_avg).  M=16 slots ->
        %% frame M*Ts = 12.8 us.  At M=1 the Doppler columns are structurally
        %% null and nu is NOT estimable, hence M>1 here.  nu_max = half a
        %% Doppler bin (sensing-scale; vehicular ~500 Hz is invisible in 12.8 us).
        B = 160e6;  N = 128;  M = 16;  deltaf = B/N;  Ts = 1/deltaf;
        nTrials = 2000;                  % pseudo-true solve + ICI kernel (~75 min);
                                         % use_crn makes 2000 enough for the 1-5% gaps
        p_list  = 0:0.002:0.1;
        store_all    = false;
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

        %% Precoders Q^H = F_N U, applied PER OFDM SYMBOL (tfun(x)=QH*x on N x M
        %% multiplies column-wise).  Doppler => ICI: use_ici = TRUE.
        N1p = 32;  N2p = 4;
        FN  = fft(eye(N))  / sqrt(N);
        FN1 = fft(eye(N1p)) / sqrt(N1p);
        QH_OFDM = eye(N);
        QH_SC   = FN;
        QH_OTFS = FN * kron(FN1', eye(N2p));
        configs = { 'OFDM', 'ofdm_direct', true, @(x) QH_OFDM*x
                    'SC',   'ofdm_direct', true, @(x) QH_SC  *x
                    'OTFS', 'ofdm_direct', true, @(x) QH_OTFS*x };
        outfile = 'compute_mcrb_mc_paperv2_results.mat';

    case 'evace'
        %% Match mcrb_ce_ici.m EXACTLY (channel + grid), but with the ABSTRACT BSC
        %% error model (swept over p) instead of a realistic detector.  Purpose:
        %% overlay the BSC MCRB against the realistic-detection MCRB of
        %% mcrb_ce_ici at matched symbol-error rate.  N=512, M=128, reduced EVA
        %% (3 taps), 4-QAM, ICI on, physical OFDM / CP-OTFS, and -- like
        %% mcrb_ce_ici -- a RANDOM channel (Doppler + phase) per trial.
        deltaf = 15e3;  Ts = 1/deltaf;
        N = 512;  M = 128;
        nTrials = 50;                        % match mcrb_ce_ici
        p_list  = 0:0.02:0.5;                % BSC crossover; SER_qpsk=1-(1-p)^2
        store_all    = false;                % reduction path + per-trial medians
        kappa_thresh = inf;                  % well-conditioned (512x128, 3 taps)
        pilot_frac   = NaN;
        fc = 4e9;  c = 3e8;  v = 120/3.6;  numax = v*fc/c;      % nu_max = 444 Hz
        delay_taps = round([0 150 370]*1e-9 * N * deltaf);      % [0 1 3] taps
        tau        = (delay_taps/(N*deltaf)).';
        amp_chan   = (10.^([0 -1.4 -0.6]/20)).';  amp_chan = amp_chan/norm(amp_chan);
        numax_chan = numax;  rand_chan = true;                  % redrawn per trial
        nu    = (numax*cos(2*pi*rand(1,numel(tau))-pi)).';      % placeholder (for nE)
        alpha = amp_chan .* exp(1j*2*pi*rand(numel(tau),1));
        chan_desc = 'reduced EVA 512x128, random Doppler (matches mcrb_ce_ici)';
        mod_type  = 'qpsk';                                     % = 4-QAM unit power
        configs = { 'OFDM ICI',    'ofdm_direct', true, @(x) x
                    'CP-OTFS ICI', 'otfs_dd',     true, @(x) isfft_u(x) };
        outfile = 'compute_mcrb_mc_evace_results.mat';

    otherwise
        error('Unknown scenario "%s" (use ''eva'', ''paper'', ''paperv2'' or ''evace'')', scenario);
end

% Channel randomization per trial (evace matches mcrb_ce_ici; others use the
% fixed channel set above).
if ~exist('rand_chan','var'), rand_chan = false;  amp_chan = [];  numax_chan = 0;  end

P  = numel(alpha);
nE = 4*P;
eta_bar = packEta(alpha, tau, nu);
idx_tau = 4*(0:P-1) + 3;
idx_nu  = 4*(0:P-1) + 4;
alpha_ref2 = abs(alpha(1))^2;
snr_db_ref = 0;
nC = size(configs,1);
np = numel(p_list);

% Delay-only channels (nu == 0, M=1) have a null Doppler block; restrict the
% conditioning check to the observable (alpha, tau) parameters.
estimate_doppler = any(nu ~= 0);
active = true(nE,1);
if ~estimate_doppler, active(idx_nu) = false; end

% ---- Pilot-only OFDM CRB (paper scenario) ------------------------------
% Known comb pilots on a fraction of subcarriers -> no data mismatch -> pure CRB
% (no floor).  Deterministic (data-independent for unit-modulus pilots).
if any(strcmpi(scenario,{'paper','paperv2'}))
    Npil      = round(pilot_frac*N);
    pilot_idx = unique(round(linspace(1, N, Npil)));
    % Comb in frequency on ALL M symbols, so pilots observe delay AND (M>1) Doppler.
    Xpil      = zeros(N, M);  Xpil(pilot_idx,:) = 1;
    if estimate_doppler                                 % Doppler => ICI kernel
        Jpil = getFIM_OFDM_multi_ICI  (Xpil, deltaf, Ts, 1, alpha, tau, nu);
    else
        Jpil = getFIM_OFDM_multi_noICI(Xpil, deltaf, Ts, 1, alpha, tau, nu);
    end
    Cpil      = stab_quad(Jpil, Jpil, active);
    trGpil    = sum(real(diag(Cpil(idx_tau, idx_tau))));
    trGpiln   = sum(real(diag(Cpil(idx_nu,  idx_nu ))));
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

%% Storage: averages/per-iteration arrays only for eva (store_all); the paper
%% path keeps per-trial tau-traces of G, F, BB (for the robust median).
if store_all
    LB_all=zeros(nE,nE,nTrials,np,nC); A_all=zeros(nE,nE,nTrials,np,nC);
    B_all =zeros(nE,nE,nTrials,np,nC); G_all=zeros(nE,nE,nTrials,np,nC);
    F_all =zeros(nE,nE,nTrials,np,nC); BB_all=zeros(nE,nE,nTrials,np,nC);
    eta0_all=zeros(nE,nTrials,np,nC);
    LB_bar=zeros(nE,nE,np,nC); A_bar=zeros(nE,nE,np,nC); B_bar=zeros(nE,nE,np,nC);
    G_bar =zeros(nE,nE,np,nC); F_bar=zeros(nE,nE,np,nC); BB_bar=zeros(nE,nE,np,nC);
    eta0_bar=zeros(nE,np,nC);
else
    % tau-block traces
    trG_all=zeros(nTrials,np,nC); trF_all=zeros(nTrials,np,nC); trBB_all=zeros(nTrials,np,nC);
    % nu-block traces (paperv2: Doppler MCRB)
    trGn_all=zeros(nTrials,np,nC); trFn_all=zeros(nTrials,np,nC); trBBn_all=zeros(nTrials,np,nC);
end
kappa_max = zeros(np,nC);
nvalid    = nTrials*ones(np,nC);
cond_warn = 1e8;

solver = 'VARPRO (tau,nu; gains closed-form)';  if ~use_varpro, solver = 'full 4P'; end
fprintf('MC pseudo-true lower bound  LB(sigma2) = sigma2*G + F + BB,  BB=(eta_bar-eta_0)(.)^T\n');
fprintf('scenario=%s (%s), mod=%s, N=%d M=%d P=%d, %d trials, p in [%s], warm-started, %s\n', ...
        scenario, chan_desc, mod_type, N, M, P, nTrials, num2str(p_list), solver);

nWork = nC*np*nTrials;   doneWork = 0;   tAll = tic;

for ci = 1:nC
    name = configs{ci,1};  ptkey = configs{ci,2};  use_ici = configs{ci,3};  tfun = configs{ci,4};
    tic;

    if store_all
        % ---- EVA path: full per-iteration storage (slabs over p, sliced by t) ----
        LB_ct=zeros(nE,nE,np,nTrials); A_ct=zeros(nE,nE,np,nTrials); B_ct=zeros(nE,nE,np,nTrials);
        G_ct =zeros(nE,nE,np,nTrials); F_ct=zeros(nE,nE,np,nTrials); BB_ct=zeros(nE,nE,np,nTrials);
        eta0_ct=zeros(nE,np,nTrials);  kappa_ct=zeros(np,nTrials);
        parfor t = 1:nTrials
            if use_crn
                rs = RandStream('Threefry','Seed',crn_seed);  rs.Substream = t;
                RandStream.setGlobalStream(rs);
            end
            [~,~,~,~,~,~,kappa_p,eta0_p,LB_p,A_p,B_p,G_p,F_p,BB_p] = pt_trial( ...
                N,M,mod_type,ptkey,tfun,p_list,deltaf,Ts,eta_bar,alpha,tau,nu,use_ici,active,idx_tau,idx_nu,use_varpro, rand_chan,amp_chan,numax_chan);
            LB_ct(:,:,:,t)=LB_p; A_ct(:,:,:,t)=A_p; B_ct(:,:,:,t)=B_p;
            G_ct(:,:,:,t)=G_p;   F_ct(:,:,:,t)=F_p; BB_ct(:,:,:,t)=BB_p;
            eta0_ct(:,:,t)=eta0_p;  kappa_ct(:,t)=kappa_p;
        end
        LB_all(:,:,:,:,ci)=permute(LB_ct,[1 2 4 3]); A_all(:,:,:,:,ci)=permute(A_ct,[1 2 4 3]);
        B_all(:,:,:,:,ci) =permute(B_ct,[1 2 4 3]);  G_all(:,:,:,:,ci)=permute(G_ct,[1 2 4 3]);
        F_all(:,:,:,:,ci) =permute(F_ct,[1 2 4 3]);  BB_all(:,:,:,:,ci)=permute(BB_ct,[1 2 4 3]);
        eta0_all(:,:,:,ci)=permute(eta0_ct,[1 3 2]);
        LB_bar(:,:,:,ci)=mean(LB_ct,4); A_bar(:,:,:,ci)=mean(A_ct,4); B_bar(:,:,:,ci)=mean(B_ct,4);
        G_bar(:,:,:,ci) =mean(G_ct,4);  F_bar(:,:,:,ci)=mean(F_ct,4); BB_bar(:,:,:,ci)=mean(BB_ct,4);
        eta0_bar(:,:,ci)=mean(eta0_ct,3);
        kappa_max(:,ci)=max(kappa_ct,[],2);
        tr0 = sum(real(diag(LB_bar(idx_tau,idx_tau,1,ci))));
        trL = sum(real(diag(LB_bar(idx_tau,idx_tau,np,ci))));
    else
        % ---- Paper path: per-trial tau-traces of G/F/BB only (median) ----
        trG_ct=zeros(np,nTrials); trF_ct=zeros(np,nTrials); trBB_ct=zeros(np,nTrials);
        trGn_ct=zeros(np,nTrials); trFn_ct=zeros(np,nTrials); trBBn_ct=zeros(np,nTrials);
        kappa_ct=zeros(np,nTrials);
        parfor t = 1:nTrials
            if use_crn
                % Key on the TRIAL only (not the config) -> every waveform sees
                % identical data and identical per-p flips: paired comparison.
                rs = RandStream('Threefry','Seed',crn_seed);  rs.Substream = t;
                RandStream.setGlobalStream(rs);
            end
            [trG_p,trF_p,trBB_p,trGn_p,trFn_p,trBBn_p,kappa_p] = pt_trial( ...
                N,M,mod_type,ptkey,tfun,p_list,deltaf,Ts,eta_bar,alpha,tau,nu,use_ici,active,idx_tau,idx_nu,use_varpro, rand_chan,amp_chan,numax_chan);
            trG_ct(:,t)=trG_p; trF_ct(:,t)=trF_p; trBB_ct(:,t)=trBB_p; kappa_ct(:,t)=kappa_p;
            trGn_ct(:,t)=trGn_p; trFn_ct(:,t)=trFn_p; trBBn_ct(:,t)=trBBn_p;
        end
        trG_all(:,:,ci)=trG_ct.'; trF_all(:,:,ci)=trF_ct.'; trBB_all(:,:,ci)=trBB_ct.';
        trGn_all(:,:,ci)=trGn_ct.'; trFn_all(:,:,ci)=trFn_ct.'; trBBn_all(:,:,ci)=trBBn_ct.';
        kappa_max(:,ci)=max(kappa_ct,[],2);
        nvalid(:,ci)=sum(kappa_ct<=kappa_thresh,2);
        tr0 = median(trG_all(:,1,ci)  + trF_all(:,1,ci)  + trBB_all(:,1,ci));   % sig2=1
        trL = median(trG_all(:,np,ci) + trF_all(:,np,ci) + trBB_all(:,np,ci));
    end

    doneWork = doneWork + np*nTrials;
    print_progress(doneWork, nWork, tAll);
    kmax = max(kappa_max(:,ci));  nbad = nnz(kappa_max(:,ci) > cond_warn);
    fprintf('  %-12s done (%.0f s) : tr{LB_tau}@SNR=0dB  p=%.3f -> %.3e ,  p=%.3f -> %.3e  [s^2]\n', ...
            name, toc, p_list(1), tr0, p_list(end), trL);
    if nbad > 0
        fprintf('    !! cond(A) up to %.1e at %d/%d p-points (> %.0e): bound UNTRUSTWORTHY there\n', ...
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

%% ----- Save ------------------------------------------------------------
common_vars = {'configs','p_list','snr_db_ref','alpha_ref2','N','M','P','nTrials', ...
     'alpha','tau','nu','deltaf','Ts','idx_tau','idx_nu','eta_bar', ...
     'kappa_max','cond_warn','nvalid','kappa_thresh','scenario','mod_type','chan_desc', ...
     'trGpil','trGpiln','pilot_idx','pilot_frac','estimate_doppler'};
if store_all
    save(fullfile(resdir, outfile), ...
         'LB_all','A_all','B_all','G_all','F_all','BB_all','eta0_all', ...
         'LB_bar','A_bar','B_bar','G_bar','F_bar','BB_bar','eta0_bar', ...
         common_vars{:}, '-v7.3');
else
    save(fullfile(resdir, outfile), 'trG_all','trF_all','trBB_all', ...
         'trGn_all','trFn_all','trBBn_all', common_vars{:}, '-v7.3');
end
fprintf('\nSaved results (store_all=%d) to results/%s\n', store_all, outfile);
warning('on','MATLAB:nearlySingularMatrix');

%% ---- helpers ----------------------------------------------------------
function [trG_p,trF_p,trBB_p,trGn_p,trFn_p,trBBn_p,kappa_p,eta0_p, LB_p,A_p,B_p,G_p,F_p,BB_p] = pt_trial( ...
        N,M,mod_type,ptkey,tfun,p_list,deltaf,Ts,eta_bar,alpha,tau,nu,use_ici,active,idx_tau,idx_nu,use_varpro, ...
        rand_chan,amp_chan,numax_chan)
% One Monte-Carlo trial: draw transmitted symbols d (fixed across p), then sweep
% p (warm-started pseudo-true) and return the tau- AND nu-block traces of G, F,
% BB per p (and, if requested via nargout>8, the full 4Px4P slabs for the EVA
% archive).  The caller seeds the RNG stream before calling (CRN pairing).
% rand_chan: redraw the channel (Doppler + phase, delays fixed) per trial from
% amp_chan / numax_chan -- matches mcrb_ce_ici's per-trial channel averaging.
    warning('off','MATLAB:nearlySingularMatrix');
    warning('off','MATLAB:singularMatrix');
    np = numel(p_list);  nE = numel(eta_bar);
    if rand_chan                           % drawn BEFORE the data, under the CRN stream
        Pn    = numel(tau);
        nu    = (numax_chan*cos(2*pi*rand(1,Pn)-pi)).';
        alpha = amp_chan .* exp(1j*2*pi*rand(Pn,1));
        eta_bar = packEta(alpha, tau, nu);
    end
    d   = gen_data(N, M, mod_type);        % transmitted symbols, FIXED across p
    Xtf = tfun(d);
    want_slabs = nargout > 8;

    trG_p=zeros(np,1); trF_p=zeros(np,1); trBB_p=zeros(np,1); kappa_p=zeros(np,1);
    trGn_p=zeros(np,1); trFn_p=zeros(np,1); trBBn_p=zeros(np,1);
    eta0_p=zeros(nE,np);
    if want_slabs
        LB_p=zeros(nE,nE,np); A_p=zeros(nE,nE,np); B_p=zeros(nE,nE,np);
        G_p =zeros(nE,nE,np); F_p=zeros(nE,nE,np); BB_p=zeros(nE,nE,np);
    end

    eta_init = [];                         % first p: default init = eta_bar
    for ip = 1:np
        dh   = corrupt_data(d, p_list(ip), mod_type);   % detected symbols ~ BSC(p)
        Xhtf = tfun(dh);
        switch ptkey
            case 'ofdm_direct'                          % OFDM / SC / OTFS-precoded
                eta0 = pseudoTrue_OFDM(Xtf, Xhtf, deltaf, Ts, alpha, tau, nu, use_ici, eta_init, use_varpro);
            case 'otfs_dd'                              % EVA physical CP-OTFS
                eta0 = pseudoTrue_CP_OTFS(d, dh, deltaf, Ts, alpha, tau, nu, use_ici, eta_init, use_varpro);
        end
        eta_init = eta0;                                % warm-start next (larger) p

        [a0, t0, n0] = unpackEta(eta0);
        [A, B, Jt, b] = getAB_OFDM_core(Xtf, Xhtf, deltaf, Ts, 1, ...
                                        a0, t0, n0, use_ici, alpha, tau, nu);
        [G, kappa] = stab_quad(A, Jt, active);
        F   = stab_quad(A, b*b.');
        bias = eta_bar - eta0;  BB = bias*bias.';

        trG_p(ip)  = sum(real(diag(G (idx_tau,idx_tau))));
        trF_p(ip)  = sum(real(diag(F (idx_tau,idx_tau))));
        trBB_p(ip) = sum(real(diag(BB(idx_tau,idx_tau))));
        trGn_p(ip)  = sum(real(diag(G (idx_nu,idx_nu))));
        trFn_p(ip)  = sum(real(diag(F (idx_nu,idx_nu))));
        trBBn_p(ip) = sum(real(diag(BB(idx_nu,idx_nu))));
        kappa_p(ip)= kappa;  eta0_p(:,ip)=eta0;

        if want_slabs
            LB = (G + F + BB);  LB = (LB + LB.')/2;
            LB_p(:,:,ip)=LB; A_p(:,:,ip)=A; B_p(:,:,ip)=B;
            G_p(:,:,ip)=G;   F_p(:,:,ip)=F; BB_p(:,:,ip)=BB;
        end
    end
end

function print_progress(done, total, tAll)
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
% Y = A^-1 C A^-1 via a Jacobi (magnitude) preconditioner + Moore-Penrose pinv
% (robust to a rank-deficient A).  active (optional) restricts the reported
% conditioning number to the observable parameters (skips a null Doppler block).
    d  = 1 ./ sqrt(abs(diag(A)));
    d(~isfinite(d)) = 1;
    D  = diag(d);
    As = D*A*D;  Cs = D*C*D;
    As = (As + As.')/2;
    Ai = pinv(As);
    Y  = real(D*(Ai*Cs*Ai)*D);
    if nargout > 1
        if nargin < 3 || isempty(active), active = true(size(A,1),1); end
        kappa = cond(As(active, active));
    end
end

function eta = packEta(alpha, tau, nu)
    eta = reshape([real(alpha(:)).'; imag(alpha(:)).'; tau(:).'; nu(:).'], [], 1);
end

% ----- modulation dispatch (scenario-selected via mod_type) -------------
function X = gen_data(N, M, mod_type)
    switch lower(mod_type)
        case 'qpsk'
            b = randi([0 3], N, M);
            X = (1/sqrt(2)) * ((2*bitand(b,1)-1) + 1j*(2*floor(b/2)-1));
        case '16psk'
            k = randi(16, N, M) - 1;
            X = exp(1j*2*pi*k/16);
        otherwise
            error('Unknown mod_type "%s" (use ''qpsk'' or ''16psk'')', mod_type);
    end
end

function Xh = corrupt_data(X, p, mod_type)
    switch lower(mod_type)
        case 'qpsk'
            sR = sign(real(X)); sI = sign(imag(X));
            fR = 1 - 2*(rand(size(X)) < p);
            fI = 1 - 2*(rand(size(X)) < p);
            Xh = (sR.*fR + 1j*(sI.*fI))/sqrt(2);
        case '16psk'
            Xh = bsc_16psk(X, p);
        otherwise
            error('Unknown mod_type "%s" (use ''qpsk'' or ''16psk'')', mod_type);
    end
end

function Xh = bsc_16psk(X, p)
% 16-PSK: independent BSC at rate p on the 4 Gray-coded bits of each symbol's
% phase index -- the direct 16-ary analogue of the QPSK 2-bit BSC.
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
