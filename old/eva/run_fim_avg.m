% run_fim_avg.m
% --------------------------------------------------------------------------
% Multi-target FIM/CRB driver with two independent switches:
%
%   param_mode    = 'eva' | 'paper' | 'paperv2' | 'paperdopp'   grid / system
%   channel_type  = 'EVA' | 'paper' | 'paperv2' | 'paperdopp'   channel
%
% 'paperdopp' is the DOPPLER-ONLY dual of 'paper': N=1 (single subcarrier, no
% delay dimension), M=128 slots, 40 targets with tau=0 and distinct Dopplers
% filling the Doppler window [0, 1/T).  Estimates nu only (tau known = 0), via the
% estimate_delay=false / estimate_doppler=true path in crb_multi.
%
% Paper-mode follows Liu et al. 2026 (crb_new/run_sim.m, eqs. 8b/8c):
%   - OFDM    : single transmission block, N=128 subcarriers, M=1.
%   - CP-OTFS : N1=32 subcarriers x N2=4 time slots,  N1*N2 = 128.
% (Both share the same bandwidth B = 160 MHz, hence different deltaf.)
% Doppler is treated as known (=0) in paper mode, matching the paper's
% theta = [tau, b_a, b_b] parameter vector.
%
% Paperv2-mode extends the paper setting to a delay-Doppler channel:
%   - channel 'paperv2': the SAME 40-path delay profile as 'paper'
%     (identical rng draw), but each path additionally carries a Doppler
%     shift drawn uniformly in [-nu_max, nu_max].
%   - param  'paperv2': same B = 160 MHz, N = 128 grid, but M = 16 time
%     slots so the frame (M*T = 12.8 us) actually spans enough time for
%     Doppler to be observable; nu_max is set to half a Doppler resolution
%     bin, 0.5/(M*T).  With M = 1 the Doppler FIM block is numerically
%     null and the joint (tau, nu) CRB is meaningless.
%
% EVA-mode keeps a single (N, M, deltaf) for both waveforms (our default
% N=512, M=128 CP-OTFS frame), so OFDM and CP-OTFS share the grid and the
% (I)SFFT-unitarity FIM equivalence applies directly.
% --------------------------------------------------------------------------

clear; clc;
thisdir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisdir, 'utils'));
resdir = fullfile(thisdir, 'results');
warning('off','MATLAB:nearlySingularMatrix');

% ====== Selectors ========================================================
channel_type = 'paperdopp';   % 'EVA' | 'paper' | 'paperv2' | 'paperdopp' (Doppler-only)
param_mode   = 'paperdopp';   % 'eva' | 'paper' | 'paperv2' | 'paperdopp'
modulation   = '16-PSK';      % 'QPSK' | '16-QAM' | '16-PSK'
                            % All three are zero-mean, unit-power and satisfy
                            % E[s^2] = 0 (paper Assumptions 1-3).  Kurtosis:
                            % QPSK = 1, 16-PSK = 1, 16-QAM = 1.32.  Per the
                            % paper, OFDM and CP-OTFS CRBs coincide when kappa=1
                            % and OFDM beats CP-OTFS when kappa>1 (16-QAM).

% ====== Per-waveform grid / system parameters ============================
c = 3e8;
switch lower(param_mode)
    case 'eva'
        % EVA-style: same grid for both waveforms.
        N_o = 512;  M_o = 128;  deltaf_o = 15e3;
        N_t = 512;  M_t = 128;  deltaf_t = 15e3;

        fc       = 4e9;
        velocity = 120 * 1e3 / 3600;

    case 'paper'
        % Paper-style:  N=128 subcarriers, M=4 OFDM symbols for both
        % waveforms.  Shared grid means the OFDM/CP-OTFS (I)SFFT-unitarity
        % equivalence applies and CP-OTFS is no longer rank-deficient.
        B   = 160e6;                         % total bandwidth [Hz]
        N_o = 128;  M_o = 1;    deltaf_o = B/N_o;
        % OTFS must resolve the L=40 delay-only targets that span [0,127]/B.
        % Our FIM resolves delay over the N_t DD rows (Tsamp = 1/(N_t*deltaf_t),
        % delay wrap at N_t bins), so the paper's N1=32 frame aliases delays
        % mod 32 and the FIM goes rank-deficient.  The paper instead resolves
        % delay over the full 128-length block (receiver F_128).  To match that
        % resolution we give OTFS the SAME 128-bin grid as OFDM; the single
        % frame CP (RCP-OTFS) then plays the role of the paper's CP operation.
        N_t = 128;  M_t = 1;    deltaf_t = B/N_t;
        % Paper's literal 32x4 frame (aliases 40 targets into 32 delay bins):
        % N_t = 32;  M_t = 4;    deltaf_t = B/N_t;

        fc       = 4e9;
        velocity = 0;                        % paper has no Doppler

    case 'paperv2'
        % Paper-style grid extended in time so Doppler is observable:
        % same bandwidth B = 160 MHz and N = 128 subcarriers, but M = 16
        % time slots (frame duration M*T = 12.8 us).  nu_max is set to
        % half a Doppler resolution bin, 0.5/(M*T) ~ 39 kHz, i.e. a
        % sensing-scale (radar target) Doppler rather than a vehicular
        % one -- with vehicular Doppler (~500 Hz) nothing is observable
        % in a 12.8 us frame.  Reported as an equivalent radial velocity.
        B   = 160e6;
        N_o = 128;  M_o = 16;   deltaf_o = B/N_o;
        N_t = 128;  M_t = 16;   deltaf_t = B/N_t;

        fc         = 4e9;
        nu_max_des = 0.5 / (M_o * N_o / B);  % half a Doppler bin [Hz]
        velocity   = nu_max_des * c / fc;    % equivalent velocity [m/s]

    case 'paperdopp'
        % Doppler-only DUAL of 'paper': single subcarrier (N=1, so there is no
        % delay/bandwidth dimension), M=128 time slots -> Doppler resolved over
        % 128 Doppler bins.  The paper resolves 40 delay targets over N=128 delay
        % bins; here 40 Doppler targets over M=128 Doppler bins.  Use with
        % channel_type='paperdopp' (tau=0, distinct nu per path).
        N_o = 1;    M_o = 128;   deltaf_o = 15e3;
        N_t = 1;    M_t = 128;   deltaf_t = 15e3;

        fc       = 4e9;
        velocity = 0;                        % Dopplers set directly in the channel

    otherwise
        error('Unknown param_mode "%s" (use ''eva'', ''paper'', ''paperv2'' or ''paperdopp'')', param_mode);
end

T_o = 1/deltaf_o;   Ts_o = T_o;             % Tcp = 0 throughout
T_t = 1/deltaf_t;   Ts_t = T_t;
sigma2   = 1;
nu_max   = velocity * fc / c;

% ====== Channel definitions (selected by channel_type) ===================
% The channel (tau, nu, alpha) is defined once and reused for OFDM and CP-OTFS.
switch lower(channel_type)
    case 'eva'
        time_delays      = [0, 150, 370] * 1e-9;
        relative_powers  = [0, -1.4, -0.6];
        rp_lin           = 10.^(relative_powers/10);
        path_amp         = sqrt(rp_lin / sum(rp_lin));

        P = numel(time_delays);

        rng(42);
        theta        = 2*pi*rand(P,1) - pi;
        nu_vec       = nu_max * cos(theta);
        if nu_max > 0
            nu_vec(1) = nu_max;
        end
        tau_vec   = time_delays(:);
        alpha_vec = path_amp(:);

    case 'paper'
        % 40-path delay-only channel.  Delays uniform on
        %   [0, (N_o - 1) * 1/B]   (= paper's "tau in [0, N-1] samples"
        % with paper's N = our N_o).  Capping at (N_o-1)*bin_t = T - bin_t
        % keeps delays inside one OFDM symbol duration T, avoiding the
        % per-symbol cyclic-shift aliasing that would otherwise collapse
        % paths with tau differing by integer multiples of T.
        L_paths      = 40;
        bin_t_ref    = 1/(N_o*deltaf_o);              % 1/B
        tau_max      = (N_o - 1) * bin_t_ref;         % within one OFDM symbol
        min_sep_samp = 1.0;
        rng(20260514);
        usable = tau_max - (L_paths-1)*min_sep_samp*bin_t_ref;
        if usable <= 0
            error('tau_max too small for %d paths with min separation 1 sample', ...
                  L_paths);
        end
        u            = sort(usable * rand(L_paths, 1));
        tau_vec      = u + min_sep_samp*bin_t_ref*(0:L_paths-1).';
        nu_vec       = zeros(L_paths, 1);
        alpha_vec    = ones(L_paths, 1);
        P            = L_paths;

    case 'paperv2'
        % Same 40-path delay profile as 'paper' (identical rng draw), but
        % each path additionally carries a Doppler shift drawn uniformly
        % in [-nu_max, nu_max].  Use with param_mode = 'paperv2' so that
        % nu_max > 0 and the frame spans several time slots; with a
        % zero-Doppler param mode this degenerates to the 'paper' channel.
        L_paths      = 40;
        bin_t_ref    = 1/(N_o*deltaf_o);              % 1/B
        tau_max      = (N_o - 1) * bin_t_ref;         % within one OFDM symbol
        min_sep_samp = 1.0;
        rng(20260514);                                % same delays as 'paper'
        usable = tau_max - (L_paths-1)*min_sep_samp*bin_t_ref;
        if usable <= 0
            error('tau_max too small for %d paths with min separation 1 sample', ...
                  L_paths);
        end
        u            = sort(usable * rand(L_paths, 1));
        tau_vec      = u + min_sep_samp*bin_t_ref*(0:L_paths-1).';
        rng(20260707);                                % independent Doppler draw
        nu_vec       = nu_max * (2*rand(L_paths, 1) - 1);
        alpha_vec    = ones(L_paths, 1);
        P            = L_paths;

    case 'paperdopp'
        % 40-path DOPPLER-ONLY channel: tau = 0 for every path; Dopplers uniform
        % over the observable range [0, (M-1)*bin_nu] with 1-bin min separation,
        % bin_nu = 1/(M*T) = deltaf/M.  The dual of the 'paper' delay-only channel
        % (same construction, delay<->Doppler).  Requires param_mode='paperdopp'
        % (M>1) so Doppler is observable.
        L_paths   = 40;
        bin_nu    = deltaf_o / M_o;                    % Doppler bin 1/(M*T)
        nu_span   = (M_o - 1) * bin_nu;                % < slot-Nyquist 1/T = deltaf
        min_sep   = 1.0;                               % bins
        rng(20260514);                                 % same draw pattern as 'paper'
        usable = nu_span - (L_paths-1)*min_sep*bin_nu;
        if usable <= 0
            error('too many Doppler paths (%d) for M=%d bins', L_paths, M_o);
        end
        u         = sort(usable * rand(L_paths, 1));
        nu_vec    = u + min_sep*bin_nu*(0:L_paths-1).';
        tau_vec   = zeros(L_paths, 1);                 % NO delay
        alpha_vec = ones(L_paths, 1);
        P         = L_paths;

    otherwise
        error('Unknown channel_type "%s" (use ''EVA'', ''paper'', ''paperv2'' or ''paperdopp'')', channel_type);
end

% Estimate a parameter iff the channel actually varies it: delay iff any tau~=0,
% Doppler iff any nu~=0.  paper -> delay-only; paperdopp -> Doppler-only;
% eva/paperv2 -> both.  crb_multi selects the estimated (tau,nu) sub-block.
estimate_doppler = any(nu_vec ~= 0);
estimate_delay   = any(tau_vec ~= 0);
assert(estimate_delay || estimate_doppler, 'channel has neither delay nor Doppler to estimate');

% ----- Echo the active configuration ------------------------------------
fprintf('param_mode  = %s\n', param_mode);
fprintf('channel     = %s   (P = %d paths)\n', upper(channel_type), P);
fprintf('modulation  = %s\n', modulation);
fprintf('OFDM    grid : N=%d, M=%d, deltaf=%.3g Hz, T=%.3g s\n', N_o, M_o, deltaf_o, T_o);
fprintf('CP-OTFS grid : N=%d, M=%d, deltaf=%.3g Hz, T=%.3g s\n', N_t, M_t, deltaf_t, T_t);
fprintf('Other     : fc=%.3g Hz, v=%.2f m/s, nu_max=%.3g Hz, sigma2=%.2f\n', ...
        fc, velocity, nu_max, sigma2);
fprintf('estimate_delay   = %d   (tau %s an unknown FIM parameter)\n', ...
        estimate_delay, ternary(estimate_delay, 'IS', 'is NOT'));
fprintf('estimate_doppler = %d   (nu %s an unknown FIM parameter)\n', ...
        estimate_doppler, ternary(estimate_doppler, 'IS', 'is NOT'));
if P <= 10
    for p = 1:P
        fprintf('  path %2d:  tau = %7.2f ns,  nu = %+7.2f Hz,  |alpha|^2 = %.3f\n', ...
                p, tau_vec(p)*1e9, nu_vec(p), abs(alpha_vec(p))^2);
    end
else
    fprintf('  tau range : [%.1f, %.1f] ns\n', min(tau_vec)*1e9, max(tau_vec)*1e9);
    fprintf('  nu  range : [%+g, %+g] Hz\n', min(nu_vec), max(nu_vec));
    fprintf('  |alpha|^2 : min=%.3f, max=%.3f, mean=%.3f\n', ...
            min(abs(alpha_vec).^2), max(abs(alpha_vec).^2), mean(abs(alpha_vec).^2));
end
fprintf('\n');

% ----- Per-realization unitarity sanity check ----------------------------
% Only meaningful when OFDM and CP-OTFS share the same grid.
if (N_o == N_t) && (M_o == M_t) && (deltaf_o == deltaf_t)
    rng(0);
    X0  = gen_symbols(N_o, M_o, modulation);
    x0  = sfft_u(X0);
    J_OF_ICI = getFIM_OFDM_multi_ICI(X0, deltaf_o, Ts_o, sigma2, alpha_vec, tau_vec, nu_vec);
    J_CPOT_ICI = getFIM_CP_OTFS_multi_ICI(x0, deltaf_t, Ts_t, sigma2, alpha_vec, tau_vec, nu_vec);
    rel = norm(J_OF_ICI - J_CPOT_ICI,'fro') / norm(J_OF_ICI,'fro');
    fprintf('Per-realization (x = SFFT(X)) ICI FIM rel. error = %.3e\n', rel);

    % RCP-OTFS is NOT (I)SFFT-equivalent to OFDM: the single frame CP adds the
    % intra-frame Doppler phase and the boundary-wrap phase.  Report the
    % per-realization RCP-vs-CP departure as a magnitude check on those terms.
    J_RCPOT_ICI = getFIM_RCP_OTFS_multi_ICI(x0, deltaf_t, Ts_t, sigma2, alpha_vec, tau_vec, nu_vec);
    rel_rcp = norm(J_RCPOT_ICI - J_CPOT_ICI,'fro') / norm(J_CPOT_ICI,'fro');
    fprintf('Per-realization RCP-vs-CP ICI FIM rel. departure = %.3e (single-CP phases)\n', rel_rcp);

    % RCP-OFDM vs RCP-OTFS: same single-frame-CP time signal and unitary
    % receivers, so the FIMs should coincide; the residual measures the
    % discrete-delay-model approximations of the RCP-OTFS derivation
    % (continuous correlator vs sampled receiver, see rcp_ofdm_fim.tex).
    J_RCPOF_ICI0 = getFIM_RCP_OFDM_multi_ICI(X0, deltaf_o, Ts_o, sigma2, alpha_vec, tau_vec, nu_vec);
    rel_rcpof = norm(J_RCPOT_ICI - J_RCPOF_ICI0,'fro') / norm(J_RCPOF_ICI0,'fro');
    fprintf('Per-realization RCP-OTFS vs RCP-OFDM ICI FIM rel. departure = %.3e (model approx.)\n\n', rel_rcpof);
else
    fprintf('OFDM and CP-OTFS use different grids; per-realization unitarity check skipped.\n\n');
end

% ----- Monte Carlo over independent QPSK draws ---------------------------
nTrials = 50;

CRBt_OF_I = zeros(nTrials, P);  CRBn_OF_I = zeros(nTrials, P);
CRBt_RCPOF_I = zeros(nTrials, P);  CRBn_RCPOF_I = zeros(nTrials, P);
CRBt_CPOT_I = zeros(nTrials, P);  CRBn_CPOT_I = zeros(nTrials, P);
CRBt_RCPOT_I = zeros(nTrials, P);  CRBn_RCPOT_I = zeros(nTrials, P);
CRBt_OF_N = zeros(nTrials, P);  CRBn_OF_N = zeros(nTrials, P);
CRBt_RCPOF_N = zeros(nTrials, P);  CRBn_RCPOF_N = zeros(nTrials, P);
CRBt_CPOT_N = zeros(nTrials, P);  CRBn_CPOT_N = zeros(nTrials, P);
CRBt_RCPOT_N = zeros(nTrials, P);  CRBn_RCPOT_N = zeros(nTrials, P);

% ----- Paper precoded family (Fig. 1 reproduction) ----------------------
% Single N=128 block, delay-only.  OFDM/SC/OTFS differ ONLY through the
% subcarrier-domain power |Q^H s|^2 = |F_N U s|^2; delay is resolved over the
% full N=128 steering inside getFIM_OFDM (decoupled from the OTFS precoder
% size), so all L=40 targets are resolvable -- unlike the physical CP/RCP
% OTFS frames whose delay bins are limited to N_t.  Precoders (paper eq. 8):
%   OFDM: U = F_N^H  -> Q^H = I ;  SC: U = I -> Q^H = F_N ;
%   OTFS: U = F_N1^H (x) I_N2 -> Q^H = F_N (F_N1^H (x) I_N2).
CRBt_POF = [];  CRBt_PSC = [];  CRBt_POT = [];  CRBt_PJEN = [];
do_precoded = strcmpi(param_mode,'paper') && (M_o == 1) && (N_o == 32*4);
if do_precoded
    N1p = 32;  N2p = 4;                        % paper OTFS Kronecker factors
    FN  = fft(eye(N_o)) / sqrt(N_o);
    FN1 = fft(eye(N1p)) / sqrt(N1p);
    QH_OFDM = eye(N_o);
    QH_SC   = FN;
    QH_OTFS = FN * kron(FN1', eye(N2p));
    CRBt_POF = zeros(nTrials, P);
    CRBt_PSC = zeros(nTrials, P);
    CRBt_POT = zeros(nTrials, P);
    % Jensen bound: data-independent, flat |Q^H s|^2 = 1.
    [CRBt_PJEN, ~] = crb_multi(getFIM_OFDM_multi_noICI(ones(N_o,1), ...
        deltaf_o, Ts_o, sigma2, alpha_vec, tau_vec, nu_vec), P, estimate_doppler, estimate_delay);
end

tic
for t = 1:nTrials
    % Independent data symbols sized to each waveform's own grid.
    X = gen_symbols(N_o, M_o, modulation);
    x = gen_symbols(N_t, M_t, modulation);

    % (CP-)OFDM FIMs use its own (N_o, M_o, deltaf_o).
    J_OF_I = getFIM_OFDM_multi_ICI  (X, deltaf_o, Ts_o, sigma2, alpha_vec, tau_vec, nu_vec);
    J_OF_N = getFIM_OFDM_multi_noICI(X, deltaf_o, Ts_o, sigma2, alpha_vec, tau_vec, nu_vec);
    % RCP-OFDM (single frame CP) on the same TF grid and data.
    J_RCPOF_I = getFIM_RCP_OFDM_multi_ICI  (X, deltaf_o, Ts_o, sigma2, alpha_vec, tau_vec, nu_vec);
    J_RCPOF_N = getFIM_RCP_OFDM_multi_noICI(X, deltaf_o, Ts_o, sigma2, alpha_vec, tau_vec, nu_vec);
    % CP-OTFS and RCP-OTFS FIMs use the OTFS grid (N_t, M_t, deltaf_t).
    J_CPOT_I = getFIM_CP_OTFS_multi_ICI  (x, deltaf_t, Ts_t, sigma2, alpha_vec, tau_vec, nu_vec);
    J_CPOT_N = getFIM_CP_OTFS_multi_noICI(x, deltaf_t, Ts_t, sigma2, alpha_vec, tau_vec, nu_vec);
    % RCP-OTFS (single frame CP): same grid, two extra deterministic phases.
    J_RCPOT_I = getFIM_RCP_OTFS_multi_ICI  (x, deltaf_t, Ts_t, sigma2, alpha_vec, tau_vec, nu_vec);
    J_RCPOT_N = getFIM_RCP_OTFS_multi_noICI(x, deltaf_t, Ts_t, sigma2, alpha_vec, tau_vec, nu_vec);

    [CRBt_OF_I(t,:), CRBn_OF_I(t,:)] = crb_multi(J_OF_I, P, estimate_doppler, estimate_delay);
    [CRBt_RCPOF_I(t,:), CRBn_RCPOF_I(t,:)] = crb_multi(J_RCPOF_I, P, estimate_doppler, estimate_delay);
    [CRBt_CPOT_I(t,:), CRBn_CPOT_I(t,:)] = crb_multi(J_CPOT_I, P, estimate_doppler, estimate_delay);
    [CRBt_RCPOT_I(t,:), CRBn_RCPOT_I(t,:)] = crb_multi(J_RCPOT_I, P, estimate_doppler, estimate_delay);
    [CRBt_OF_N(t,:), CRBn_OF_N(t,:)] = crb_multi(J_OF_N, P, estimate_doppler, estimate_delay);
    [CRBt_RCPOF_N(t,:), CRBn_RCPOF_N(t,:)] = crb_multi(J_RCPOF_N, P, estimate_doppler, estimate_delay);
    [CRBt_CPOT_N(t,:), CRBn_CPOT_N(t,:)] = crb_multi(J_CPOT_N, P, estimate_doppler, estimate_delay);
    [CRBt_RCPOT_N(t,:), CRBn_RCPOT_N(t,:)] = crb_multi(J_RCPOT_N, P, estimate_doppler, estimate_delay);

    % Paper precoded waveforms: same single-block delay-only OFDM FIM, fed the
    % subcarrier-domain symbols Q^H s for each waveform's precoder.
    if do_precoded
        s = X(:);                              % paper single-block symbols (N_o x 1)
        [CRBt_POF(t,:), ~] = crb_multi(getFIM_OFDM_multi_noICI(QH_OFDM*s, deltaf_o, Ts_o, sigma2, alpha_vec, tau_vec, nu_vec), P, estimate_doppler, estimate_delay);
        [CRBt_PSC(t,:), ~] = crb_multi(getFIM_OFDM_multi_noICI(QH_SC  *s, deltaf_o, Ts_o, sigma2, alpha_vec, tau_vec, nu_vec), P, estimate_doppler, estimate_delay);
        [CRBt_POT(t,:), ~] = crb_multi(getFIM_OFDM_multi_noICI(QH_OTFS*s, deltaf_o, Ts_o, sigma2, alpha_vec, tau_vec, nu_vec), P, estimate_doppler, estimate_delay);
    end
end
fprintf('MC time: '); toc

% ----- Per-path / trace averaged delay CRBs -----------------------------
if estimate_delay
if P <= 10
    fprintf('\n========== Per-path averaged CRBs ==========\n');
    for p = 1:P
        fprintf('--- Path %d  (tau=%7.2f ns) ---\n', p, tau_vec(p)*1e9);
        fprintf('  CP-OFDM (ICI)         CRB(tau)=%.4e\n', mean(CRBt_OF_I(:,p)));
        fprintf('  RCP-OFDM (ICI)        CRB(tau)=%.4e\n', mean(CRBt_RCPOF_I(:,p)));
        fprintf('  CP-OTFS (OFDM-ICI)    CRB(tau)=%.4e\n', mean(CRBt_CPOT_I(:,p)));
        fprintf('  RCP-OTFS (OFDM-ICI)   CRB(tau)=%.4e\n', mean(CRBt_RCPOT_I(:,p)));
        fprintf('  CP-OFDM (no-ICI)      CRB(tau)=%.4e\n', mean(CRBt_OF_N(:,p)));
        fprintf('  RCP-OFDM (no-ICI)     CRB(tau)=%.4e\n', mean(CRBt_RCPOF_N(:,p)));
        fprintf('  CP-OTFS (no-OFDM-ICI) CRB(tau)=%.4e\n', mean(CRBt_CPOT_N(:,p)));
        fprintf('  RCP-OTFS (no-OFDM-ICI)CRB(tau)=%.4e\n\n', mean(CRBt_RCPOT_N(:,p)));
    end
else
    fprintf('\n========== Trace averaged CRBs (P=%d paths) ==========\n', P);
    fprintf('  CP-OFDM (ICI)         tr CRB(tau)=%.4e\n', sum(mean(CRBt_OF_I,1)));
    fprintf('  RCP-OFDM (ICI)        tr CRB(tau)=%.4e\n', sum(mean(CRBt_RCPOF_I,1)));
    fprintf('  CP-OTFS (OFDM-ICI)    tr CRB(tau)=%.4e\n', sum(mean(CRBt_CPOT_I,1)));
    fprintf('  RCP-OTFS (OFDM-ICI)   tr CRB(tau)=%.4e\n', sum(mean(CRBt_RCPOT_I,1)));
    fprintf('  CP-OFDM (no-ICI)      tr CRB(tau)=%.4e\n', sum(mean(CRBt_OF_N,1)));
    fprintf('  RCP-OFDM (no-ICI)     tr CRB(tau)=%.4e\n', sum(mean(CRBt_RCPOF_N,1)));
    fprintf('  CP-OTFS (no-OFDM-ICI) tr CRB(tau)=%.4e\n', sum(mean(CRBt_CPOT_N,1)));
    fprintf('  RCP-OTFS (no-OFDM-ICI)tr CRB(tau)=%.4e\n', sum(mean(CRBt_RCPOT_N,1)));
end
end   % if estimate_delay
if estimate_doppler
    fprintf('\n---------- Trace averaged Doppler CRBs ----------\n');
    fprintf('  CP-OFDM (ICI)         tr CRB(nu)=%.4e\n', sum(mean(CRBn_OF_I,1)));
    fprintf('  RCP-OFDM (ICI)        tr CRB(nu)=%.4e\n', sum(mean(CRBn_RCPOF_I,1)));
    fprintf('  CP-OTFS (OFDM-ICI)    tr CRB(nu)=%.4e\n', sum(mean(CRBn_CPOT_I,1)));
    fprintf('  RCP-OTFS (OFDM-ICI)   tr CRB(nu)=%.4e\n', sum(mean(CRBn_RCPOT_I,1)));
    fprintf('  CP-OFDM (no-ICI)      tr CRB(nu)=%.4e\n', sum(mean(CRBn_OF_N,1)));
    fprintf('  RCP-OFDM (no-ICI)     tr CRB(nu)=%.4e\n', sum(mean(CRBn_RCPOF_N,1)));
    fprintf('  CP-OTFS (no-OFDM-ICI) tr CRB(nu)=%.4e\n', sum(mean(CRBn_CPOT_N,1)));
    fprintf('  RCP-OTFS (no-OFDM-ICI)tr CRB(nu)=%.4e\n', sum(mean(CRBn_RCPOT_N,1)));
end

% ----- Save -------------------------------------------------------------
% Use a representative |alpha_1|^2 (same in both waveforms) for SNR axis
% in the plot script.
save(fullfile(resdir, 'eva_fim_avg_results.mat'), ...
     'CRBt_OF_I','CRBn_OF_I','CRBt_CPOT_I','CRBn_CPOT_I','CRBt_RCPOT_I','CRBn_RCPOT_I', ...
     'CRBt_OF_N','CRBn_OF_N','CRBt_CPOT_N','CRBn_CPOT_N','CRBt_RCPOT_N','CRBn_RCPOT_N', ...
     'CRBt_RCPOF_I','CRBn_RCPOF_I','CRBt_RCPOF_N','CRBn_RCPOF_N', ...
     'N_o','M_o','deltaf_o','Ts_o', ...
     'N_t','M_t','deltaf_t','Ts_t', ...
     'sigma2', ...
     'alpha_vec','tau_vec','nu_vec','P','nTrials', ...
     'CRBt_POF','CRBt_PSC','CRBt_POT','CRBt_PJEN','do_precoded', ...
     'channel_type','param_mode','estimate_doppler','estimate_delay','modulation');

warning('on','MATLAB:nearlySingularMatrix');

% ------------------------------------------------------------------------
function S = gen_symbols(N, M, mod_type)
% i.i.d. zero-mean, unit-power, zero-pseudovariance data symbols on the
% N x M grid.  Matches the paper's Assumptions 1-3 for QPSK, 16-QAM, and
% 16-PSK constellations (cf. crb_new/run_sim.m gen_symbols).
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

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end

function [crb_tau, crb_nu] = crb_multi(J, P, estimate_doppler, estimate_delay)
% Extract per-path CRB(tau_p) and/or CRB(nu_p) from a 4P x 4P multi-target FIM
% via Schur complement w.r.t. the (alphaR, alphaI) block.  Which of (tau, nu) are
% treated as unknown is set by (estimate_delay, estimate_doppler): delay-only
% ('paper'), Doppler-only ('paperdopp'), or joint ('eva'/'paperv2').
    if nargin < 4, estimate_delay = true; end   % backward compatible (delay-only default)

    ncomp = estimate_delay + estimate_doppler;
    idx_a = zeros(1, 2*P);
    idx_k = zeros(1, ncomp*P);
    for p = 1:P
        base = 4*(p-1);
        idx_a(2*p-1:2*p) = [base+1, base+2];
        ks = [];
        if estimate_delay,   ks = [ks, base+3]; end   % tau
        if estimate_doppler, ks = [ks, base+4]; end   % nu
        idx_k((p-1)*ncomp + (1:ncomp)) = ks;          % per-path order: [tau, nu]
    end

    Jaa = J(idx_a, idx_a);
    Jak = J(idx_a, idx_k);
    Jkk = J(idx_k, idx_k);
    Jschur = Jkk - Jak' * (Jaa \ Jak);
    C      = pinv_equilibrated((Jschur + Jschur')/2);
    D      = reshape(diag(C), ncomp, P).';            % P x ncomp, one row per path

    crb_tau = nan(1, P);  crb_nu = nan(1, P);  col = 0;
    if estimate_delay,   col = col + 1;  crb_tau = D(:, col).'; end
    if estimate_doppler, col = col + 1;  crb_nu  = D(:, col).'; end
end

function C = pinv_equilibrated(J)
% pinv of a symmetric FIM after symmetric diagonal equilibration.
% pinv alone is robust to rank deficiency (e.g. paper-mode CP-OTFS with 40
% paths in only N=32 delay bins), but its single absolute tolerance fails
% when parameters have wildly different natural scales: joint (tau, nu)
% blocks mix ~1e21 (delay, 1/s^2) with ~1e-5 (Doppler, 1/Hz^2) entries and
% the whole Doppler subspace gets truncated as "numerically zero".
% Scaling J -> D*J*D with D = diag(1./sqrt(diag(J))) puts all parameters
% on a unit diagonal, so the pinv tolerance only removes genuinely
% deficient directions; genuinely null rows/cols (zero diagonal) keep
% scale 1 and are still zeroed by pinv.
    d = sqrt(abs(diag(J)));
    d(d <= 0) = 1;
    Dinv = diag(1 ./ d);
    C = Dinv * pinv(Dinv * J * Dinv) * Dinv;
end
