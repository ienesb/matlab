% mcrb_ce_ici.m
% Misspecified CRB (MCRB) at the PSEUDO-TRUE parameter for OFDM-ICI and
% CP-OTFS-ICI sensing, where the misspecification is the REALISTIC channel-
% estimation + data-demodulation error pattern from the eva CE receiver
% (ber_ce_ici.m) -- NOT an abstract BSC.  Swept over SNR (the BSC rate p in
% compute_mcrb_mc.m is replaced by the SNR-dependent detector error rate).
%
% Mechanism (same as compute_mcrb_mc): both configs are scored in the SAME
% OFDM-ICI sensing model (getAB_OFDM_core), so the ONLY thing that differs is the
% detected-symbol error pattern Xhat.  At equal symbol-error rate the two MCRBs
% coincide; since OTFS has far lower BER at a given SNR, its Xhat is closer to X
% -> smaller pseudo-true bias (BB) and score-bias floor (F) -> LOWER MCRB.
%
% Per (SNR, trial, config), monostatic (sensing noise = data noise, sigma2 = N0):
%   1. channel eta_bar=(alpha,tau,nu): fixed EVA delays, random Doppler+phase.
%   2. transmit frame X (data + pilots), detect with the realistic CE receiver
%      (OFDM: periodic-pilot LS+interp -> 1-tap MMSE; OTFS: embedded-pilot
%      threshold CE -> MP), HARD-DECIDE data cells -> Xhat (pilots/guard kept).
%   3. eta_0 = pseudoTrue_{OFDM,CP_OTFS}(X, Xhat, ...)   (warm cold-start eta_bar).
%   4. [A,B,Jt,b] = getAB_OFDM_core(Xtf, Xhat_tf, ..., sigma2=1, eta_0, eta_bar).
%   5. affine pieces at sigma2=1: G=A^-1 Jt A^-1, F=A^-1 bb^T A^-1,
%      BB=(eta_bar-eta_0)(.)^T ;  LB(N0) = N0*G + F + BB.
%      Genie CRB reference (Xhat=X): CRB(N0) = N0 * J^-1 (b=0,F=0,BB=0).
%
% Consistency caveat: OTFS DETECTION uses the exact delay-Doppler channel
% (needed for MP), while the MCRB is scored in the ofdm_mean model -- both configs
% share that sensing model, so the comparison is fair; only the error rate differs.
clear; clc;
thisdir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisdir, 'utils'));
warning('off','MATLAB:nearlySingularMatrix');

%% ---- toggles / sweep --------------------------------------------------------
use_full_eva = false;   % false = reduced EVA (P=3, well-conditioned A); true = full EVA (P=9)
doppler_scale = 1;      % 1 = nominal EVA Doppler ~444 Hz; 0 = NO Doppler.  Sets nu_max below.
use_ici      = false;   % false = ASSUME NO ICI: no-ICI sensing model (getAB/pseudoTrue) AND
                        %         no-ICI OFDM channel (ofdm_mean).  Doppler stays (inter-symbol
                        %         phase only).  OTFS keeps its native DD channel (Doppler in DD
                        %         is not "ICI"); only its MCRB scoring is no-ICI.  true = full ICI.
ofdm_only    = false;   % true = skip the CP-OTFS config
SNR_dbs   = 0:1:25;     % sensing = data SNR sweep [dB]
nTrials   = 1000;       % MC trials per (config, SNR)   [channel+data+noise]
seed_base = 9000;       % RNG seed base; change per run to get independent trials for merging
% Auto-name so this run does not overwrite the nominal / no-Doppler files.
dtag = ''; if doppler_scale==0, dtag='_nodopp'; elseif doppler_scale~=1, dtag=sprintf('_dopp%gx',doppler_scale); end
if ~use_ici, dtag = [dtag '_noici']; end
results_name = sprintf('mcrb_ce_ici%s_%dt_results.mat', dtag, nTrials);
n_iter_mp = 30;  damping = 0.3;  E_dopp = 10;  W_ici = 8;

% ---- OFDM channel estimation (see ber_ce_channels.m / diag_ofdm_floor.m) ----
% 'interp' : the ORIGINAL LS + linear interpolation.  It floors OFDM's detection
%            BER at ~1e-2 because the channel rotates ~78 deg between frequency
%            pilots and linear interp chords that arc -- a receiver artifact, not
%            physics, which then inflates the pseudo-true bias BB.
% 'dft'    : LS fit of the (l_tau+1) delay taps across the pilot comb.
% 'dftdd'  : 'dft' + decision-directed refinement (re-fit the taps from ALL N
%            subcarriers per symbol).  Measured to REACH the genie-CE bound, so
%            OFDM's Xhat is no longer handicapped by estimation error.
ce_mode = 'dftdd';
n_dd    = 2;            % decision-directed refinement passes

% ---- BSC-matched control arm ------------------------------------------------
% Isolates the ERROR PATTERN from the error RATE.  Per trial, take the realistic
% detector's error COUNT and re-scatter that many symbol errors i.i.d. uniform
% over the data cells (each to a random different symbol) -- the abstract-BSC
% model of compute_mcrb_mc.  Same channel/data/SNR/error-count as the real arm;
% only the PLACEMENT differs (real errors cluster in deep fades / MP failures,
% BSC errors do not).  If the two MCRBs match, the BSC abstraction is justified;
% if not, clustering matters.  Costs a 2nd pseudo-true solve per trial (~2x).
do_bsc = true;

% ---- pilot budget: data-referenced and matched across waveforms -------------
% Total pilot energy is identical for OFDM and OTFS, and the non-data (guard+
% pilot) cell count is matched by match_pilot_spacing -> same data:overhead
% ratio.  Data-referenced (not the old SNRp_dB=50 noise-referenced convention,
% which made the pilots shrink as 1/SNR while data stayed unit-power and let
% data leakage swamp the OFDM pilots at high SNR).
pilot_boost_dB = 0;     % total pilot energy vs the same #cells of unit-power data

%% ---- system / grid (same as ber_ce_ici.m) ----------------------------------
N = 512;  M = 128;  deltaf = 15e3;  Ts = 1/deltaf;
fc = 4e9;  c = 3e8;  v = 120/3.6;  nu_max = v*fc/c * doppler_scale;   % 0 -> no Doppler, no ICI
M_order = 4;  constellation = qammod((0:M_order-1).', M_order, UnitAveragePower=true);
nbits   = log2(M_order);

if use_full_eva
    eva_delays_ns = [0 30 150 310 370 710 1090 1730 2510];
    eva_pow_db    = [0 -1.5 -1.4 -3.6 -0.6 -9.1 -7.0 -12.0 -16.9];
else
    eva_delays_ns = [0 150 370];
    eva_pow_db    = [0 -1.4 -0.6];
end
evaname    = 'full EVA';  if ~use_full_eva, evaname = 'reduced EVA'; end
P          = numel(eva_delays_ns);
delay_taps = round(eva_delays_ns*1e-9 * N * deltaf);
tau        = (delay_taps / (N*deltaf)).';                 % P x 1 rounded delays [s]
amp        = (10.^(eva_pow_db/20)).';  amp = amp / norm(amp);   % P x 1, sum|a|^2=1
l_tau      = ceil(max(eva_delays_ns*1e-9) * N * deltaf);
k_nu       = ceil(nu_max * M / deltaf);

nE = 4*P;  idx_tau = 4*(0:P-1) + 3;  idx_nu = 4*(0:P-1) + 4;
bin_tau = 1/(N*deltaf);   bin_nu = deltaf/M;              % resolution bins (for reference)

%% ---- OTFS embedded pilot + full guard --------------------------------------
l_p = l_tau;  k_p = M/2;
guard_rows = 1 : (2*l_tau + 1);
ce_rows    = (l_p+1) : (2*l_tau + 1);
data_mask_otfs = true(N, M);  data_mask_otfs(guard_rows, :) = false;

%% ---- OFDM periodic pilots: overhead-matched to OTFS guard+pilot -------------
otfs_nondata = numel(guard_rows) * M;
[dN_p, dM_p, n_pilot_cells] = match_pilot_spacing(N, M, otfs_nondata, l_tau, k_nu);
n_pilots = 1:dN_p:N;  m_pilots = 1:dM_p:M;
pilot_mask_ofdm = false(N, M);  pilot_mask_ofdm(n_pilots, m_pilots) = true;
data_mask_ofdm  = ~pilot_mask_ofdm;

% ---- matched, data-referenced pilot budget (independent of N0) --------------
E_pil = n_pilot_cells * db2pow(pilot_boost_dB);   % total pilot energy, BOTH waveforms
a_p   = sqrt(E_pil / n_pilot_cells);              % OFDM: spread over the lattice
x_p   = sqrt(E_pil);                              % OTFS: concentrated in one cell
assert(abs(n_pilot_cells*a_p^2 - x_p^2) < 1e-9*E_pil, 'pilot energy mismatch');
assert(otfs_nondata == n_pilot_cells, 'non-data overhead mismatch');
fprintf(['FAIRNESS: non-data cells OTFS=%d == OFDM=%d (%.2f%%) | ' ...
         'total pilot E OTFS=%.4e == OFDM=%.4e\n'], ...
        otfs_nondata, n_pilot_cells, 100*n_pilot_cells/(N*M), x_p^2, n_pilot_cells*a_p^2);

configs = { 'OFDM ICI', 'ofdm'; 'CP-OTFS ICI', 'otfs' };
if ofdm_only, configs = configs(1,:); end
nC = size(configs, 1);
nSNR = numel(SNR_dbs);

% Per-config SNR cap: OFDM to 25 dB, CP-OTFS to 15 dB (with no Doppler the OTFS
% detector is error-free past ~15 dB so BB->0 there anyway).  Capped (uncomputed)
% blocks are left as NaN so plots/analysis skip them.
snr_cap = 25*ones(1,nC);
for cc = 1:nC, if strcmp(configs{cc,2},'otfs'), snr_cap(cc) = 15; end, end

%% ---- storage (trial-averaged affine pieces at sigma2=1, per SNR & config) ---
G_bar  = zeros(nE, nE, nSNR, nC);   F_bar  = zeros(nE, nE, nSNR, nC);
BB_bar = zeros(nE, nE, nSNR, nC);   Gc_bar = zeros(nE, nE, nSNR, nC);
% BSC-matched arm (same error count, i.i.d. placement):
G_bar_b  = zeros(nE, nE, nSNR, nC);  F_bar_b = zeros(nE, nE, nSNR, nC);
BB_bar_b = zeros(nE, nE, nSNR, nC);  ser_data_b = zeros(nSNR, nC);
kappa_max = zeros(nSNR, nC);
ser_data  = zeros(nSNR, nC);        % detector symbol-error rate (drives the mismatch)
ber_data  = zeros(nSNR, nC);        % detector BIT-error rate (true BER, Gray 4-QAM)

fprintf('MCRB at pseudo-true from REALISTIC CE detection (%s, N=%d M=%d P=%d)\n', evaname, N, M, P);
fprintf('SNR in [%s] dB, pilot boost %+d dB (data-referenced), %d trials, OFDM CE = %s\n', ...
        num2str(SNR_dbs), pilot_boost_dB, nTrials, ce_mode);
tAll = tic;

for ci = 1:nC
    modn = configs{ci,2};
    for si = 1:nSNR
        if SNR_dbs(si) > snr_cap(ci)             % skip capped blocks -> leave NaN
            G_bar(:,:,si,ci)=NaN; F_bar(:,:,si,ci)=NaN; BB_bar(:,:,si,ci)=NaN; Gc_bar(:,:,si,ci)=NaN;
            G_bar_b(:,:,si,ci)=NaN; F_bar_b(:,:,si,ci)=NaN; BB_bar_b(:,:,si,ci)=NaN;
            ser_data(si,ci)=NaN; ber_data(si,ci)=NaN; ser_data_b(si,ci)=NaN; kappa_max(si,ci)=NaN;
            continue;
        end
        N0  = 1/db2pow(SNR_dbs(si));   % a_p / x_p are data-referenced -> N0-independent
        rng(7000 + si);

        G_t  = zeros(nE, nE, nTrials);  F_t  = zeros(nE, nE, nTrials);
        BB_t = zeros(nE, nE, nTrials);  Gc_t = zeros(nE, nE, nTrials);
        ka_t = zeros(1, nTrials);       be_t = zeros(1, nTrials);   % SER
        bb_t = zeros(1, nTrials);                                   % BER (bits)
        Gb_t = zeros(nE, nE, nTrials);  Fb_t = zeros(nE, nE, nTrials);  % BSC arm
        BBb_t = zeros(nE, nE, nTrials); beb_t = zeros(1, nTrials);

        parfor t = 1:nTrials
            warning('off','MATLAB:nearlySingularMatrix');
            % Independent, reproducible per-(SNR,trial) stream (counter-based Threefry).
            % seed_base distinguishes runs so separate runs can be merged as new trials.
            rs = RandStream('Threefry', 'Seed', seed_base);
            rs.Substream = 10000*si + t;                % 10000 spacing -> up to 10000 trials/point
            RandStream.setGlobalStream(rs);
            % ---- channel realization (true parameter) ----
            nu      = (nu_max * cos(2*pi*rand(1,P) - pi)).';
            alpha   = amp .* exp(1j*2*pi*rand(P,1));
            eta_bar = packEta(alpha, tau, nu);

            % ---- transmit frame + realistic detection -> (d, dh) native domain
            sym = randi(M_order, N, M) - 1;
            Xd  = qammod(sym, M_order, UnitAveragePower=true);
            Z   = sqrt(N0/2) * (randn(N,M) + 1j*randn(N,M));

            if strcmp(modn, 'ofdm')
                d = Xd;  d(pilot_mask_ofdm) = a_p;                 % TF frame + pilots
                Y = ofdm_mean(d, deltaf, Ts, alpha, tau, nu, use_ici) + Z;
                % CE: 'dftdd' reaches the genie-CE bound (ber_ce_channels.m), so
                % Xhat is limited by ICI/diversity, not by estimation error.
                H_hat  = ce_ofdm(Y, d, pilot_mask_ofdm, n_pilots, m_pilots, ...
                                 N, M, l_tau, N0, M_order, ce_mode, n_dd);
                Xsoft  = conj(H_hat) .* Y ./ (abs(H_hat).^2 + N0);   % 1-tap MMSE
                dh = d;                                              % pilots preserved
                bh = qamdemod(Xsoft(data_mask_ofdm), M_order, OutputType="bit", UnitAveragePower=true);
                dh(data_mask_ofdm) = qammod(bh, M_order, InputType="bit", UnitAveragePower=true);
                errmask = data_mask_ofdm;
                Xtf = d;  Xhtf = dh;
                eta0 = pseudoTrue_OFDM(d, dh, deltaf, Ts, alpha, tau, nu, use_ici, [], true);
                if do_bsc
                    dh_b     = bsc_match(d, dh, sym, errmask, M_order);
                    Xhtf_b   = dh_b;
                    eta0_b   = pseudoTrue_OFDM(d, dh_b, deltaf, Ts, alpha, tau, nu, use_ici, [], true);
                end
            else
                d = zeros(N,M);  d(data_mask_otfs) = Xd(data_mask_otfs);
                d(l_p+1, k_p+1) = x_p;                              % DD frame + pilot
                Gtrue = build_effective_dd_eva(alpha, delay_taps, nu, N, M, deltaf, Inf);
                Y = apply_dd_channel_eva(d, Gtrue) + Z;
                h_ce = zeros(N,M);
                h_ce(ce_rows,:) = Y(ce_rows,:) / x_p;
                h_ce(abs(h_ce) < 3*sqrt(N0)/x_p) = 0;
                [rc,cc,gc] = find(h_ce);
                dh = d;                                             % start from truth (pilot/guard)
                if ~isempty(gc)
                    dl_c = mod(rc-(l_p+1), N);  dk_c = mod(cc-(k_p+1), M);
                    known = zeros(N,M);  known(l_p+1,k_p+1) = x_p;
                    xhat = mp_detector(Y, dl_c, dk_c, gc, N0, constellation, ...
                                       n_iter_mp, data_mask_otfs, known, damping);
                    dh(data_mask_otfs) = xhat(data_mask_otfs);      % detected data only
                end
                errmask = data_mask_otfs;
                Xtf = isfft_u(d);  Xhtf = isfft_u(dh);
                eta0 = pseudoTrue_CP_OTFS(d, dh, deltaf, Ts, alpha, tau, nu, use_ici, [], true);
                if do_bsc
                    dh_b     = bsc_match(d, dh, sym, errmask, M_order);
                    Xhtf_b   = isfft_u(dh_b);
                    eta0_b   = pseudoTrue_CP_OTFS(d, dh_b, deltaf, Ts, alpha, tau, nu, use_ici, [], true);
                end
            end
            be_t(t) = mean(dh(errmask) ~= d(errmask));              % symbol-error rate
            % True BER from the SAME detected symbols (dh holds the detected
            % constellation points on the data cells for both waveforms).
            sh_idx  = qamdemod(dh(errmask), M_order, UnitAveragePower=true);
            bb_t(t) = mean(int2bit(sym(errmask),nbits) ~= int2bit(sh_idx,nbits), 'all');

            % ---- MCRB affine pieces at sigma2 = 1 (evaluated at pseudo-true) ----
            [a0, t0, n0] = unpackEta(eta0);
            [A, ~, Jt, b] = getAB_OFDM_core(Xtf, Xhtf, deltaf, Ts, 1, a0, t0, n0, use_ici, alpha, tau, nu);
            [G, kap] = stab_quad(A, Jt);
            F  = stab_quad(A, b*b.');
            bias = eta_bar - eta0;   BB = bias*bias.';
            % ---- genie CRB (Xhat = X): b=0, F=0, BB=0, G = J^-1 ----
            [Ac, ~, Jtc, ~] = getAB_OFDM_core(Xtf, Xtf, deltaf, Ts, 1, alpha, tau, nu, use_ici, alpha, tau, nu);
            Gc = stab_quad(Ac, Jtc);

            % ---- BSC-matched arm: same error count, i.i.d. placement ----
            if do_bsc
                [ab0,tb0,nb0] = unpackEta(eta0_b);
                [Ab,~,Jtb,bvb] = getAB_OFDM_core(Xtf, Xhtf_b, deltaf, Ts, 1, ab0, tb0, nb0, use_ici, alpha, tau, nu);
                Gb = stab_quad(Ab, Jtb);
                Fb = stab_quad(Ab, bvb*bvb.');
                biasb = eta_bar - eta0_b;   BBb = biasb*biasb.';
                beb_t(t) = mean(dh_b(errmask) ~= d(errmask));
            else
                Gb = zeros(nE); Fb = zeros(nE); BBb = zeros(nE); beb_t(t) = NaN;
            end

            G_t(:,:,t)=G;  F_t(:,:,t)=F;  BB_t(:,:,t)=BB;  Gc_t(:,:,t)=Gc;  ka_t(t)=kap;
            Gb_t(:,:,t)=Gb;  Fb_t(:,:,t)=Fb;  BBb_t(:,:,t)=BBb;
        end

        G_bar(:,:,si,ci)=mean(G_t,3);  F_bar(:,:,si,ci)=mean(F_t,3);
        BB_bar(:,:,si,ci)=mean(BB_t,3);  Gc_bar(:,:,si,ci)=mean(Gc_t,3);
        kappa_max(si,ci)=max(ka_t);  ser_data(si,ci)=mean(be_t);  ber_data(si,ci)=mean(bb_t);
        G_bar_b(:,:,si,ci)=mean(Gb_t,3);  F_bar_b(:,:,si,ci)=mean(Fb_t,3);
        BB_bar_b(:,:,si,ci)=mean(BBb_t,3);  ser_data_b(si,ci)=mean(beb_t);

        N0 = 1/db2pow(SNR_dbs(si));
        LB = N0*G_bar(:,:,si,ci) + F_bar(:,:,si,ci) + BB_bar(:,:,si,ci);
        rt = sqrt(sum(real(diag(LB(idx_tau,idx_tau)))));
        % Split the delay bound into its pieces so BB (the pseudo-true parameter
        % bias) is visible next to the total: it dominates once SER > 0.
        tG = N0*sum(real(diag(G_bar (idx_tau,idx_tau,si,ci))));
        tF =    sum(real(diag(F_bar (idx_tau,idx_tau,si,ci))));
        tB =    sum(real(diag(BB_bar(idx_tau,idx_tau,si,ci))));
        rc = sqrt(N0*sum(real(diag(Gc_bar(idx_tau,idx_tau,si,ci)))));   % genie CRB
        % BSC-matched delay RMSE (same SER, i.i.d. placement) for comparison.
        LBb = N0*G_bar_b(:,:,si,ci) + F_bar_b(:,:,si,ci) + BB_bar_b(:,:,si,ci);
        rtb = sqrt(sum(real(diag(LBb(idx_tau,idx_tau)))));
        fprintf(['  %-12s SNR=%2d dB | SER=%.2e (bsc %.2e) BER=%.2e | RMSE_tau: real=%.3e ' ...
                 'BSC=%.3e CRB=%.3e s | BB %.0f%% | cond<=%.1e (%.0fs)\n'], ...
                configs{ci,1}, SNR_dbs(si), ser_data(si,ci), ser_data_b(si,ci), ber_data(si,ci), ...
                rt, rtb, rc, 100*tB/max(tG+tF+tB,realmin), kappa_max(si,ci), toc(tAll));

        % --- incremental checkpoint (raw pieces) so a long run is crash-safe ---
        save(fullfile(thisdir,'results',['CKPT_' results_name]), ...
             'SNR_dbs','snr_cap','configs','nTrials','N','M','P','idx_tau','idx_nu', ...
             'G_bar','F_bar','BB_bar','Gc_bar','G_bar_b','F_bar_b','BB_bar_b', ...
             'ser_data','ber_data','ser_data_b','kappa_max','ci','si');
    end
end

%% ---- assemble curves (physical RMSE), MCRB and BB kept SEPARATE -------------
% Convention: MCRB = A^-1 B A^-1 = N0*G + F ;  BB = (eta_bar-eta_0)(eta_bar-eta_0)^T ;
%             LB = MCRB + BB.  Capped (NaN) blocks stay NaN.
rmse_tau_mcrb = nan(nSNR,nC);  rmse_tau_bb = nan(nSNR,nC);  rmse_tau_lb = nan(nSNR,nC);  rmse_tau_crb = nan(nSNR,nC);
rmse_nu_mcrb  = nan(nSNR,nC);  rmse_nu_bb  = nan(nSNR,nC);  rmse_nu_lb  = nan(nSNR,nC);  rmse_nu_crb  = nan(nSNR,nC);
trd = @(Xm,idx) sqrt(sum(real(diag(Xm(idx,idx)))));
for ci = 1:nC
    for si = 1:nSNR
        N0 = 1/db2pow(SNR_dbs(si));
        MCRB = N0*G_bar(:,:,si,ci) + F_bar(:,:,si,ci);      % A^-1 B A^-1
        BB   = BB_bar(:,:,si,ci);                            % pseudo-true bias
        LB   = MCRB + BB;                                    % total bound
        CRB  = N0*Gc_bar(:,:,si,ci);                         % genie (no mismatch)
        rmse_tau_mcrb(si,ci)=trd(MCRB,idx_tau); rmse_tau_bb(si,ci)=trd(BB,idx_tau);
        rmse_tau_lb(si,ci)  =trd(LB,idx_tau);   rmse_tau_crb(si,ci)=trd(CRB,idx_tau);
        rmse_nu_mcrb(si,ci) =trd(MCRB,idx_nu);  rmse_nu_bb(si,ci) =trd(BB,idx_nu);
        rmse_nu_lb(si,ci)   =trd(LB,idx_nu);    rmse_nu_crb(si,ci)=trd(CRB,idx_nu);
    end
end

outpath = fullfile(thisdir, 'results', results_name);
save(outpath, ...
     'SNR_dbs','snr_cap','doppler_scale','use_ici','pilot_boost_dB','ce_mode','configs','use_full_eva','evaname','N','M','P','nTrials', ...
     'G_bar','F_bar','BB_bar','Gc_bar','G_bar_b','F_bar_b','BB_bar_b','ser_data_b','do_bsc', ...
     'kappa_max','ser_data','ber_data','idx_tau','idx_nu', ...
     'rmse_tau_mcrb','rmse_tau_bb','rmse_tau_lb','rmse_tau_crb', ...
     'rmse_nu_mcrb','rmse_nu_bb','rmse_nu_lb','rmse_nu_crb','bin_tau','bin_nu');

if doppler_scale==1 && use_ici && ~ofdm_only, plot_mcrb_ce_ici(outpath); end   % avoid clobbering nominal figures
fprintf('Done in %.0fs. Saved %s\n', toc(tAll), results_name);
warning('on','MATLAB:nearlySingularMatrix');

%% ================= local functions =================
function dh_b = bsc_match(d, dh, sym, errmask, M_order)
% Re-scatter the SAME number of symbol errors the real detector made, i.i.d.
% uniform over the data cells, each to a uniformly-random DIFFERENT symbol
% (the abstract-BSC model).  Pilots/guard are kept exactly (= d), so only the
% data-cell error PLACEMENT and wrong-symbol choice are randomized; the count
% (hence SER) matches the real arm per trial.
    n_err = nnz(dh(errmask) ~= d(errmask));
    dh_b  = d;
    if n_err > 0
        di = find(errmask);
        pe = di(randperm(numel(di), n_err));                       % random positions
        w  = mod(sym(pe) + randi(M_order-1, n_err, 1), M_order);   % != true symbol
        dh_b(pe) = qammod(w, M_order, UnitAveragePower=true);
    end
end

function H_hat = ce_ofdm(Y, d, pilot_mask, n_pilots, m_pilots, N, M, l_tau, ...
                         N0, M_order, ce_mode, n_dd)
% OFDM channel estimation from the periodic pilot lattice (see ber_ce_channels.m).
%   'interp' : LS at pilots -> linear interpolation (original; chords the ~78 deg
%              phase rotation between frequency pilots -> floors OFDM at ~1e-2).
%   'dft'    : LS fit of the (l_tau+1) delay taps across the comb -> exact
%              reconstruction at every subcarrier (valid: l_tau < N/dN_p).
%   'dftdd'  : + decision-directed refinement, re-fitting the taps from ALL N
%              subcarriers per symbol (pilots known + data hard-decisions).
%              Measured to reach the genie-CE bound.
    H_ls = Y(n_pilots, m_pilots) ./ d(n_pilots, m_pilots);
    [Nq, Mq] = ndgrid(1:N, 1:M);
    Afull = exp(-1j*2*pi * (0:N-1).' * (0:l_tau) / N);
    switch ce_mode
        case 'interp'
            F_H   = griddedInterpolant({n_pilots, m_pilots}, H_ls, 'linear','nearest');
            H_hat = F_H(Nq, Mq);
            return
        case {'dft','dftdd'}
            Abas  = exp(-1j*2*pi * (n_pilots(:)-1) * (0:l_tau) / N);
            Hp    = Afull * (Abas \ H_ls);
            F_H   = griddedInterpolant({(1:N).', m_pilots}, Hp, 'linear','nearest');
            H_hat = F_H(Nq, Mq);
        otherwise
            error('Unknown ce_mode "%s"', ce_mode);
    end
    if strcmp(ce_mode,'dftdd')
        for it = 1:n_dd
            Xs   = conj(H_hat).*Y ./ (abs(H_hat).^2 + N0);
            dhat = qammod(qamdemod(Xs, M_order, UnitAveragePower=true), ...
                          M_order, UnitAveragePower=true);
            dhat(pilot_mask) = d(pilot_mask);
            H_hat = Afull * (Afull \ (Y ./ dhat));
        end
    end
end

function eta = packEta(alpha, tau, nu)
    eta = reshape([real(alpha(:)).'; imag(alpha(:)).'; tau(:).'; nu(:).'], [], 1);
end

function [Y, kappa] = stab_quad(A, C)
% Y = A^-1 C A^-1 with a Jacobi (magnitude) preconditioner, robust to a rank-
% deficient A via the pseudoinverse (copied from compute_mcrb_mc.m).
    d  = 1 ./ sqrt(abs(diag(A)));
    d(~isfinite(d)) = 1;
    D  = diag(d);
    As = D*A*D;  Cs = D*C*D;
    As = (As + As.')/2;
    Ai = pinv(As);
    Y  = real(D*(Ai*Cs*Ai)*D);
    if nargout > 1, kappa = cond(As); end
end
