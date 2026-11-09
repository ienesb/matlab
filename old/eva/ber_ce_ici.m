% ber_ce_ici.m
% BER vs SNR_d for the EVA channel (ICI, fractional Doppler, full guard) with
% REALISTIC channel estimation + data demodulation -- the receiver stage that
% will later feed detected symbols Xhat into the misspecified CRB (compute_mcrb_*).
%
% This reproduces the qualitative message of Fig. 9 of Raviteja et al. (Embedded
% Pilot-Aided CE for OTFS, IEEE TVT 2019): OTFS's delay-Doppler diversity gives a
% steep BER waterfall, while OFDM stays diversity-1 (shallow).  Curves:
%
%   OFDM (ideal)  : known ICI channel, banded MMSE equalisation (paper reference).
%   OFDM (CE)     : periodic TF pilots -> LS + 2D interp -> 1-tap MMSE.
%   OTFS (ideal)  : known DD taps + MP detector.
%   OTFS (CE)     : embedded single pilot + full guard, threshold CE + MP,
%                   at SNR_p in SNRp_dB.
%
% Pilot-power matching (as requested): the ONE OTFS pilot carries the same total
% energy as ALL the OFDM pilot cells combined.  With E_p = |x_p|^2 = SNRp_lin*N0,
% the OFDM per-pilot amplitude is sqrt(E_p / n_pilots).
%
% Self-contained in eva/: the parametric ICI channel (ofdm_mean.m), the MP
% detector (mp_detector.m), the (I)SFFT (isfft_u/sfft_u) and the ICI kernel
% (get_ici_coeffs) all live under eva/utils.  The OFDM banded MMSE is a local
% function below.  Only ICI configs are simulated.
clear; clc;
thisdir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisdir, 'utils'));

%% ---- System / grid (delay = rows N, Doppler = cols M; paper Fig. 9 grid) ----
N = 512;  M = 128;                       % N subcarriers/delay, M symbols/Doppler
deltaf = 15e3;  Ts = 1/deltaf;           % subcarrier spacing / symbol duration
fc = 4e9;  c = 3e8;  v = 120/3.6;        % 4 GHz, 120 km/h
nu_max = v*fc/c;                         % max Doppler ~444 Hz
M_order = 4;                             % 4-QAM
constellation = qammod((0:M_order-1).', M_order, UnitAveragePower=true);

%% ---- EVA channel (integer delay taps).  Toggle reduced/full (cf. compute_mcrb_mc) --
use_full_eva = true;   % false = reduced EVA (3 strongest taps), true = full EVA (9 paths)
if use_full_eva
    % Full EVA, 9 paths (3GPP TS 36.104, Annex B)
    eva_delays_ns = [0 30 150 310 370 710 1090 1730 2510];
    eva_pow_db    = [0 -1.5 -1.4 -3.6 -0.6 -9.1 -7.0 -12.0 -16.9];
else
    % Reduced EVA, 3 strongest taps (0, -0.6, -1.4 dB)
    eva_delays_ns = [0 150 370];
    eva_pow_db    = [0 -1.4 -0.6];
end
evaname       = 'full EVA';  if ~use_full_eva, evaname = 'reduced EVA'; end
P             = numel(eva_delays_ns);
delay_taps    = round(eva_delays_ns*1e-9 * N * deltaf);   % full: [0 0 1 2 3 5 8 13 19]
tau           = delay_taps / (N*deltaf);                  % rounded delays [s]
amp           = 10.^(eva_pow_db/20);
amp           = amp / norm(amp);                          % sum|alpha|^2 = 1
l_tau         = ceil(max(eva_delays_ns*1e-9) * N * deltaf);

%% ---- Embedded pilot (OTFS): single impulse + full Doppler guard (eq. 8) -----
l_p = l_tau;               % pilot delay row
k_p = M/2;                 % pilot Doppler col
guard_rows = 1 : (2*l_tau + 1);            % delay rows 0..2*l_tau (MATLAB idx)
ce_rows    = (l_p+1) : (2*l_tau + 1);      % rows the pilot response occupies
data_mask_otfs = true(N, M);
data_mask_otfs(guard_rows, :) = false;     % guard band carries no data

%% ---- OFDM periodic pilots: OVERHEAD-matched to the OTFS guard+pilot ----------
% Same non-data fraction as OTFS ((2*l_tau+1)*M cells): the OFDM pilot sublattice
% (dN_p,dM_p) is the closest scattered pattern to that cell count under the CE
% oversampling limits (dN_p < N/(2 l_tau), dM_p < M/(2 k_nu)).  Pilot POWER is
% matched too -- the pilots share the OTFS single-pilot energy (see a_p_of below).
k_nu = ceil(nu_max * M / deltaf);
otfs_nondata = numel(guard_rows) * M;
[dN_p, dM_p, n_pilot_cells] = match_pilot_spacing(N, M, otfs_nondata, l_tau, k_nu);
n_pilots = 1:dN_p:N;   m_pilots = 1:dM_p:M;
pilot_mask_ofdm = false(N, M);
pilot_mask_ofdm(n_pilots, m_pilots) = true;
data_mask_ofdm  = ~pilot_mask_ofdm;
W_ici = 8;                                 % banded MMSE half-bandwidth (>= few*k_nu)

%% ---- Sweep / detector -------------------------------------------------------
SNR_dbs      = 10:2:18;
SNRp_dB      = [40 50];                    % OTFS CE pilot SNRs
SNRp_ofdm_dB = 40;                         % OFDM CE aggregate pilot SNR (matched)
nMonteCarlo  = 20;
n_iter_mp    = 30;
damping      = 0.3;                        % weight on the fresh MP message; 0.7 diverges on
                                           % some EVA realizations (-> BER floor), 0.3 is stable
E_dopp       = 10;                         % Doppler-spread half-window kept by MP detector

nSNR = numel(SNR_dbs);  nP = numel(SNRp_dB);
ber_ofdm_ideal = zeros(nSNR, nMonteCarlo);
ber_ofdm_ce    = zeros(nSNR, nMonteCarlo);
ber_otfs_ideal = zeros(nSNR, nMonteCarlo);
ber_otfs_ce    = zeros(nSNR, nMonteCarlo, nP);

fprintf('BER vs SNR_d  (%s ICI, N=%d M=%d, P=%d, %d trials)  full guard, frac. Doppler\n', ...
        evaname, N, M, P, nMonteCarlo);
fprintf('OFDM pilots dN_p=%d dM_p=%d -> %d cells (%.2f%%); OTFS non-data %.2f%% (power-matched)\n', ...
        dN_p, dM_p, n_pilot_cells, 100*n_pilot_cells/(N*M), 100*otfs_nondata/(N*M));

for si = 1:nSNR
    SNR_lin = db2pow(SNR_dbs(si));
    N0 = 1/SNR_lin;
    rng(2026);
    tic;
    parfor mc = 1:nMonteCarlo
        % ---- channel realization: fixed delays/powers, random Doppler+phase ----
        nu    = nu_max * cos(2*pi*rand(1,P) - pi);
        alpha = amp .* exp(1j*2*pi*rand(1,P));

        % ---- shared noise across OTFS variants (DD white) and OFDM (TF white) ----
        Z_DD = sqrt(N0/2) * (randn(N,M) + 1j*randn(N,M));
        Z_TF = sqrt(N0/2) * (randn(N,M) + 1j*randn(N,M));

        % ===================== OFDM (ideal, banded MMSE) =====================
        sym_of   = randi(M_order, N, M) - 1;
        X_of     = qammod(sym_of, M_order, UnitAveragePower=true);
        Y_of     = ofdm_mean(X_of, deltaf, Ts, alpha, tau, nu, true) + Z_TF;
        Xhat_ofi = ofdm_banded_mmse(Y_of, deltaf, Ts, alpha, tau, nu, N0, W_ici);
        txb = qamdemod(X_of(:),     M_order, OutputType="bit", UnitAveragePower=true);
        rxb = qamdemod(Xhat_ofi(:), M_order, OutputType="bit", UnitAveragePower=true);
        ber_ofdm_ideal(si,mc) = mean(txb ~= rxb);

        % ===================== OFDM (CE, periodic pilots, 1-tap) =============
        E_p_of   = db2pow(SNRp_ofdm_dB) * N0;              % OTFS-equivalent pilot energy
        a_p_of   = sqrt(E_p_of / n_pilot_cells);           % per-pilot amplitude (matched)
        X_ofc    = X_of;
        X_ofc(pilot_mask_ofdm) = a_p_of;                   % overwrite pilot cells
        Y_ofc    = ofdm_mean(X_ofc, deltaf, Ts, alpha, tau, nu, true) + Z_TF;
        H_ls     = Y_ofc(n_pilots, m_pilots) ./ X_ofc(n_pilots, m_pilots);
        F_H      = griddedInterpolant({n_pilots, m_pilots}, H_ls, 'linear', 'nearest');
        [Nq, Mq] = ndgrid(1:N, 1:M);
        H_hat    = F_H(Nq, Mq);
        Xhat_ofc = conj(H_hat) .* Y_ofc ./ (abs(H_hat).^2 + N0);   % 1-tap MMSE
        txb = qamdemod(X_ofc(data_mask_ofdm),    M_order, OutputType="bit", UnitAveragePower=true);
        rxb = qamdemod(Xhat_ofc(data_mask_ofdm), M_order, OutputType="bit", UnitAveragePower=true);
        ber_ofdm_ce(si,mc) = mean(txb ~= rxb);

        % ===================== Effective DD channel (exact OTFS I/O) ==========
        % ofdm_mean's TF-ICI model is NOT a pure DD circular convolution, so a
        % plain-circshift MP cannot match it (~8% irreducible residual -> the BER
        % floor/divergence).  OTFS therefore uses the exact delay-Doppler response
        % (Raviteja WCNC'18 eq. 13): Gtrue applies the channel (full Doppler
        % spread), Gdet is the E-truncated tap list the detector works with.
        Gtrue = build_effective_dd_eva(alpha, delay_taps, nu, N, M, deltaf, Inf);
        Gdet  = build_effective_dd_eva(alpha, delay_taps, nu, N, M, deltaf, E_dopp);

        % ===================== OTFS (ideal, known taps + MP) =================
        sym_ot = randi(M_order, N, M) - 1;
        X_ot   = qammod(sym_ot, M_order, UnitAveragePower=true);
        Y_ot   = apply_dd_channel_eva(X_ot, Gtrue) + Z_DD;

        [ri,ci,gi] = find(Gdet);
        dl_i = ri - 1;  dk_i = ci - 1;                 % absolute DD shifts (origin ref.)
        Xhat_oti = mp_detector(Y_ot, dl_i, dk_i, gi, N0, constellation, ...
                               n_iter_mp, true(N,M), zeros(N,M), damping);
        txb = qamdemod(X_ot(:),      M_order, OutputType="bit", UnitAveragePower=true);
        rxb = qamdemod(Xhat_oti(:),  M_order, OutputType="bit", UnitAveragePower=true);
        ber_otfs_ideal(si,mc) = mean(txb ~= rxb);

        % ===================== OTFS (CE, embedded pilot + MP) ================
        ber_ce_row = zeros(1,nP);
        for pj = 1:nP
            x_p = sqrt(db2pow(SNRp_dB(pj)) * N0);
            X_ce = zeros(N,M);
            X_ce(data_mask_otfs) = X_ot(data_mask_otfs);   % reuse data symbols
            X_ce(l_p+1, k_p+1)   = x_p;                    % embedded pilot
            Y_ce = apply_dd_channel_eva(X_ce, Gtrue) + Z_DD;

            % threshold CE on the guard-protected pilot rows (|y| >= 3 sigma)
            h_ce = zeros(N,M);
            h_ce(ce_rows,:) = Y_ce(ce_rows,:) / x_p;
            h_ce(abs(h_ce) < 3*sqrt(N0)/x_p) = 0;
            [rc,cc,gc] = find(h_ce);
            if isempty(gc), ber_ce_row(pj) = NaN; continue; end
            dl_c = mod(rc - (l_p+1), N);  dk_c = mod(cc - (k_p+1), M);

            known = zeros(N,M);  known(l_p+1,k_p+1) = x_p;
            Xhat_ce = mp_detector(Y_ce, dl_c, dk_c, gc, N0, constellation, ...
                                  n_iter_mp, data_mask_otfs, known, damping);
            txb = qamdemod(X_ce(data_mask_otfs),    M_order, OutputType="bit", UnitAveragePower=true);
            rxb = qamdemod(Xhat_ce(data_mask_otfs), M_order, OutputType="bit", UnitAveragePower=true);
            ber_ce_row(pj) = mean(txb ~= rxb);
        end
        ber_otfs_ce(si,mc,:) = ber_ce_row;
    end
    fprintf('  SNR_d=%2d dB | OFDM=%.3e OFDMce=%.3e OTFS=%.3e', SNR_dbs(si), ...
        mean(ber_ofdm_ideal(si,:)), mean(ber_ofdm_ce(si,:)), mean(ber_otfs_ideal(si,:)));
    for pj = 1:nP
        fprintf(' OTFSce%d=%.3e', SNRp_dB(pj), mean(ber_otfs_ce(si,:,pj),'omitnan'));
    end
    fprintf(' | %.1fs\n', toc);
end

%% ---- aggregate + plot -------------------------------------------------------
BER_ofdm_ideal = mean(ber_ofdm_ideal, 2);
BER_ofdm_ce    = mean(ber_ofdm_ce,    2);
BER_otfs_ideal = mean(ber_otfs_ideal, 2);
BER_otfs_ce    = squeeze(mean(ber_otfs_ce, 2, 'omitnan'));   % nSNR x nP

save(fullfile(thisdir,'results','ber_ce_ici_results.mat'), ...
     'SNR_dbs','SNRp_dB','SNRp_ofdm_dB','BER_ofdm_ideal','BER_ofdm_ce', ...
     'BER_otfs_ideal','BER_otfs_ce','N','M','nMonteCarlo','use_full_eva','evaname');

figure; hold on; grid on;
semilogy(SNR_dbs, BER_ofdm_ideal, '-d', 'LineWidth', 2);
semilogy(SNR_dbs, BER_ofdm_ce,    '-x', 'LineWidth', 2);
semilogy(SNR_dbs, BER_otfs_ideal, '-o', 'LineWidth', 2);
mk = {'--s','--^','--v'};
for pj = 1:nP
    semilogy(SNR_dbs, BER_otfs_ce(:,pj), mk{pj}, 'LineWidth', 2);
end
set(gca,'YScale','log');
xlabel('SNR_d [dB]');  ylabel('BER');
ce_lbl = arrayfun(@(x) sprintf('OTFS (CE, SNR_p=%d dB)', x), SNRp_dB, 'UniformOutput', false);
legend([{'OFDM (ideal, banded MMSE)', ...
         sprintf('OFDM (CE, SNR_p=%d dB)', SNRp_ofdm_dB), ...
         'OTFS (ideal, MP)'}, ce_lbl], 'Location', 'southwest');
title(sprintf('%s ICI channel — CE + demodulation BER', evaname));

% Detection helpers (match_pilot_spacing, build_effective_dd_eva,
% apply_dd_channel_eva, ofdm_banded_mmse) now live in eva/utils and are shared
% with mcrb_ce_ici.m.
