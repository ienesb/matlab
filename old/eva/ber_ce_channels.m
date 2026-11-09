% ber_ce_channels.m
% --------------------------------------------------------------------------
% BER vs SNR for OFDM-ICI and CP-OTFS-ICI with the realistic CE + demodulation
% chain of mcrb_ce_ici.m -- NO MCRB, BER only.
%
% FAIRNESS (asserted at run time, printed per channel):
%   * identical NON-DATA overhead: OTFS guard+pilot cells == OFDM pilot cells
%     (match_pilot_spacing), so the data:overhead ratio is the same.
%   * identical TOTAL pilot energy: E_pil for both.  OTFS concentrates it in the
%     single embedded pilot cell (the DD channel is sparse -> one pilot suffices);
%     OFDM spreads it over its pilot lattice (it must sample the TF grid).  That
%     split is the waveforms' own design, not a budget asymmetry.
%   * pilot power is DATA-referenced (E_pil = n_pilot_cells * boost, data cells
%     are unit power), NOT noise-referenced.  The old SNRp_dB=50 convention made
%     the pilots shrink as 1/SNR while data stayed at unit power, so data leakage
%     progressively swamped the OFDM pilots -> the BER curve turned UP past 32 dB.
%     Data-referencing removes that artifact and keeps both waveforms matched.
%
% OFDM channel estimation (ce_mode column of `configs`):
%   'interp' : LS at pilots -> LINEAR interpolation.  Original mcrb_ce_ici chain.
%              Floors hard: the channel rotates ~78 deg between frequency pilots
%              and linear interp chords that arc (see diag_ofdm_floor.m).
%   'dft'    : LS at pilots -> LS fit of the (l_tau+1) DELAY TAPS across the comb
%              -> exact reconstruction at every subcarrier -> linear interp in
%              time.  Valid: the comb meets the delay-Nyquist limit l_tau < N/dN_p.
%   'dftdd'  : 'dft' followed by DECISION-DIRECTED refinement -- re-fit the taps
%              using ALL N subcarriers (pilots known + data hard-decisions), per
%              OFDM symbol.  Fitting l_tau+1 taps from N equations instead of
%              from n_pilots gives a large processing gain, and each symbol gets
%              its own fit so no time interpolation is needed.
% OTFS: embedded-pilot threshold CE -> message-passing (MP) detector.
%
% NOTE the diversity wall: uncoded OFDM puts one symbol on one subcarrier = one
% flat fade -> diversity order 1 -> BER ~ 1/(4*SNR) (~2.5e-4 at 30 dB).  No CE or
% equalizer fix beats that; only coding/spreading would.  ber_vs_doppler.m shows
% the OFDM-vs-OTFS gap survives even at ZERO Doppler, i.e. it is diversity, not ICI.
% --------------------------------------------------------------------------
clear; clc; close all;
thisdir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisdir, 'utils'));
resdir = fullfile(thisdir, 'results');
warning('off','MATLAB:nearlySingularMatrix');

%% ---- toggles ---------------------------------------------------------------
channel_set    = {'eva','paperv2comms'};
SNR_dbs        = -10:1:35;
% Scale applied to each channel's nu_max:  0 = NO Doppler, 0.01 = very small,
% 1 = nominal (vehicular EVA / paperv2).  At nu=0 the ICI kernel collapses to a
% delta, so no-Doppler is automatically ICI-free -> OFDM's 1-tap MMSE is exact
% and any residual gap vs OTFS is pure diversity, not ICI (cf. ber_vs_doppler.m).
doppler_scale  = 0;
pilot_boost_dB = 0;      % total pilot energy vs the same #cells of unit-power data
nTrials        = 20;
seed_base      = 9100;
n_iter_mp      = 30;   damping = 0.3;
n_dd           = 2;      % decision-directed CE refinement passes ('dftdd')
M_order        = 4;
% Auto-name so a low-Doppler run does not overwrite the nominal one.
if     doppler_scale == 1, dtag = '';
elseif doppler_scale == 0, dtag = '_nodopp';
else,  dtag = sprintf('_dopp%gx', doppler_scale);  end
results_name   = ['ber_ce_channels' dtag '_results.mat'];

% {display name, waveform, OFDM ce_mode}
configs = { 'OFDM (LS+interp)', 'ofdm', 'interp'
            'OFDM (DFT CE)',    'ofdm', 'dft'
            'OFDM (DFT+DD CE)', 'ofdm', 'dftdd'
            'CP-OTFS (MP)',     'otfs', ''      };

nbits         = log2(M_order);
constellation = qammod((0:M_order-1).', M_order, UnitAveragePower=true);
nC = size(configs,1);  nSNR = numel(SNR_dbs);  nCh = numel(channel_set);
ber = nan(nSNR, nC, nCh);
n_data_bits = zeros(nC, nCh);
tAll = tic;

for chi = 1:nCh
    chname = channel_set{chi};
    switch lower(chname)
        case 'eva'
            N = 512;  M = 128;  deltaf = 15e3;  Ts = 1/deltaf;
            fc = 4e9;  c = 3e8;  v = 120/3.6;  nu_max = v*fc/c;
            delay_taps = round([0 150 370]*1e-9 * N * deltaf);
            amp = (10.^([0 -1.4 -0.6]/20)).';  amp = amp/norm(amp);
            dopp_mode = 'jakes';
        case 'paperv2comms'
            B = 160e6;  N = 128;  M = 16;  deltaf = B/N;  Ts = 1/deltaf;
            delay_taps = [0 2 5 8];
            amp = (10.^([0 -1.4 -3.6 -6.0]/20)).';  amp = amp/norm(amp);
            nu_max = 0.5/(M*Ts);
            dopp_mode = 'uniform';
        otherwise
            error('Unknown channel "%s"', chname);
    end
    nu_max = nu_max * doppler_scale;      % 0 = no Doppler (and hence no ICI)
    tau   = (delay_taps/(N*deltaf)).';
    P     = numel(delay_taps);
    l_tau = max(delay_taps);
    k_nu  = max(1, ceil(nu_max * M / deltaf));

    l_p = l_tau;  k_p = M/2;
    guard_rows = 1:(2*l_tau+1);
    ce_rows    = (l_p+1):(2*l_tau+1);
    assert(max(guard_rows) <= N, ...
        'OTFS guard (2*l_tau+1=%d) exceeds N=%d for "%s".', max(guard_rows), N, chname);
    data_mask_otfs = true(N,M);  data_mask_otfs(guard_rows,:) = false;

    otfs_nondata = numel(guard_rows) * M;
    [dN_p, dM_p, n_pilot_cells] = match_pilot_spacing(N, M, otfs_nondata, l_tau, k_nu);
    n_pilots = 1:dN_p:N;  m_pilots = 1:dM_p:M;
    pilot_mask_ofdm = false(N,M);  pilot_mask_ofdm(n_pilots, m_pilots) = true;
    data_mask_ofdm  = ~pilot_mask_ofdm;

    % ---- matched pilot budget (data-referenced, independent of N0) ----------
    E_pil = n_pilot_cells * db2pow(pilot_boost_dB);   % total pilot energy, both
    a_p   = sqrt(E_pil / n_pilot_cells);              % OFDM: spread over lattice
    x_p   = sqrt(E_pil);                              % OTFS: concentrated, 1 cell
    assert(abs(n_pilot_cells*a_p^2 - x_p^2) < 1e-9*E_pil, 'pilot energy mismatch');

    n_data_bits(strcmp(configs(:,2),'ofdm'), chi) = nnz(data_mask_ofdm) * nbits;
    n_data_bits(strcmp(configs(:,2),'otfs'), chi) = nnz(data_mask_otfs) * nbits;

    fprintf('\n[%s] N=%d M=%d df=%.3g Hz P=%d | l_tau=%d k_nu=%d nu_max=%.3g Hz\n', ...
            chname, N, M, deltaf, P, l_tau, k_nu, nu_max);
    fprintf('   FAIRNESS: non-data cells  OTFS=%d (%.2f%%)  OFDM=%d (%.2f%%)  [match=%d]\n', ...
            otfs_nondata, 100*otfs_nondata/(N*M), n_pilot_cells, 100*n_pilot_cells/(N*M), ...
            otfs_nondata == n_pilot_cells);
    fprintf('             total pilot E   OTFS=%.4e   OFDM=%.4e (%d x %.4e)\n', ...
            x_p^2, n_pilot_cells*a_p^2, n_pilot_cells, a_p^2);
    fprintf('             data cells      OTFS=%d   OFDM=%d\n', ...
            nnz(data_mask_otfs), nnz(data_mask_ofdm));
    fprintf('   pilot lattice (%d,%d); delay-Nyquist l_tau=%d < N/dN_p=%.1f -> %s\n', ...
            dN_p, dM_p, l_tau, N/dN_p, string(l_tau < N/dN_p));

    for ci = 1:nC
        modn = configs{ci,2};  ce_mode = configs{ci,3};
        for si = 1:nSNR
            N0 = 1/db2pow(SNR_dbs(si));
            be_t = zeros(1, nTrials);
            parfor t = 1:nTrials
                warning('off','MATLAB:nearlySingularMatrix');
                rs = RandStream('Threefry','Seed',seed_base);
                rs.Substream = 1000*si + t;      % paired across configs & channels
                RandStream.setGlobalStream(rs);

                switch dopp_mode
                    case 'jakes',   nu = (nu_max * cos(2*pi*rand(1,P) - pi)).';
                    case 'uniform', nu = (nu_max * (2*rand(1,P) - 1)).';
                end
                alpha = amp .* exp(1j*2*pi*rand(P,1));
                sym = randi(M_order, N, M) - 1;
                Xd  = qammod(sym, M_order, UnitAveragePower=true);
                Z   = sqrt(N0/2) * (randn(N,M) + 1j*randn(N,M));

                if strcmp(modn,'ofdm')
                    d = Xd;  d(pilot_mask_ofdm) = a_p;
                    Y = ofdm_mean(d, deltaf, Ts, alpha, tau, nu, true) + Z;
                    H_hat = ce_ofdm(Y, d, pilot_mask_ofdm, n_pilots, m_pilots, ...
                                    N, M, l_tau, N0, M_order, ce_mode, n_dd);
                    Xsoft = conj(H_hat) .* Y ./ (abs(H_hat).^2 + N0);
                    sh = qamdemod(Xsoft(data_mask_ofdm), M_order, UnitAveragePower=true);
                    be_t(t) = mean(int2bit(sym(data_mask_ofdm),nbits) ~= int2bit(sh,nbits),'all');
                else
                    d = zeros(N,M);  d(data_mask_otfs) = Xd(data_mask_otfs);
                    d(l_p+1, k_p+1) = x_p;
                    Gtrue = build_effective_dd_eva(alpha, delay_taps, nu, N, M, deltaf, Inf);
                    Y = apply_dd_channel_eva(d, Gtrue) + Z;
                    h_ce = zeros(N,M);
                    h_ce(ce_rows,:) = Y(ce_rows,:) / x_p;
                    h_ce(abs(h_ce) < 3*sqrt(N0)/x_p) = 0;
                    [rc,cc,gc] = find(h_ce);
                    if isempty(gc)
                        be_t(t) = 0.5;
                    else
                        dl_c = mod(rc-(l_p+1),N);  dk_c = mod(cc-(k_p+1),M);
                        known = zeros(N,M);  known(l_p+1,k_p+1) = x_p;
                        xhat = mp_detector(Y, dl_c, dk_c, gc, N0, constellation, ...
                                           n_iter_mp, data_mask_otfs, known, damping);
                        sh = qamdemod(xhat(data_mask_otfs), M_order, UnitAveragePower=true);
                        be_t(t) = mean(int2bit(sym(data_mask_otfs),nbits) ~= int2bit(sh,nbits),'all');
                    end
                end
            end
            ber(si,ci,chi) = mean(be_t);
            if mod(si,10)==1 || si==nSNR
                fprintf('   %-18s SNR=%3d dB | BER=%.3e   (%.0fs)\n', ...
                        configs{ci,1}, SNR_dbs(si), ber(si,ci,chi), toc(tAll));
            end
        end
    end
end

outpath = fullfile(resdir, results_name);
SNRp_dB = pilot_boost_dB;   % kept for plot-script compatibility
save(outpath,'ber','SNR_dbs','configs','channel_set','SNRp_dB','pilot_boost_dB', ...
     'nTrials','M_order','n_data_bits');
plot_ber_ce_channels(outpath);
fprintf('\nDone in %.0fs. Saved %s\n', toc(tAll), results_name);
warning('on','MATLAB:nearlySingularMatrix');

%% ================= local functions =================
function H_hat = ce_ofdm(Y, d, pilot_mask, n_pilots, m_pilots, N, M, l_tau, ...
                         N0, M_order, ce_mode, n_dd)
% OFDM channel estimation from the periodic pilot lattice.
    H_ls = Y(n_pilots, m_pilots) ./ d(n_pilots, m_pilots);      % Np x Mp
    [Nq, Mq] = ndgrid(1:N, 1:M);
    Afull = exp(-1j*2*pi * (0:N-1).' * (0:l_tau) / N);          % N x (l_tau+1)

    switch ce_mode
        case 'interp'                       % LS + linear interp (original)
            F_H   = griddedInterpolant({n_pilots, m_pilots}, H_ls, 'linear','nearest');
            H_hat = F_H(Nq, Mq);
            return
        case {'dft','dftdd'}                % delay-tap LS fit across the comb
            Abas  = exp(-1j*2*pi * (n_pilots(:)-1) * (0:l_tau) / N);   % Np x (l_tau+1)
            Hp    = Afull * (Abas \ H_ls);                             % N x Mp
            F_H   = griddedInterpolant({(1:N).', m_pilots}, Hp, 'linear','nearest');
            H_hat = F_H(Nq, Mq);
        otherwise
            error('Unknown ce_mode "%s"', ce_mode);
    end

    if strcmp(ce_mode,'dftdd')
        % Decision-directed refinement: re-fit the (l_tau+1) taps from ALL N
        % subcarriers per symbol (pilots known + data hard-decisions).  4-QAM and
        % the pilots are unit-modulus, so Y./dhat is well conditioned.  Fitting
        % l_tau+1 taps from N equations (vs numel(n_pilots)) is a large
        % processing gain, and every symbol gets its own fit -> no time interp.
        for it = 1:n_dd
            Xs   = conj(H_hat).*Y ./ (abs(H_hat).^2 + N0);
            dhat = qammod(qamdemod(Xs, M_order, UnitAveragePower=true), ...
                          M_order, UnitAveragePower=true);
            dhat(pilot_mask) = d(pilot_mask);            % pilots are known
            H_hat = Afull * (Afull \ (Y ./ dhat));       % LS tap fit, per symbol
        end
    end
end
