function run_sim(varargin)
% RUN_SIM  Reproduce Fig. 9 of Raviteja et al., "Embedded Pilot-Aided
% Channel Estimation for OTFS in Delay-Doppler Channels," IEEE TVT 2019.
%
%   Fig. 9: BER vs SNR_d, fractional Doppler, FULL guard symbols.
%   Curves: OTFS with proposed channel estimation @ SNR_p = {40,45,50} dB,
%           OTFS-Ideal (known channel), OFDM-Ideal (known channel).
%
% Usage:
%   run_sim                       % full run (params_fig9 defaults)
%   run_sim('n_trials', 8)        % override any parameter field
%   run_sim('quick', true)        % fast smoke run (few trials, small E)
%
% Trials run in parallel (parfor) with one reproducible RNG substream each.
% Outputs: results.mat (BER tables + params) and fig_9.png in this folder.

p = params_fig9();

% ---- optional overrides ----------------------------------------------
quick = false; nworkers = 8;
for a = 1:2:numel(varargin)
    switch varargin{a}
        case 'quick',    quick    = varargin{a+1};
        case 'nworkers', nworkers = varargin{a+1};
        otherwise,       p.(varargin{a}) = varargin{a+1};
    end
end
if quick
    p.n_trials = 3; p.E = 6; p.mp_iter = 15;
end

N = p.N; M = p.M;
nSNRd = numel(p.SNRd_dB);
nSNRp = numel(p.SNRp_dB);

% ---- embedded-pilot frame masks (eq. 8, full guards) -----------------
guardMask = false(N, M);
guardMask(:, (p.l_p-p.l_tau:p.l_p+p.l_tau)+1) = true;  % full Doppler guard band
pilotMask = false(N, M); pilotMask(p.k_p+1, p.l_p+1) = true;
guardMask(pilotMask) = false;
knownMask = guardMask | pilotMask;                     % non-data positions
dataMask  = ~knownMask;
nDataEmb  = nnz(dataMask);
nDataFull = N*M;

% ---- parallel pool ----------------------------------------------------
if nworkers > 1 && isempty(gcp('nocreate'))
    parpool(min(nworkers, feature('numcores')));
end

% ---- per-trial error accumulators (sliced for parfor) ----------------
ei = zeros(p.n_trials, nSNRd);          bi = zeros(p.n_trials, nSNRd);          % OTFS ideal
ece = zeros(p.n_trials, nSNRp, nSNRd);  bce = zeros(p.n_trials, nSNRp, nSNRd);  % OTFS CE
eof = zeros(p.n_trials, nSNRd);         bof = zeros(p.n_trials, nSNRd);         % OFDM ideal

fprintf('Fig.9 sim: N=%d M=%d, %d trials, E=%d, iters=%d\n', ...
        N, M, p.n_trials, p.E, p.mp_iter);
t0 = tic;

parfor tr = 1:p.n_trials
    % reproducible independent RNG per trial
    s = RandStream('Philox','Seed',p.seed); s.Substream = tr;
    RandStream.setGlobalStream(s);

    ch    = gen_eva_channel(p);
    Gtrue = build_effective_dd(ch, p, Inf);     % exact channel to apply
    Gdet  = build_effective_dd(ch, p, p.E);     % detector's truncated view
    taps_ideal = dd_to_taps(Gdet, p.tap_tol);

    le = zeros(1,nSNRd);  lb = zeros(1,nSNRd);              % local OTFS ideal
    lce = zeros(nSNRp,nSNRd); lbce = zeros(nSNRp,nSNRd);    % local OTFS CE

    for si = 1:nSNRd
        SNRd = 10^(p.SNRd_dB(si)/10);
        dataSyms = p.alphabet*sqrt(SNRd);
        sigma2 = 1;

        % ---------- OTFS Ideal (full data frame, known channel) -------
        bits = randi([0 1], nDataFull, p.Qbits);
        X = reshape(qam4_mod(bits, dataSyms, p.bitmap), N, M);
        Y = apply_dd_channel(X, Gtrue, sigma2);
        Xhat = mp_detect_otfs(Y, taps_ideal, dataSyms, sigma2, ...
                  false(N,M), zeros(N,M), p.mp_iter, p.mp_damp);
        bhat = qam4_demod(Xhat(:), dataSyms, p.bitmap);
        le(si) = le(si) + sum(bits(:) ~= bhat(:));
        lb(si) = lb(si) + numel(bits);

        % ---------- OTFS with channel estimation (per SNR_p) ----------
        for pj = 1:nSNRp
            xp = sqrt(10^(p.SNRp_dB(pj)/10));
            bitsD = randi([0 1], nDataEmb, p.Qbits);
            X = zeros(N, M);
            X(dataMask)  = qam4_mod(bitsD, dataSyms, p.bitmap);
            X(pilotMask) = xp;
            Y = apply_dd_channel(X, Gtrue, sigma2);

            Ghat = estimate_channel(Y, p, xp, p.thr_sigma);  % T = 3*sigma
            taps_est = dd_to_taps(Ghat);
            if isempty(taps_est.gain), continue; end
            knownVal = zeros(N,M); knownVal(pilotMask) = xp;
            Xhat = mp_detect_otfs(Y, taps_est, dataSyms, sigma2, ...
                      knownMask, knownVal, p.mp_iter, p.mp_damp);
            bhat = qam4_demod(Xhat(dataMask), dataSyms, p.bitmap);
            lce(pj,si)  = lce(pj,si)  + sum(bitsD(:) ~= bhat(:));
            lbce(pj,si) = lbce(pj,si) + numel(bitsD);
        end
    end

    % ---------- OFDM Ideal (all SNR_d at once; H_n reused) ------------
    [lbf, lof] = ofdm_ideal(ch, p, 10.^(p.SNRd_dB/10));

    ei(tr,:) = le;   bi(tr,:) = lb;
    ece(tr,:,:) = lce; bce(tr,:,:) = lbce;
    eof(tr,:) = lof; bof(tr,:) = lbf;
    fprintf('  trial %d/%d done (%.1f s elapsed)\n', tr, p.n_trials, toc(t0));
end

% ---- aggregate BER ----------------------------------------------------
BER_ideal = sum(ei,1)  ./ max(sum(bi,1),1);
BER_ofdm  = sum(eof,1) ./ max(sum(bof,1),1);
err_ce = squeeze(sum(ece,1)); bit_ce = squeeze(sum(bce,1));
BER_ce = err_ce ./ max(bit_ce,1);

save(fullfile(fileparts(mfilename('fullpath')),'results.mat'), ...
     'p','BER_ideal','BER_ce','BER_ofdm','err_ce','bit_ce');

plot_fig9(p, BER_ideal, BER_ce, BER_ofdm);
fprintf('Done in %.1f s. Saved results.mat and fig_9.png\n', toc(t0));
end
