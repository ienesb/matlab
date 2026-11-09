function p = params_fig9()
% PARAMS_FIG9  System / channel / simulation parameters for reproducing
% Fig. 9 of Raviteja et al., "Embedded Pilot-Aided Channel Estimation for
% OTFS in Delay-Doppler Channels," IEEE TVT 2019 (the "Embedded" paper).
%
% Fig. 9: BER vs SNR_d, fractional Doppler, FULL guard symbols.
%   N = 128, M = 512, l_tau = 20, 4-QAM, threshold T = 3*sigma_p,
%   UE speed 120 km/h (k_nu = 4), EVA channel.
%   Curves: OTFS-CE @ SNR_p = {40,45,50} dB, OTFS-Ideal, OFDM-Ideal.
%
% All equation/section references below are to the Embedded paper unless
% prefixed [MP], which refers to Raviteja et al., "Low-complexity iterative
% detection for OTFS," IEEE WCNC 2018 (the MP detector reference).

% ---- OTFS grid (Sec. IV) ----------------------------------------------
p.N      = 128;          % Doppler bins (time slots), index k = 0..N-1
p.M      = 512;          % delay bins (subcarriers), index l = 0..M-1
p.deltaf = 15e3;         % subcarrier spacing [Hz]
p.T      = 1/p.deltaf;   % symbol duration [s]
p.fc     = 4e9;          % carrier frequency [Hz]
p.c      = 299792458;    % speed of light [m/s]

% ---- Modulation -------------------------------------------------------
p.Qbits  = 2;            % 4-QAM -> 2 bits/symbol
p.alphabet = ([-1-1i, -1+1i, 1-1i, 1+1i].')/sqrt(2);  % unit-power 4-QAM
% Gray map for the four constellation points (bit pairs), used for BER.
p.bitmap = [0 0; 0 1; 1 0; 1 1];

% ---- Mobility / Doppler (Sec. IV) -------------------------------------
p.v_kmh  = 120;                       % UE speed [km/h]
v        = p.v_kmh/3.6;               % [m/s]
p.nu_max = v*p.fc/p.c;                % max Doppler shift [Hz] (~444 Hz)
p.k_nu   = round(p.nu_max*p.N*p.T);   % max Doppler tap (= 4 at 120 km/h)

% ---- EVA channel power-delay profile (3GPP TS 36.104, Annex B) --------
% 9 taps. Delays rounded to nearest delay tap l = round(tau*M*deltaf).
p.eva_delays_ns = [0 30 150 310 370 710 1090 1730 2510];
p.eva_pow_dB    = [0 -1.5 -1.4 -3.6 -0.6 -9.1 -7.0 -12.0 -16.9];
p.l_tau         = 20;                 % max delay tap used for guard sizing

% ---- Embedded pilot arrangement: fractional, FULL guards (eq. 8) ------
%   x[k,l] = x_p              at (k_p,l_p)
%          = 0                for 0<=k<=N-1, l_p-l_tau <= l <= l_p+l_tau
%          = x_d              otherwise
p.l_p = p.l_tau;                      % pilot delay index (>= l_tau, <= M-1-l_tau)
p.k_p = p.N/2;                        % pilot Doppler index (Doppler is circular)

% ---- Detector settings ------------------------------------------------
p.E         = 10;     % Doppler-spread half-window kept by detector ([MP] Fig.3
                      % shows E=10 saturates performance for k_nu=4).
p.tap_tol   = 0.02;   % prune ideal-detector taps below this fraction of the
                      % peak |G| (drops negligible Dirichlet IDI tails; the
                      % estimated channel is pruned by the 3*sigma threshold).
p.mp_iter   = 20;     % MP iterations ([MP] Sec. IV: converges within ~20)
p.mp_damp   = 0.7;    % damping factor Delta ([MP] Fig. 4 optimum)

% ---- Channel estimation threshold (Sec. IV-B, Fig. 9 uses T = 3*sigma_p)
% sigma^2 = 1 in all sims; |x_p|^2 = SNR_p. The estimate is hhat = y/x_p,
% thresholded at 3*sigma_p with sigma_p = sigma/|x_p|, which is identical to
% thresholding |y| at 3*sigma = 3 (since sigma = 1). We threshold |y|.
p.thr_sigma = 3;      % T = 3*sigma, sigma = 1

% ---- Sweep ------------------------------------------------------------
p.SNRd_dB   = 10:2:18;          % data SNR sweep (x-axis)
p.SNRp_dB   = [40 45 50];       % pilot SNR curves
p.n_trials  = 50;               % Monte-Carlo frames per point
p.seed      = 2019;             % RNG seed for reproducibility

% ---- OFDM reference ---------------------------------------------------
p.cp_len = round(2.6e-6*p.M*p.deltaf);  % CP length (~20 samples, paper: 2.6 us)
end
