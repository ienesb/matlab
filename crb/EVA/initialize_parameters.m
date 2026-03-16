% EVA Channel model    1, 5, 3, 2, 4, 7, 6, 8, 9

nMonteCarlo = 1;

N = 512; % number of subcarriers
M = 128; % number of symbols

M_order = 4;

c = 3e8;
fc = 4*1e9; % 4 GHz
lambda = c / fc;

deltaf = 15*1e3; % 15 kHz
Tcp = 0;
T = 1/deltaf;
Ts = T + Tcp;

sigma2 = 1;

SNR_dbs = 10:2:18;
SNR_lins = db2pow(SNR_dbs);

time_delays = [0, 30 ,150, 310, 370, 710, 1090, 1730, 2510] .* 1e-9; % ns
relative_powers = [0, -1.5, -1.4, -3.6, -0.6, -9.1, -7, -12, -16.9]; % dB

P = length(time_delays);

velocity = 120 * 1e3 / 3600; % kmph
nu_max = velocity * fc / c;
% doppler_shifts = nu_max * cos(2*pi*rand(1, 9) - pi);

% Convert relative dB powers to linear power, then normalize
relative_powers_lin = 10.^(relative_powers/10);
normalized_power_profile = relative_powers_lin / sum(relative_powers_lin);

% Convert to voltage amplitudes for the channel taps
path_gains_amp = sqrt(normalized_power_profile);

% SNR_idx = 0;
% SNR_db = SNR_dbs(SNR_idx);
% path_gains_db = SNR_db + relative_powers;
% path_gains = db2pow(path_gains_db); % add random phase shift !!!!!!

% Delay-Doppler tap indices (integer, worst-case)
l_tau = ceil(max(time_delays) * N * deltaf);  % max delay tap  = 20
k_nu  = ceil(nu_max * M / deltaf);            % max Doppler tap = 4

% Per-path integer delay taps: l_i = round(tau_i * N * deltaf)
% Bin width in delay = 1/(N*deltaf).  EVA: [0 0 1 2 3 5 8 13 19]
delay_taps = round(time_delays * N * deltaf);

% Per-path integer Doppler taps: k_i = round(nu_i * M / deltaf)
% Bin width in Doppler = deltaf/M.  nu_i is random per MC iteration,
% drawn as nu_i = nu_max * cos(theta_i).  Worst-case magnitude = k_nu.
% doppler_taps = round(current_doppler * M / deltaf);  % computed inside MC loop
