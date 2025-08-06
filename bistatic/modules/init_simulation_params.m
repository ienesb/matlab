
function params = init_simulation_params()
% Simulation parameters for bistatic ISAC

params.N = 400;                   % Number of OFDM subcarriers
params.M = 60;                    % Number of OFDM symbols
params.fc = 28e9;                 % Carrier frequency
params.delta_f = 120e3;           % Subcarrier spacing
params.BW = params.N * params.delta_f;
params.Tsym = 1 / params.delta_f; % Symbol duration
params.Tcp = 1.16e-6;             % Cyclic prefix duration
params.PT = 10^(20/10) / 1000;    % Transmit power (20 dBm)
params.NT = 8;                    % Number of TX antennas
params.noise_fig = 8;             % Noise figure (dB)
params.N0 = 10^(-174/10) / 1000;  % Noise PSD (W/Hz)
params.snr = 10;                  % SNR in dB
params.sigma2 = 10^(-params.snr/10);

params.K = 2;                     % Number of targets
params.known_data = false;       % Genie-aided or not
params.pilot_ratio = 0.05;       % Fraction of pilots
params.modulation = 'QPSK';      % Data modulation

% Locations (x,y): TX, RX, and 2 targets
params.pT = [0, 0];
params.pR = [50, 0];
params.targets = [56.9, 10; 79.4, 7];
params.velocities = [1.4, -2.2; 2.2, -13.7];
params.rcs_dB = [4.9, 1.5];       % Target RCS in dBsm

% Derived parameters
c = 3e8;                          % Speed of light
lambda = c / params.fc;
params.lambda = lambda;

% Delay/Doppler resolution
params.delay_res = 1 / params.BW;
params.doppler_res = 1 / (params.M * params.Tsym);

% Allocate arrays
params.delays = zeros(params.K, 1);
params.dopplers = zeros(params.K, 1);
params.delay_idxs = zeros(params.K, 1);
params.doppler_idxs = zeros(params.K, 1);

for k = 1:params.K
    pk = params.targets(k, :);
    vk = params.velocities(k, :);

    d1 = norm(params.pT - pk);
    d2 = norm(pk - params.pR);
    tau = (d1 + d2) / c;
    params.delays(k) = tau;

    los_vector = (params.pR - pk) / norm(params.pR - pk);
    v_rel = dot(vk, los_vector);
    nu = 2 * v_rel / lambda;
    params.dopplers(k) = nu;

    % Find closest bin indices
    delay_bins = (0:params.N-1) * params.delay_res;
    doppler_bins = (-params.M/2:params.M/2-1) * params.doppler_res;

    [~, delay_idx] = min(abs(delay_bins - tau));
    [~, doppler_idx] = min(abs(doppler_bins - nu));

    params.delay_idxs(k) = delay_idx;
    params.doppler_idxs(k) = doppler_idx + params.M/2; % to match fftshifted layout
end

end
