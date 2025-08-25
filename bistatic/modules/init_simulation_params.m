function params = init_simulation_params()
    c = physconst('LightSpeed');
    
    N = 400;                   % Number of OFDM subcarriers
    M = 60;                    %#ok<*NASGU> % Number of OFDM symbols
    N_fft = 4096;
    M_fft = 256;

    fc = 28e9;                 % Carrier frequency
    delta_f = 120e3;           % Subcarrier spacing
    BW = N * delta_f;
    T = 1 / delta_f;            % Symbol duration
    Tcp = T * 0.07;
    Tsym = Tcp + T;

    PT_db = 20;                 % Transmit power (20 dBm)
    PT = db2pow(PT_db) * 1e-3;  % Transmit power (0.1 W)
    % NT = 8;                   % Number of TX antennas
    noise_fig = 8;              % Noise figure (dB)
    N0 = -174;                  % Noise PSD (dbm/Hz)
    sigma2 = db2pow(N0 + noise_fig) * 1e-3 * N * delta_f;
    
    K = 2;                     % Number of targets
    known_data = false;       % Genie-aided or not
    pilot_ratio = 0.02;       % Fraction of pilots
    modulation = 'QPSK';      % Data modulation
    
    % Locations (x,y): TX, RX, and 2 targets
    pT = [0, 0];
    pR = [50, 0];
    targets = [56.9, 10; 79.4, 7]; % targets(:, k) k = 1,2,...,K
    velocities = [1.4, -2.2; 2.2, -13.7];
    rcs_dB = [4.9, 15];       % Target RCS in dBsm
    rcs = db2pow(rcs_dB) * 1e-3;
    
    % Derived parameters
    lambda = c / fc;

    taus = [0.13281, 0.39953] * 1e-6;
    nus = [375.8984, 377.3253];
    alphas = [5.3963e-7, 5.4666e-8];

    tau_res = 1/BW;
    nu_res = 1/(M*Tsym);

    tau_idx = tau2idx(taus, Tsym, N_fft);
    nu_idx = nu2idx(nus, delta_f, M_fft);

    % % Allocate arrays
    % delays = zeros(K, 1);
    % dopplers = zeros(K, 1);
    % delay_idxs = zeros(K, 1);
    % doppler_idxs = zeros(K, 1);
    % 
    % for k = 1:K
    %     pk = targets(k, :);
    %     vk = velocities(k, :);
    % 
    %     d1 = norm(pT - pk);
    %     d2 = norm(pk - pR);
    %     tau = (d1 + d2) / c;
    %     delays(k) = tau;
    % 
    %     los_vector = (pR - pk) / norm(pR - pk);
    %     v_rel = dot(vk, los_vector);
    %     nu = 2 * v_rel / lambda;
    %     dopplers(k) = nu;
    % 
    %     % Find closest bin indices
    %     delay_bins = (0:N-1) * delay_res;
    %     doppler_bins = (-M/2:M/2-1) * doppler_res;
    % 
    %     [~, delay_idx] = min(abs(delay_bins - tau));
    %     [~, doppler_idx] = min(abs(doppler_bins - nu));
    % 
    %     delay_idxs(k) = delay_idx;
    %     doppler_idxs(k) = doppler_idx + M/2; % to match fftshifted layout
    % end
    params = v2struct();
end
