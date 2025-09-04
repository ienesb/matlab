function params = init_simulation_params(rho)
    c = physconst('LightSpeed');
    
    N = 400;                   % Number of OFDM subcarriers
    M = 60;                    %#ok<*NASGU> % Number of OFDM symbols
    N_fft = 4096;
    M_fft = 256;

    fc = 28e9;                 % Carrier frequency
    delta_f = 120e3;           % Subcarrier spacing
    BW = N * delta_f;
    T = 1 / delta_f;            % Symbol duration
    % Tcp = T * 0.07; 
    Tcp = 0;
    Tsym = Tcp + T;

    is_fftshifted = 1;

    NT = 8;                   % Number of TX antennas
    PT_db = 20;                 % Transmit power (20 dBm)
    PT = db2pow(PT_db) * 1e-3 * NT;  % Transmit power (0.1 W)
    noise_fig = 8;              % Noise figure (dB)
    N0 = -174;                  % Noise PSD (dbm/Hz)
    sigma2 = db2pow(N0 + noise_fig) * 1e-3 * N * delta_f;
    
    K = 2;                     % Number of targets
    known_data = false;       % Genie-aided or not
    pilot_ratio = 0.05;       % Fraction of pilots
    modulation = 'QPSK';      % Data modulation
    
    % Locations (x,y): TX, RX, and 2 targets
    pT = [0, 0];
    pR = [50, 0];
    positions = [56.9, 10; 79.4, 7]; % targets(k, :) k = 1,2,...,K
    velocities = [1.4, -2.2; 2.2, -13.7];
    rcs_dB = [4.9, 2];       % Target RCS in dBsm
    rcs = db2pow(rcs_dB);

    rcs(2) = rho;

    ref_target_idx = 2;

    % lambda * sqrt(rcs) / ((4*pi)^1.5 * norm(targets(1, :) - pT) * norm(targets(1, :) - pR))
    
    % Derived parameters
    lambda = c / fc;

    % taus = [0.13281, 0.39953] * 1e-6;
    % nus = [375.8984, 377.3253];
    % alphas = [5.3963e-7, 5.4666e-8];

    taus = zeros(K, 1);
    nus = zeros(K, 1);
    alphas = zeros(K, 1);

    for k = 1:K
        pk = positions(k, :);
        vk = velocities(k, :);

        tau = pos2tau(pk, pT, pR);
        taus(k) = tau;

        % los_vector = (pR - pk) / norm(pR - pk);
        % v_rel = dot(vk, los_vector);
        % nu = 2 * v_rel / lambda;
        nu = velocity2nu(vk.', pk.', pT.', pR.', lambda);
        nus(k) = nu;

        alpha = lambda * sqrt(rcs(k)) / ((4*pi)^1.5 * norm(pk - pT) * norm(pk - pR));
        alpha = alpha * sqrt(PT);
        alphas(k) = alpha;
    end

    % alphas = alphas .* PT;

    tau_idx = tau2idx(taus.', Tsym, N_fft);
    nu_idx = nu2idx(nus.', delta_f, M_fft);

    tau_res = 1/BW;
    nu_res = 1/(M*Tsym);

    delay_array = getDelayArray(Tsym, N_fft);
    doppler_array = getDopplerArray(delta_f, M_fft, is_fftshifted);

    monteCarlo = 500;

    params = v2struct(c, N, M, N_fft, M_fft, fc, delta_f, BW, T, Tcp, Tsym, is_fftshifted, PT_db, PT, noise_fig, N0, sigma2, K, known_data, pilot_ratio, modulation, pT, pR, positions, velocities, rcs_dB, rcs, lambda, taus, nus, alphas, tau_idx, nu_idx, tau_res, nu_res, delay_array, doppler_array, ref_target_idx, monteCarlo);
end
