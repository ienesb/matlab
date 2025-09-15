function params = init_simulation_params(rcs2_dB, rho, nIter, monteCarlo, is_genie, is_data_only)
    % c = physconst('LightSpeed');
    c = 3e8;
    
    N = 400;                   % Number of OFDM subcarriers
    M = 60;                    % Number of OFDM symbols (words)
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

    NT = 8;                     % Number of TX antennas
    PT_db = 20;                 % Transmit power (20 dBm)
    PT = db2pow(PT_db) * 1e-3;  % Transmit power (0.1 W)
    noise_fig = 8;              % Noise figure (dB)
    N0 = -174;                  % Noise PSD (dbm/Hz)
    sigma2 = db2pow(N0 + noise_fig) * 1e-3 * N * delta_f;
    
    if is_genie
        pilot_ratio = 1;       % Fraction of pilots
    else
        pilot_ratio = rho;
    end

    % modulation = 'QPSK';      % Data modulation
    M_order = 4;
    
    % Locations (x,y): TX, RX, and 2 targets
    pos_Tx = [0, 0];
    pos_Rx = [50, 0];
    positions = [56.9, 10; 79.4, 7]; % targets(k, :) k = 1,2,...,K
    velocities = [1.4, -2.2; 2.2, -13.7];
    rcs_dB = [4.9, rcs2_dB];       % Target RCS in dBsm
    rcs = db2pow(rcs_dB);

    K = length(positions);

    angles = [10, 5];

    ref_target_idx = 2;
    
    lambda = c / fc;
    ft = conj(geta(10, lambda, NT));
    ft = ft ./ norm(ft);

    taus = zeros(K, 1);
    nus = zeros(K, 1);
    alphas = zeros(K, 1);

    for k = 1:K
        pk = positions(k, :);
        vk = velocities(k, :);

        tau = pos2tau(pk, pos_Tx, pos_Rx);
        taus(k) = tau;

        % los_vector = (pR - pk) / norm(pR - pk);
        % v_rel = dot(vk, los_vector);
        % nu = 2 * v_rel / lambda;
        nu = velocity2nu(vk.', pk.', pos_Tx.', pos_Tx.', lambda);
        nus(k) = nu;

        alpha = lambda * sqrt(rcs(k)) / ((4*pi)^1.5 * norm(pk - pos_Tx) * norm(pk - pos_Rx));
        alpha = alpha * sqrt(PT);
        alpha = alpha * abs(geta(angles(k), lambda, NT).' * ft);
        alphas(k) = alpha;
    end


    tau_idx = tau2idx(taus.', Tsym, N_fft);
    nu_idx = nu2idx(nus.', delta_f, M_fft);

    tau_res = 1/BW;
    nu_res = 1/(M*Tsym);

    delay_array = getDelayArray(Tsym, N_fft);
    doppler_array = getDopplerArray(delta_f, M_fft, is_fftshifted);

    params = v2struct(c, N, M, N_fft, M_fft, fc, delta_f, BW, T, Tcp, ...
        Tsym, is_fftshifted, NT, PT_db, PT, noise_fig, N0, sigma2, nIter, ...
        pilot_ratio, M_order, pos_Tx, pos_Rx, positions, velocities, ...
        rcs_dB, rcs, K, angles, lambda, taus, nus, alphas, tau_idx, nu_idx, tau_res, nu_res, ...
        delay_array, doppler_array, ref_target_idx, monteCarlo, is_data_only, is_genie);
end
