function params = init_params()
    c = physconst('LightSpeed');
    
    N = 400;                   % Number of OFDM subcarriers
    M = 60;                    % Number of OFDM symbols
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
    PT = db2pow(PT_db) * 1e-3;  % Transmit power (0.1 W)
    noise_fig = 8;              % Noise figure (dB)
    N0 = -174;                  % Noise PSD (dbm/Hz)
    sigma2 = db2pow(N0 + noise_fig) * 1e-3 * N * delta_f;
    
    K = 1;                   
    taus = 2.3323e-07;
    nus = 0;
    alphas = 1;
    modulation = 'QPSK';      % Data modulation
    params = 0;
    % params = v2struct(c, N, M, N_fft, M_fft, fc, delta_f, BW, T, Tcp, Tsym, is_fftshifted, NT, PT_db, PT, noise_fig, N0, sigma2, K, modulation, , rcs, lambda, taus, nus, alphas, tau_idx, nu_idx, tau_res, nu_res, delay_array, doppler_array, ref_target_idx);
end
