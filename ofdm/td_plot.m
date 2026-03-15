%% OFDM Baseband Signal - Time Domain Plot with AWGN Channel & Demodulation
% Parameters
N      = 64;          % Number of subcarriers
M      = 4;           % Number of OFDM symbols
deltaf = 15e3;        % Subcarrier spacing (Hz)
T      = 1/deltaf;    % OFDM symbol duration (s)
CP     = T/4;         % Cyclic prefix duration

% Sampling parameters
L    = 10;             % Oversampling factor (quasi-continuous)
Nfft = L * N;          % = 640 — FFT size must be >= N
Fs   = Nfft / T;       % Sampling frequency derived from Nfft
Ts   = 1/Fs;           % Sampling period
Ncp  = round(CP/Ts);  % Samples per cyclic prefix

% SNR
SNR_dB = -20;          % SNR per subcarrier in dB (try 5, 10, 20)
SNR    = 10^(SNR_dB/10);

fprintf('=== OFDM System Parameters ===\n');
fprintf('Subcarriers N         = %d\n',    N);
fprintf('OFDM Symbols M        = %d\n',    M);
fprintf('Subcarrier spacing    = %.1f kHz\n', deltaf/1e3);
fprintf('Symbol duration T     = %.4f ms\n',  T*1e3);
fprintf('Samples per symbol    = %d\n',    Nfft);
fprintf('CP samples            = %d\n',    Ncp);
fprintf('SNR                   = %d dB\n', SNR_dB);

%% QPSK Modulation
qpsk_map = (1/sqrt(2)) * [1+1j, -1+1j, -1-1j, 1-1j]; % Gray-coded

rng(42);
data_idx  = randi([0 3], N, M);        % Random symbol indices [N x M]
qpsk_syms = qpsk_map(data_idx + 1);    % QPSK symbols [N x M]

%% --- TRANSMITTER: OFDM Modulation (IFFT + CP) ---
s_total   = [];
sym_start = zeros(1, M);   % Track where each symbol starts (for plotting)

for m = 1:M
    freq_domain  = qpsk_syms(:, m);
    time_symbol  = sqrt(N) * ifft(freq_domain, Nfft);
    cp           = time_symbol(end - Ncp + 1 : end);
    ofdm_symbol  = [cp; time_symbol];

    sym_start(m) = length(s_total);
    s_total      = [s_total; ofdm_symbol];
end

%% --- AWGN CHANNEL ---
% Signal power (per sample)
sig_power  = mean(abs(s_total).^2);

% Noise variance: sigma^2 = sig_power / (SNR * N)
% Factor N because SNR is defined per subcarrier
noise_var  = sig_power / (SNR * N);
noise      = sqrt(noise_var/2) * (randn(size(s_total)) + 1j*randn(size(s_total)));

r_total    = s_total + noise;   % Received signal

fprintf('\nSignal power          = %.4f\n', sig_power);
fprintf('Noise variance        = %.6f\n',  noise_var);

%% --- RECEIVER: CP Removal + FFT Demodulation ---
rx_syms    = zeros(N, M);   % Recovered QPSK symbols
sym_len    = Nfft + Ncp;

for m = 1:M
    % Extract the m-th OFDM symbol from received signal
    seg = r_total(sym_start(m)+1 : sym_start(m)+sym_len);

    % Remove cyclic prefix
    seg_no_cp = seg(Ncp+1 : end);       % Keep only the useful part [Nfft x 1]

    % FFT: time domain → frequency domain
    freq_rx = (1/sqrt(N)) * fft(seg_no_cp, Nfft);
    rx_syms(:, m) = freq_rx(1:N);   % First N bins = your N subcarriers
end

%% --- QPSK DEMODULATION (Nearest Neighbour Decision) ---
% For each received symbol, find the closest QPSK constellation point
rx_flat  = rx_syms(:);                  % Flatten to vector
tx_flat  = qpsk_syms(:);

% Compute distances to all 4 constellation points
dist     = abs(rx_flat - qpsk_map).^2; % [N*M x 4] distance matrix
[~, dec_idx] = min(dist, [], 2);        % Index of minimum distance

% Detected symbols and original symbols
det_syms = qpsk_map(dec_idx).';        % Detected QPSK symbols
tx_syms  = tx_flat;

% Symbol Error Rate
SER = mean(det_syms ~= tx_syms);
fprintf('\n=== Demodulation Results ===\n');
fprintf('Symbol Error Rate     = %.4f (%d errors out of %d)\n', ...
        SER, sum(det_syms ~= tx_syms), N*M);

%% ========== PLOTS ==========
t = (0 : length(s_total)-1) * Ts * 1e3;  % Time in ms

figure('Name','OFDM Full Chain','Position',[50 30 1200 900]);

%-- 1: Transmitted Signal (Real) --%
subplot(4,2,1);
plot(t, real(s_total), 'b', 'LineWidth', 0.7);
hold on;
for m = 0:M-1
    t_s = sym_start(m+1)*Ts*1e3;
    t_c = t_s + Ncp*Ts*1e3;
    patch([t_s t_c t_c t_s], [-3 -3 3 3], 'y', ...
          'FaceAlpha',0.2,'EdgeColor','none');
end
xlabel('Time (ms)'); ylabel('Amplitude');
title('Tx Signal — Real Part'); grid on; xlim([t(1) t(end)]);
legend('Re\{s(t)\}','CP','Location','northeast'); hold off;

%-- 2: Received Signal (Real) --%
subplot(4,2,2);
plot(t, real(r_total), 'r', 'LineWidth', 0.7);
xlabel('Time (ms)'); ylabel('Amplitude');
title(sprintf('Rx Signal — Real Part (SNR = %d dB)', SNR_dB));
grid on; xlim([t(1) t(end)]);

%-- 3: Transmitted Signal (Imag) --%
subplot(4,2,3);
plot(t, imag(s_total), 'b', 'LineWidth', 0.7);
xlabel('Time (ms)'); ylabel('Amplitude');
title('Tx Signal — Imaginary Part'); grid on; xlim([t(1) t(end)]);

%-- 4: Received Signal (Imag) --%
subplot(4,2,4);
plot(t, imag(r_total), 'r', 'LineWidth', 0.7);
xlabel('Time (ms)'); ylabel('Amplitude');
title('Rx Signal — Imaginary Part'); grid on; xlim([t(1) t(end)]);

%-- 5: Instantaneous Power Comparison --%
subplot(4,2,[5 6]);
plot(t, abs(s_total).^2, 'b', 'LineWidth', 0.7); hold on;
plot(t, abs(r_total).^2, 'Color', [1 0 0 0.6], 'LineWidth', 0.7);
yline(mean(abs(s_total).^2),'b--','Tx Mean','LineWidth',1.2);
yline(mean(abs(r_total).^2),'r--','Rx Mean','LineWidth',1.2);
xlabel('Time (ms)'); ylabel('|s(t)|^2');
title('Instantaneous Power — Tx vs Rx');
legend('Tx','Rx','Location','northeast'); grid on; xlim([t(1) t(end)]); hold off;

%-- 6: Tx Constellation --%
subplot(4,2,7);
scatter(real(tx_syms), imag(tx_syms), 40, 'b', 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
% scatter(real(qpsk_map), imag(qpsk_map), 120, 'k', 'p', 'LineWidth', 1.5);
xlabel('In-Phase'); ylabel('Quadrature');
title('Tx Constellation (QPSK)');
grid on; axis equal; xlim([-1.5 1.5]); ylim([-1.5 1.5]);
xline(0,'k--','LineWidth',0.5); yline(0,'k--','LineWidth',0.5);
legend('Tx Symbols','Ideal Points'); hold off;

%-- 7: Rx Constellation (after FFT demod) --%
subplot(4,2,8);
% Color by error: green = correct, red = error
correct_mask = (det_syms == tx_syms);
scatter(real(rx_flat(correct_mask)),  imag(rx_flat(correct_mask)),  ...
        40, [0 0.7 0], 'filled', 'MarkerFaceAlpha', 0.5); hold on;
scatter(real(rx_flat(~correct_mask)), imag(rx_flat(~correct_mask)), ...
        60, 'r', 'filled', 'MarkerFaceAlpha', 0.8);
% scatter(real(qpsk_map), imag(qpsk_map), 120, 'k', 'p', 'LineWidth', 1.5);
xlabel('In-Phase'); ylabel('Quadrature');
title(sprintf('Rx Constellation — SER = %.4f', SER));
grid on; axis equal; xlim([-1.5 1.5]); ylim([-1.5 1.5]);
xline(0,'k--','LineWidth',0.5); yline(0,'k--','LineWidth',0.5);
legend('Correct','Error','Ideal','Location','northeast'); hold off;

sgtitle(sprintf('OFDM Full Chain: N=%d, M=%d, \\Deltaf=%.0f kHz, SNR=%d dB, QPSK', ...
        N, M, deltaf/1e3, SNR_dB), 'FontSize',13,'FontWeight','bold');