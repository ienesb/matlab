%% OFDM with Multipath Channel, Pilot Estimation, Equalization & Y Comparison
clear; clc;

%% ── System Parameters ─────────────────────────────────────────────────────
N      = 64;          % Number of subcarriers
M      = 8;           % Total OFDM symbols (data + pilot)
deltaf = 15e3;        % Subcarrier spacing (Hz)
T      = 1/deltaf;    % Symbol duration (s)
CP     = T/4;         % Cyclic prefix duration

L      = 10;          % Oversampling factor (quasi-continuous)
Nfft   = L * N;       % FFT size = 640
Fs     = Nfft / T;    % Sampling frequency
Ts     = 1/Fs;        % Sampling period
Ncp    = round(CP/Ts);% CP samples

SNR_dB = 20;          % SNR per subcarrier (dB)
SNR    = 10^(SNR_dB/10);

%% ── Pilot Configuration ───────────────────────────────────────────────────
pilot_sym_idx = [1, 4, 8];                        % Pilot OFDM symbol positions
data_sym_idx  = setdiff(1:M, pilot_sym_idx);      % Data symbol positions
n_data        = length(data_sym_idx);

fprintf('=== OFDM System Parameters ===\n');
fprintf('Subcarriers N         = %d\n',    N);
fprintf('OFDM Symbols M        = %d\n',    M);
fprintf('Oversampling L        = %d\n',    L);
fprintf('FFT size              = %d\n',    Nfft);
fprintf('Subcarrier spacing    = %.1f kHz\n', deltaf/1e3);
fprintf('Symbol duration T     = %.4f ms\n',  T*1e3);
fprintf('Samples per symbol    = %d\n',    Nfft);
fprintf('CP samples            = %d\n',    Ncp);
fprintf('SNR                   = %d dB\n', SNR_dB);
fprintf('Pilot symbols         : %s\n',    num2str(pilot_sym_idx));
fprintf('Data  symbols         : %s\n\n',  num2str(data_sym_idx));

%% ── Multipath Channel Definition ─────────────────────────────────────────
% Each path: [complex gain,  Doppler shift (Hz),  delay (s)]
% Constraint: max delay < CP duration to avoid ISI
paths = [ ...
    1.0 + 0.0j,   10,   0;           % Path 1: LOS, 10 Hz Doppler, no delay
    0.6 - 0.3j,   25,   T/20;        % Path 2: 25 Hz Doppler, T/20 delay
    0.3 + 0.4j,  -15,   T/10;        % Path 3: -15 Hz Doppler, T/10 delay
    0.2 - 0.1j,   40,   T/8;         % Path 4: 40 Hz Doppler, T/8 delay
];
n_paths   = size(paths, 1);
gains     = paths(:,1);
dopplers  = paths(:,2);
delays    = paths(:,3);
delay_smp = round(delays / Ts);

fprintf('─── Channel Paths ───\n');
fprintf('Path | Gain (mag) | Doppler (Hz) | Delay (samples)\n');
for p = 1:n_paths
    fprintf('  %d  |   %.3f     |     %+.0f      |      %d\n', ...
        p, abs(gains(p)), dopplers(p), delay_smp(p));
end
fprintf('CP length = %d samples  |  Max delay = %d samples\n\n', ...
    Ncp, max(delay_smp));

%% ── QPSK Constellation ────────────────────────────────────────────────────
qpsk_map = (1/sqrt(2)) * [1+1j, -1+1j, -1-1j, 1-1j];  % Gray-coded

% Fixed pilot pattern (same across all pilot symbols)
rng(0);
pilot_idx_pattern = randi([0 3], N, 1);
pilot_freq        = qpsk_map(pilot_idx_pattern + 1).';  % [N x 1]

% Random data symbols
rng(42);
data_idx  = randi([0 3], N, n_data);
data_freq = qpsk_map(data_idx + 1);                     % [N x n_data]

% Assemble full frequency-domain frame [N x M]
X = zeros(N, M);
X(:, pilot_sym_idx) = repmat(pilot_freq, 1, length(pilot_sym_idx));
X(:, data_sym_idx)  = data_freq;

%% ── TRANSMITTER: Zero-Padded IFFT + Cyclic Prefix ────────────────────────
s_total   = [];
sym_start = zeros(1, M);

for m = 1:M
    % Zero-pad N subcarriers into Nfft-point IFFT for oversampling
    freq_zp      = [X(:,m); zeros(Nfft-N, 1)];
    time_sym     = sqrt(Nfft) * ifft(freq_zp, Nfft);
    cp           = time_sym(end-Ncp+1 : end);
    sym_start(m) = length(s_total);
    s_total      = [s_total; cp; time_sym];
end

total_smp = length(s_total);
t_ax      = (0:total_smp-1) * Ts * 1e3;   % Time axis in ms

%% ── MULTIPATH CHANNEL (Time-Varying) ─────────────────────────────────────
% r[n] = sum_p gain_p * exp(j*2pi*f_dp*n*Ts) * s[n - delay_smp_p]
n_vec  = (0:total_smp-1).';
r_chan = zeros(total_smp, 1);

for p = 1:n_paths
    d             = delay_smp(p);
    doppler_phase = exp(1j * 2*pi * dopplers(p) * n_vec * Ts);
    s_delayed     = [zeros(d,1); s_total(1:end-d)];
    r_chan        = r_chan + gains(p) * doppler_phase .* s_delayed;
end

%% ── AWGN NOISE ────────────────────────────────────────────────────────────
sig_power = mean(abs(r_chan).^2);
noise_var = sig_power / (SNR * N);
noise     = sqrt(noise_var/2) * (randn(total_smp,1) + 1j*randn(total_smp,1));
r_total   = r_chan + noise;

fprintf('=== Power ===\n');
fprintf('Signal power  = %.4f\n',  sig_power);
fprintf('Noise variance= %.6f\n\n',noise_var);

%% ── RECEIVER: CP Removal + FFT ────────────────────────────────────────────
sym_len = Nfft + Ncp;
Y       = zeros(N, M);   % Received freq-domain frame [N x M]

for m = 1:M
    seg      = r_total(sym_start(m)+1 : sym_start(m)+sym_len);
    seg_nocp = seg(Ncp+1 : end);                    % Remove CP
    freq_rx  = (1/sqrt(Nfft)) * fft(seg_nocp, Nfft);
    Y(:, m)  = freq_rx(1:N);                        % Keep N subcarrier bins
end

%% ── FREQUENCY-DOMAIN CHANNEL MATRIX H [N x M] ────────────────────────────
% H[k,m] = sum_p gain_p * exp(j*2pi*f_dp*t_m) * exp(-j*2pi*k*delay_p/Nfft)
% where t_m = centre time of useful part (after CP) of symbol m
H_matrix = zeros(N, M);

for m = 1:M
    t_m = (sym_start(m) + Ncp + Nfft/2) * Ts;   % Centre time of symbol m
    for k = 1:N
        H_km = 0;
        for p = 1:n_paths
            doppler_phase = exp(1j  * 2*pi * dopplers(p) * t_m);
            delay_phase   = exp(-1j * 2*pi * (k-1) * delay_smp(p) / Nfft);
            H_km          = H_km + gains(p) * doppler_phase * delay_phase;
        end
        H_matrix(k, m) = H_km;
    end
end

%% ── GENERATE Y_fd = X .* H ────────────────────────────────────────────────
Y_fd = X .* H_matrix;    % Ideal frequency-domain channel model [N x M]

%% ── COMPARISON: Y vs Y_fd ─────────────────────────────────────────────────
diff_matrix = abs(Y - Y_fd);

fprintf('=== Y Comparison: Time-Domain Pipeline vs Freq-Domain Model ===\n');
fprintf('Mean absolute difference : %.6f\n',   mean(diff_matrix(:)));
fprintf('Max  absolute difference : %.6f\n',   max(diff_matrix(:)));
fprintf('Mean |Y|                 : %.6f\n',   mean(abs(Y(:))));
fprintf('Relative error           : %.4f%%\n\n', ...
        100 * mean(diff_matrix(:)) / mean(abs(Y(:))));

%% ── CHANNEL ESTIMATION: LS at Pilot Symbols ───────────────────────────────
H_pilot = zeros(N, length(pilot_sym_idx));

for pi = 1:length(pilot_sym_idx)
    m             = pilot_sym_idx(pi);
    H_pilot(:,pi) = Y(:,m) ./ X(:,m);    % Per-subcarrier LS estimate
end

%% ── CHANNEL INTERPOLATION: Linear in Time ────────────────────────────────
H_est = zeros(N, M);

for k = 1:N
    H_est(k,:) = interp1(pilot_sym_idx, H_pilot(k,:), ...
                         1:M, 'linear', 'extrap');
end

%% ── EQUALIZATION: Zero-Forcing ────────────────────────────────────────────
X_eq = Y ./ H_est;       % ZF: divide out estimated channel [N x M]

%% ── QPSK DECISION ─────────────────────────────────────────────────────────
X_eq_data = X_eq(:, data_sym_idx);   % [N x n_data]
X_tx_data = X(:,  data_sym_idx);     % [N x n_data]

rx_flat  = X_eq_data(:);
tx_flat  = X_tx_data(:);

dist     = abs(rx_flat - qpsk_map).^2;
[~, dec_idx] = min(dist, [], 2);
det_syms = qpsk_map(dec_idx).';

SER      = mean(det_syms ~= tx_flat);
n_errors = sum(det_syms ~= tx_flat);

fprintf('=== Demodulation Results ===\n');
fprintf('Data symbols    : %d\n',     N * n_data);
fprintf('Symbol errors   : %d\n',     n_errors);
fprintf('SER             : %.5f\n\n', SER);

%% ══════════════════ FIGURE 1: Time-Domain Signals ═════════════════════════
figure('Name','OFDM: Time-Domain Signals','Position',[30 30 1300 500]);

subplot(1,3,1);
plot(t_ax, real(s_total), 'b', 'LineWidth', 0.6); hold on;
sym_len_ms = sym_len * Ts * 1e3;
for m = 0:M-1
    ts = sym_start(m+1)*Ts*1e3;
    tc = ts + Ncp*Ts*1e3;
    te = ts + sym_len_ms;
    patch([ts tc tc ts],[-3 -3 3 3],'y','FaceAlpha',0.2,'EdgeColor','none');
    if ismember(m+1, pilot_sym_idx)
        patch([ts te te ts],[-3 -3 3 3],'g','FaceAlpha',0.12,'EdgeColor','none');
    end
end
xlabel('Time (ms)'); ylabel('Amplitude');
title('Tx Signal — Real Part');
legend('s(t)','CP','Pilot Sym','Location','ne');
grid on; xlim([t_ax(1) t_ax(end)]); hold off;

subplot(1,3,2);
plot(t_ax, real(r_chan), 'Color',[0.8 0.4 0], 'LineWidth', 0.6);
xlabel('Time (ms)'); ylabel('Amplitude');
title('After Multipath Channel — Real Part');
grid on; xlim([t_ax(1) t_ax(end)]);

subplot(1,3,3);
plot(t_ax, real(r_total), 'r', 'LineWidth', 0.6);
xlabel('Time (ms)'); ylabel('Amplitude');
title(sprintf('Rx Signal + AWGN (SNR=%d dB) — Real Part', SNR_dB));
grid on; xlim([t_ax(1) t_ax(end)]);

sgtitle('OFDM: Tx → Multipath → AWGN', 'FontSize', 12, 'FontWeight', 'bold');

%% ══════════════════ FIGURE 2: Y vs Y_fd Comparison ══════════════════════
figure('Name','Y Comparison','Position',[30 30 1300 750]);

subplot(2,3,1);
imagesc(1:M, 1:N, abs(Y));
colorbar; colormap(gca,'jet');
xlabel('OFDM Symbol'); ylabel('Subcarrier');
title('|Y| — Time-Domain Pipeline');

subplot(2,3,2);
imagesc(1:M, 1:N, abs(Y_fd));
colorbar; colormap(gca,'jet');
xlabel('OFDM Symbol'); ylabel('Subcarrier');
title('|Y_{fd}| — Freq-Domain Model (X .* H)');

subplot(2,3,3);
imagesc(1:M, 1:N, diff_matrix);
colorbar; colormap(gca,'hot');
xlabel('OFDM Symbol'); ylabel('Subcarrier');
title('|Y - Y_{fd}| — Absolute Difference');

subplot(2,3,4);
scatter(real(Y(:)), real(Y_fd(:)), 20, 'b', 'filled', 'MarkerFaceAlpha', 0.4);
hold on;
ref = linspace(min(real(Y(:))), max(real(Y(:))), 100);
plot(ref, ref, 'r--', 'LineWidth', 1.5);
xlabel('Re\{Y\} — Time-Domain'); ylabel('Re\{Y_{fd}\} — Freq-Domain');
title('Scatter: Real Parts'); legend('Subcarriers','Ideal y=x'); grid on; hold off;

subplot(2,3,5);
scatter(imag(Y(:)), imag(Y_fd(:)), 20, [0.1 0.6 0.2], 'filled', 'MarkerFaceAlpha', 0.4);
hold on;
ref = linspace(min(imag(Y(:))), max(imag(Y(:))), 100);
plot(ref, ref, 'r--', 'LineWidth', 1.5);
xlabel('Im\{Y\} — Time-Domain'); ylabel('Im\{Y_{fd}\} — Freq-Domain');
title('Scatter: Imaginary Parts'); legend('Subcarriers','Ideal y=x'); grid on; hold off;

subplot(2,3,6);
bar(1:M, mean(diff_matrix,1), 'FaceColor',[0.3 0.5 0.8]); hold on;
for pi = pilot_sym_idx
    bar(pi, mean(diff_matrix(:,pi)), 'FaceColor',[0.9 0.6 0.1]);
end
xlabel('OFDM Symbol Index'); ylabel('Mean |Y - Y_{fd}|');
title('Mean Difference per Symbol');
legend('Data Symbol','Pilot Symbol'); grid on; hold off;

sgtitle(sprintf('Y Comparison | Rel. Error = %.4f%%', ...
    100*mean(diff_matrix(:))/mean(abs(Y(:)))), ...
    'FontSize', 12, 'FontWeight', 'bold');

%% ══════════════════ FIGURE 3: Channel Estimation ════════════════════════
% True channel at each symbol centre
H_true = zeros(N, M);
for m = 1:M
    t_m = (sym_start(m) + Ncp + Nfft/2) * Ts;
    for k = 1:N
        H_km = 0;
        for p = 1:n_paths
            H_km = H_km + gains(p) ...
                * exp(1j  * 2*pi * dopplers(p) * t_m) ...
                * exp(-1j * 2*pi * (k-1) * delay_smp(p) / Nfft);
        end
        H_true(k,m) = H_km;
    end
end

figure('Name','Channel Estimation','Position',[30 30 1300 750]);

subplot(2,3,1);
imagesc(1:M, 1:N, abs(H_est));
colorbar; colormap(gca,'jet');
hold on;
for pi = pilot_sym_idx, xline(pi,'w--','LineWidth',1.5); end
hold off;
xlabel('OFDM Symbol'); ylabel('Subcarrier');
title('|H_{est}| — Estimated Magnitude');

subplot(2,3,2);
imagesc(1:M, 1:N, angle(H_est)*180/pi);
colorbar; colormap(gca,'hsv');
hold on;
for pi = pilot_sym_idx, xline(pi,'w--','LineWidth',1.5); end
hold off;
xlabel('OFDM Symbol'); ylabel('Subcarrier');
title('\angle H_{est} (deg) — Estimated Phase');

subplot(2,3,3);
plot(1:N, abs(H_true(:,1)), 'b-',  'LineWidth', 1.5); hold on;
plot(1:N, abs(H_est(:,1)),  'r--', 'LineWidth', 1.5);
xlabel('Subcarrier Index'); ylabel('|H[k]|');
title('True vs Estimated Channel (Symbol 1)');
legend('H_{true}','H_{est}','Location','best');
grid on; hold off;

est_err = abs(H_true - H_est);

subplot(2,3,4);
imagesc(1:M, 1:N, est_err);
colorbar; colormap(gca,'hot');
xlabel('OFDM Symbol'); ylabel('Subcarrier');
title('|H_{true} - H_{est}| — Estimation Error');

subplot(2,3,5);
bar(1:N, mean(est_err,2), 'FaceColor',[0.2 0.5 0.8]);
xlabel('Subcarrier Index'); ylabel('Mean |Error|');
title('Mean Ch. Est. Error per Subcarrier'); grid on;

subplot(2,3,6);
bar(1:M, mean(est_err,1), 'FaceColor',[0.8 0.3 0.2]); hold on;
for pi = pilot_sym_idx
    bar(pi, mean(est_err(:,pi)), 'FaceColor',[0.1 0.7 0.3]);
end
xlabel('OFDM Symbol Index'); ylabel('Mean |Error|');
title('Mean Ch. Est. Error per Symbol');
legend('Data Symbol','Pilot Symbol'); grid on; hold off;

sgtitle('Channel Estimation: LS + Linear Interpolation', ...
    'FontSize', 12, 'FontWeight', 'bold');

%% ══════════════════ FIGURE 4: Constellation & Equalization ══════════════
figure('Name','Constellation & Equalization','Position',[30 30 1100 400]);

m_show = data_sym_idx(1);   % First data symbol for per-symbol plots

subplot(1,3,1);
scatter(real(Y(:,m_show)), imag(Y(:,m_show)), 30, 'r', 'filled', ...
        'MarkerFaceAlpha', 0.5); hold on;
scatter(real(qpsk_map), imag(qpsk_map), 100,'k','p','LineWidth',1.5);
xline(0,'k--','LineWidth',0.5); yline(0,'k--','LineWidth',0.5);
xlabel('In-Phase'); ylabel('Quadrature');
title(sprintf('Before Equalization (Sym %d)', m_show));
grid on; axis equal; xlim([-2 2]); ylim([-2 2]); hold off;

subplot(1,3,2);
scatter(real(X_eq(:,m_show)), imag(X_eq(:,m_show)), 30, [0 0.6 0.9], ...
        'filled','MarkerFaceAlpha', 0.5); hold on;
scatter(real(qpsk_map), imag(qpsk_map), 100,'k','p','LineWidth',1.5);
xline(0,'k--','LineWidth',0.5); yline(0,'k--','LineWidth',0.5);
xlabel('In-Phase'); ylabel('Quadrature');
title(sprintf('After ZF Equalization (Sym %d)', m_show));
grid on; axis equal; xlim([-2 2]); ylim([-2 2]); hold off;

subplot(1,3,3);
correct_mask = (det_syms == tx_flat);
scatter(real(rx_flat(correct_mask)),  imag(rx_flat(correct_mask)),  ...
        25, [0 0.75 0], 'filled','MarkerFaceAlpha',0.4); hold on;
scatter(real(rx_flat(~correct_mask)), imag(rx_flat(~correct_mask)), ...
        50, 'r','filled','MarkerFaceAlpha',0.9);
scatter(real(qpsk_map), imag(qpsk_map), 120,'k','p','LineWidth',1.5);
xline(0,'k--','LineWidth',0.5); yline(0,'k--','LineWidth',0.5);
xlabel('In-Phase'); ylabel('Quadrature');
title(sprintf('Full Rx Constellation — SER=%.4f (%d err)', SER, n_errors));
legend('Correct','Error','Ideal','Location','ne');
grid on; axis equal; xlim([-2 2]); ylim([-2 2]); hold off;

sgtitle(sprintf('OFDM Full Chain | N=%d, M=%d, %d Paths | SNR=%d dB | SER=%.4f', ...
    N, M, n_paths, SNR_dB, SER), 'FontSize', 12, 'FontWeight', 'bold');