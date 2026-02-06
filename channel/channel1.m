clear;
close all;
clc;

%% parameters
nMonteCarlo = 1;

N = 70; % number of subcarriers
M = 50; % number of OFDM symbols

M_order = 4; % QPSK

c = 3e8;
fc = 30*1e9; % 30 GHz
lambda = c / fc;

deltaf = 200*1e3; % 200 kHz
% Tcp = 1*1e-6; % 1 us
Tcp = 0; % 1 us
T = 1/deltaf;
Tsym = T + Tcp; % Symbol duration
Ts = Tsym / 5000; % sampling rate for illustration

p_tx = [-40, 0];
p_rx = [0, 40];

p = [rand*20+80, rand*20-100]; % m

v = rand*60-30; % m/s
delta = rand*10-5; % degree

d_tx = norm(p - p_tx);
d_rx = norm(p - p_rx);
D = norm(p_tx - p_rx);

d_bis = (d_tx + d_rx);
v_bis = v * cosd(delta);

beta = acos((d_tx^2 + d_rx^2 - D^2)/(2*d_tx*d_rx));

tau_gt = d_bis / c;
tau_gt = round(tau_gt/Ts)*Ts;
nu_gt = 2*v_bis*cos(beta/2)/lambda;

% SNR_dbs = -30:30;
SNR_dbs = 0;
SNR_lins = db2pow(SNR_dbs);


%% Time domain signal generation

symbols = randi(M_order, N, M)-1;
X = qammod(symbols, M_order, UnitAveragePower=true);

ns = (0:(N-1)).';

t = [];
s = [];

for m = 0:(M-1)
    tm = (0:Ts:(Tsym-Ts)) + m*Tsym;
    sm = 1/sqrt(N) * sum(X(:, m+1) .* exp(1j*2*pi*ns*deltaf*tm), 1);

    t = [t, tm];
    s = [s, sm];
end

% figure;
% plot(t, real(s)); hold on;
% plot(t, imag(s));
% xlabel("Time (s)");
% ylabel("s(t)");
% grid on;

%% Apply doppler shift and time delay

y = s;

doppler_shift = exp(1j*2*pi*nu_gt*t);

y = doppler_shift .* y;

tau_gt = 0;
time_shift = round(tau_gt/Ts);

t_a = (Ts:Ts:((time_shift)*Ts)) + (t(end));
t = [t, t_a];
s = [s, zeros(1, length(t_a))];

y_new = zeros(1, length(s));
y_new((time_shift+1):length(y)+time_shift) = y;

y = y_new;

% figure;
% plot(t, real(s)); hold on;
% plot(t, real(y));
% xlabel("Time (s)");
% ylabel("s(t)");
% grid on;


Y = zeros(N, M);
for m = 0:(M-1)

    tm_idx = (1:(Tsym/Ts)) + m*(Tsym/Ts);
    tm = t(tm_idx);
    ym = y(tm_idx);

    % Correlate against subcarrier basis
    E = exp(-1j*2*pi*(ns*deltaf)*tm);  % N x Nuse

    % Y(:, m+1) = sqrt(1/N) * (Ts/Tsym) * E * (ym.');
    Y(:, m+1) = E * (ym.') / 600;

end

ns = (0:(N-1)).';
ms = 0:(M-1);

H = exp(-1j*2*pi*deltaf*ns*0) * exp(1j*2*pi*T*ms*nu_gt);

Y2 = H .* X;

diff = mean(abs(Y2(:) - Y(:)));


%% sadasdasdad

b = exp(-1j*2*pi*ns*deltaf*tau_gt);
c = exp(-1j*2*pi*ms*Tsym*nu_gt);

D = exp(1j*2*pi*fc*T*ns/N*nu_gt);
D = diag(D);
% D = eye(N);

F = dftmtx(N) / sqrt(N);
% F = eye(N);

Y3 = D * F' * (X .* (b * conj(c)));
% Y4 = F * Y3;
Y3 = Y3(:);
Y3 = Y3.';

t3 = linspace(0, M*Tsym, N*M+1);
t3 = t3(1:end-1);

figure;
plot(t, imag(y)); hold on;
plot(t3, imag(Y3), "*");
grid on;
