%% mcrb.m -> MCRB and MLE for pilot-only, data aided, and genie (MLE parts not complete!!!)
clear;
% close all;
clc;

options = optimoptions('fmincon', ...
              'Algorithm','interior-point', ...
              'Display','off', ...
              'MaxFunctionEvaluations',1e5, ...
              'FiniteDifferenceType','central', ...
              'FiniteDifferenceStepSize',[1e-2; 1e-10], ...  % (Hz step, seconds step)
              'StepTolerance',1e-14, ...
              'OptimalityTolerance',1e-12);

nMonteCarlo = 100;

N = 70;
M = 50;

c = 3e8;
fc = 30*1e9; % 30 GHz
lambda = c / fc;

deltaf = 200*1e3; % 200 kHz
Tcp = 1*1e-6; % 1 us
T = 1/deltaf + Tcp;

delta_tau = 1/(N*deltaf);
delta_nu = 1/(M*T);

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
nu_gt = 2*v_bis*cos(beta/2)/lambda;

k_tau = tau_gt/delta_tau;
l_nu = nu_gt/delta_nu;

% SNR_dbs = -30:60;
SNR_dbs = 0;
SNR_lins = db2pow(SNR_dbs);

MCRBs = zeros(4, 4, length(SNR_dbs), nMonteCarlo);
tau_errors = zeros(length(SNR_dbs), nMonteCarlo);

% pilot_periods = [10, 2, 2, 1; 5, 5, 1, 1];
pilot_period = [2; 5];

np = pilot_period(1);
mp = pilot_period(2);

tic
for SNR_idx = 1:length(SNR_dbs)
    SNR_db = SNR_dbs(SNR_idx);
    SNR_lin = SNR_lins(SNR_idx);
    % rng(0);
    sigma2 = 1;
    alpha_gt = sqrt(SNR_lin * sigma2) * exp(1j*rand*2*pi);
    
    symbols = randi(4, N, M)-1;
    Xdd = qammod(symbols, 4, UnitAveragePower=true);
    Xdd(N/2:N/2+17, M/2-3:M/2+3) = zeros(18, 7);
    Xdd(N/2, M/2) = 1;

    Xtf = isfft(Xdd);
    
    ns = (0:(N-1)).';
    ms = 0:(M-1);
    
    H = alpha_gt * exp(-1j*2*pi*deltaf*ns*tau_gt) * exp(1j*2*pi*T*ms*nu_gt);
    
    % Z = (randn(N, M) + 1j*randn(N, M)) * sqrt(sigma2/2);
    Z = 0;
    Ytf = H .* Xtf + Z;
    
    Ydd = sfft(Ytf);
end
toc

figure;
mesh(abs(Xdd));

figure;
mesh(abs(Ydd));