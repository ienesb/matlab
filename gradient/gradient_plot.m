clear;
close all;
clc;

%% Parameters
is_load = 1;

N = 400;
M = 60;

deltaf = 120000;
Tsym = 1/deltaf*1.07; % 8.9167e-06;

% sigmaY2 = 0;
sigmaY2 = 1.2057e-12; % 0.01; % 1.2057e-12;

sigmaX2 = 1;
% sigmaX2 = 0;

%% Target Parameters
taus = [0.1328, 0.3995] * 1e-6;
nus =   [375.8984, 377.3253];
% gains = [5.3963e-07 + 1.9816e-24i, 5.4666e-08 - 7.4480e-08i];
gains = [5.3963e-07, 5.4666e-08];
% gains = [5.3963e-07 0];
% K = length(gains);
K = 2;

%% Ground Truth Generation
H = zeros(N, M);
for idx = 1:K
    b_tau = getb(taus(idx), deltaf, N);
    c_nu = getc(nus(idx), Tsym, M);
    H = H + gains(idx) * (b_tau * c_nu.');
end

X = qammod(randi([0, 3], N, M), 4, UnitAveragePower=true);
% X = ones(N, M);
Y = X .* H;
Y = Y + sqrt(sigmaY2/2) * (randn(N, M) + 1j * randn(N, M));


%% Error Curve Generate/Load
if ~is_load
    tau_trials = taus(1) + (-100:0.5:200) ./ 100 .* 0.2e-6;
    nu_trials = nus(1) + (-100:0.5:100) ./ 100 .* 1000;
    
    errors = zeros(length(tau_trials), length(nu_trials));
    
    for tau_idx = 1:length(tau_trials)
        for nu_idx = 1:length(nu_trials)
            tau = tau_trials(tau_idx);
            nu = nu_trials(nu_idx);

            errors(tau_idx, nu_idx) = gradient_error(tau, nu, gains(1), X_hat, Y);
        end
    end
    save("errors.mat", "errors", "tau_trials", "nu_trials")
else
    load("errors.mat");
end

%% Initial Estimation
% grids
[nu_grid, tau_grid] = meshgrid(nu_trials, tau_trials);   % NU: T x N, TAU: T x N

idx = randi([1, numel(nu_grid)], 1, 1);
tau_hat = taus(1); % taus(1) + 1.5e-8; % tau_grid(idx); % taus(1) + 1e-8;
nu_hat = nus(1); % nus(1) + 400; % nu_grid(idx); % nus(1) + 0.01;
alpha_hat = gains(1);

b_tau = getb(tau_hat, deltaf, N);
c_nu = getc(nu_hat, Tsym, M);
H_hat = alpha_hat * (b_tau * c_nu.');
% H_hat = H;

X_hat = X;
% X_hat = (Y .* conj(H_hat)) ./ (abs(H_hat).^2 + sigmaY2);
% X_hat = (Y .* conj(H_hat)) ./ (abs(H_hat).^2 + sigmaY2) + sqrt(sigmaX2/2) * (randn(size(X)) + 1j*randn(size(X)));
X_hat = sqrt(sigmaX2/2) * (randn(size(X)) + 1j*randn(size(X)));
% X_hat = X + sqrt(sigmaX2/2) * (randn(size(X)) + 1j*randn(size(X)));
[RMSE, ber] = getSer(X, X_hat);

X_hat_LMMSE = (Y .* conj(H_hat)) ./ (abs(H_hat).^2 + sigmaY2);
[RMSE_LMMSE, ber_LMMSE] = getSer(X, X_hat_LMMSE);

figure;
scatter(real(X_hat(:)), imag(X_hat(:)), "x", "red");

J = gradient_error(tau_hat, nu_hat, alpha_hat, X_hat, Y);

%% Epochs
for e = 1:1000
    [tau_hat(e+1), nu_hat(e+1), X_hat(:, :, e+1)] = gradient_step(tau_hat(e), nu_hat(e), alpha_hat, X_hat(:, :, e), Y);
    J(e+1) = gradient_error(tau_hat(e+1), nu_hat(e+1), alpha_hat, X_hat(:, :, e+1), Y);
    [RMSE(e+1), ber(e+1)] = getSer(X, X_hat(:, : ,e+1));
end


%% Results
% plot
figure;
mesh(nu_grid, tau_grid, errors); hold on;
plot3(nu_hat, tau_hat, J, 'r', 'LineWidth', 2);
xlabel('\nu (Hz)');
ylabel('\tau (s)');
zlabel('Error');
title('Error over (\tau,\nu)');
grid on;

figure;
plot(J);
title("J");

figure;
plot(RMSE, 'LineWidth', 2); hold on;
plot(RMSE_LMMSE*ones(1, length(RMSE)), 'LineWidth', 2);
xlabel("Iterations");
ylabel("RMSE");
legend("Gradient", "LMMSE");
grid on;

figure;
plot(ber, 'LineWidth', 2); hold on;
plot(ber_LMMSE*ones(1, length(ber)), 'LineWidth', 2);
xlabel("Iterations");
ylabel("BER");
legend("Gradient", "LMMSE");
grid on;

figure; hold on;
plot(tau_hat);
plot(ones(size(tau_hat)) * taus(1));
title("tau hat");

figure; hold on;
plot(nu_hat);
plot(ones(size(nu_hat)) * nus(1));
title("nu hat");

figure;
scatter(real(X_hat(:, :, end)), imag(X_hat(:, :, end)), "x", "red");

c = linspace(1,10,length(squeeze(real(X_hat(1, 1, :)))));
figure;
scatter(squeeze(real(X_hat(20, 30, :))), squeeze(imag(X_hat(20, 30, :))), [], c, "x");