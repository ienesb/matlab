clear;
close all;
clc;

%% Parameters
N = 400;

deltaf = 120000;
% deltaf = 1;
Tsym = 1/deltaf*1.07; % 8.9167e-06;

%% Target Parameters
tau_gt = 0.1328 * 1e-6;
ns = 0:N-1; % Define the sample indices

gt = exp(-1j*2*pi*ns*deltaf*tau_gt);

tau_trials = tau_gt + (-100:0.5:100) ./ 100 .* 0.2e-6;

errors = zeros(length(tau_trials), 1);
errors_log = zeros(length(tau_trials), 1);

for tau_idx = 1:length(tau_trials)
    tau_hat = tau_trials(tau_idx);

    est = exp(-1j*2*pi*ns*deltaf*tau_hat);

    errors(tau_idx) = sum(abs(gt - est).^2);
    errors_log(tau_idx) = sum(abs(log(gt) - log(est)).^2);
    % figure
    % plot(abs(log(gt)).^2)
end

max_idx = islocalmax(errors);
min_idx = islocalmin(errors);

figure; hold on;
plot(tau_trials, errors);
plot(tau_trials(max_idx), errors(max_idx), "*");
plot(tau_trials(min_idx), errors(min_idx), "*");
% plot(tau_trials, errors_log);
grid on;

in = -1j*2*pi*deltaf*tau_trials;
myexp = exp(in);
logexp = log(myexp);

figure;
plot(abs(logexp));

figure;
plot(abs(in));

figure
plot(angle(gt))


