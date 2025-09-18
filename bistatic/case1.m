clear;      % sim1: data aided, sim2 genie (100% pilot), sim3 pilot only (only first stage)
% close all;
clc;

addpath('modules');

rcs2_dB = -5:0.5:5;
rcs2 = db2pow(rcs2_dB);

nIter = 1;
pilot_ratio = 0.05;
monteCarlo = 500;
is_genie = 0;

Pds = zeros(length(rcs2), nIter);
sers = zeros(length(rcs2), 1);
rmse_dBs = zeros(length(rcs2), 1);

tic
for idx = 1:length(rcs2)    
    params = init_simulation_params(rcs2_dB(idx), pilot_ratio, nIter, monteCarlo);

    Pd = zeros(nIter, 1);
    ser = 0;
    rmse = 0;
    parfor m = 1:monteCarlo % parfor
        % is_ref_detected = pilot_only(params, is_genie);
        [is_ref_detected, s, r] = data_aided_sim(params);
        Pd = Pd + is_ref_detected;
        ser = ser + s;
        rmse = rmse + r;
    end

    Pd = Pd ./ monteCarlo;
    ser = ser ./ monteCarlo;
    rmse = rmse ./ monteCarlo;
    Pds(idx, :) = Pd;
    sers(idx) = ser;
    rmse_dBs(idx) = pow2db(rmse);
end
toc

colors = orderedcolors("gem");
figure; hold on;
plot(rcs2_dB, Pds(:, 1), "Color", colors(2, :), LineWidth=2);
% plot(rcs2_dB, Pds(:, 2), "Color", colors(3, :), LineWidth=2);
% plot(rcs2_dB, Pds(:, 3), "Color", colors(4, :), LineWidth=2);
xlabel("Target RCS (dBsm)");
ylabel("Probability of Detection");
% legend("Iteration-1", "Iteration-2", "Iteration-3");
grid on;
% 
% figure;
% plot(rcs2_dB, sers, LineWidth=2);
% grid on;
% xlabel("Target RCS (dBsm)");
% ylabel("Symbol Error Rate");
% 
% figure;
% plot(rcs2_dB, rmse_dBs);

% save("results/sim1v6hard.mat")
