clear;
close all;
clc

load("results\sim11v7.mat");

colors = orderedcolors("gem");
figure; hold on;
plot(rcs2_dB, Pds(:, 1), "Color", colors(2, :), LineWidth=2);
plot(rcs2_dB, Pds(:, 2), "Color", colors(3, :), LineWidth=2);
plot(rcs2_dB, Pds(:, 3), "Color", colors(4, :), LineWidth=2);
legend("Iteration-1", "Iteration-2", "Iteration-3");
xlabel("Target RCS (dBsm)");
ylabel("Probability of Detection");
grid on;


% colors = orderedcolors("gem");
% figure; hold on;
% plot(pilot_ratio.*100, Pds, "Color", colors(2, :), LineWidth=2);
% xlabel("Pilot Percentage");
% ylabel("Probability of Detection");
% ylim([0, 1])
% grid on;
