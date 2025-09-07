clear;
close all;
clc

load("sim7v5.mat");

figure; hold on;
plot(rcs2_dB, Pds(:, 1), LineWidth=2, Color="red")
plot(rcs2_dB, Pds(:, 2), LineWidth=2, Color="yellow")
plot(rcs2_dB, Pds(:, 3), LineWidth=2, Color="magenta")


clear;
load("sim2v5.mat");

plot(rcs2_dB, Pds, LineWidth=2, Color="green")

clear;
load("sim3v5.mat");

plot(rcs2_dB, Pds, LineWidth=2, Color="blue")


xlabel("Target RCS (dBsm)");
ylabel("Probability of Detection");
legend("Data Aided (Iteration-1)", "Data Aided (Iteration-1)", "Data Aided (Iteration-1)", "Genie", "Pilot Only");
grid on;
