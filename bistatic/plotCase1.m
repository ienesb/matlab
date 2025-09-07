clear;
close all;
clc

load("sim1v5.mat");

figure; hold on;
plot(rcs2_dB, Pds, LineWidth=2, Color="magenta")

clear;
load("sim2v5.mat");

plot(rcs2_dB, Pds, LineWidth=2, Color="green")

clear;
load("sim3v5.mat");

plot(rcs2_dB, Pds, LineWidth=2, Color="blue")

xlabel("Target RCS (dBsm)");
ylabel("Probability of Detection");
legend("Data Aided", "Genie", "Pilot Only");
grid on;
