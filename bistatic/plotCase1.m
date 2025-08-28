clear;
close all;
clc

load("sim1.mat");

figure; hold on;
plot(rcs2_dB, Pds)

clear;
load("sim2.mat");

plot(rcs2_dB, Pds)

clear;
load("sim3.mat");

plot(rcs2_dB, Pds)

xlabel("Target RCS (dBsm)");
ylabel("Probability of Detection");
legend("Data Aided", "Genie", "Pilot Only");
grid on;
