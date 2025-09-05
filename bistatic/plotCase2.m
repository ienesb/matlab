clear;
close all;
clc

load("sim4v5.mat");

figure(1); hold on;
plot(pilot_ratios, Pds)

figure(2); hold on;
plotRangeProfile(HDD, params, 0);


clear;
load("sim5v5.mat");

figure(1); hold on;
plot(pilot_ratios, Pds)
figure(2); hold on;
plotRangeProfile(HDD, params, 0);


clear;
load("sim6v5.mat");

figure(1); hold on;
plot(pilot_ratios, Pds)
figure(2); hold on;
plotRangeProfile(HDD, params, 0);


figure(1);
xlabel("Pilot Ratio (%)");
ylabel("Probability of Detection");
legend("Data Aided", "Genie", "Pilot Only");
grid on;

figure(2);
xlabel("Delay");
legend("Data Aided", "Genie", "Pilot Only");
grid on;
