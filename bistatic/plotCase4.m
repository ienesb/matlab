clear;
close all;
clc

load("sim8v5.mat");

figure; hold on;

delays = getDelayArray(Tsym, N_fft);
delays = delays .* 1e6;
HDD_db = pow2db(HDD);
range_profile = HDD_db(:, 129); 
range_profile = range_profile - max(range_profile);
plot(delays, range_profile, LineWidth=2, Color="red")


clear;
load("sim9v5.mat");

delays = getDelayArray(Tsym, N_fft);
delays = delays .* 1e6;
HDD_db = pow2db(HDD);
range_profile = HDD_db(:, 129); 
range_profile = range_profile - max(range_profile);
plot(delays, range_profile, LineWidth=2, Color="green")


clear;
load("sim10v5.mat");

delays = getDelayArray(Tsym, N_fft);
delays = delays .* 1e6;
HDD_db = pow2db(HDD);
range_profile = HDD_db(:, 129); 
range_profile = range_profile - max(range_profile);
plot(delays, range_profile, LineWidth=2, Color="blue")

xline(taus.*1e6, LineWidth=2);

xlim([0, taus(2)*2*1e6])
xlabel("Delay (\mus)");
ylabel("Magnitude [dB]");
legend("Data Aided", "Genie", "Pilot Only", "Targets");
grid on;

