clear;
close all;
clc;

ofdm_ber;
figure(100);
semilogy(SNR_dbs, bers, "LineStyle", '-', 'LineWidth', 2);

otfs_ber;
figure(100);
hold on;
semilogy(SNR_dbs, bers, "LineStyle", '-', 'LineWidth', 2);
legend("OFDM", "OTFS");
grid on;
xlabel("SNR [dB]");
ylabel("BER");
title("BER vs SNR");
theme(gcf, "light");