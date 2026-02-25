%% crb_pilot_only.m -> CRB for pilot only with different pilot ratios
clear;
% close all;
clc;

plot_count = 1;

N = 70;
M = 50;

deltaf = 200*1e3;   % 200 kHz
Tcp = 1e-6;         % 1 us
T = 1 / deltaf;
Ts = T + Tcp;

fc = 30e9;          % 30 GHz
c = 3e8;
lambda = c / fc;

beta = atan(130/90) - atan(50/90);

figure(1);
% figure(2);

for p = [10, 2, 2, 1; 5, 5, 1, 1]
    np = p(1);
    mp = p(2);

    K = ceil((N - 1) / np);
    L = ceil((M - 1) / mp);
    
    absP = (K + 1) * (L + 1);
    pilot_percentage = absP/(N*M);
    
    SNR_db = -30:30;
    
    SNR = db2pow(SNR_db);
    
    CRB_range = 12 / (K * (K + 2) * absP * np^2) * c^2 ./ (8 * pi^2 * SNR * deltaf^2);
    CRB_range = sqrt(CRB_range);

    CRB_velocity = 12 / (L * (L + 2) * absP * mp^2) * lambda^2 ./ (32 * pi^2 * SNR * Ts^2 * cos(beta/2)^2);
    CRB_velocity = sqrt(CRB_velocity);
    
    if plot_count == 1
        figure(1);
        semilogy(SNR_db, CRB_range, "Color", "blue", "LineWidth", 2);
        figure(2);
        semilogy(SNR_db, CRB_velocity, "Color", "blue", "LineWidth", 2);
    elseif plot_count == 2
        figure(1); hold on;
        semilogy(SNR_db, CRB_range, "Color", "red", "LineWidth", 2);
        figure(2); hold on;
        semilogy(SNR_db, CRB_velocity, "Color", "red", "LineWidth", 2);
    elseif plot_count == 3
        figure(1); hold on;
        semilogy(SNR_db, CRB_range, "Color", "green", "LineWidth", 2);
        figure(2); hold on;
        semilogy(SNR_db, CRB_velocity, "Color", "green", "LineWidth", 2);
    elseif plot_count == 4
        figure(1); hold on;
        semilogy(SNR_db, CRB_range, "Color", "#ff5500", "LineWidth", 2);
        figure(2); hold on;
        semilogy(SNR_db, CRB_velocity, "Color", "#ff5500", "LineWidth", 2);
    end
    plot_count = plot_count+1;
end

figure(1);
xlabel('SNR (dB)');
ylabel('CRB Range');
title('Cramér-Rao Bound vs SNR');
legend("0.02", "0.1", "0.5", "1");
grid on;

figure(2);
xlabel('SNR (dB)');
ylabel('CRB Velocity');
title('Cramér-Rao Bound vs SNR');
legend("0.02", "0.1", "0.5", "1");
grid on;
