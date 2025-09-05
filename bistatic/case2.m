clear;      % sim1: data aided, sim2 genie (100% pilot), sim3 pilot only (only first stage)
close all;
clc;

addpath('modules');  % Optional folder for functions

pilot_ratios = (2:10)./100;
Pds = zeros(size(pilot_ratios));

tic
for r_idx = 1:length(pilot_ratios)
    r_idx
    clearvars -except r_idx pilot_ratios Pds;
    rcs2_dB = 2;
    nIter = 1;
    monteCarlo = 500;
    is_genie = 1;
    is_data_only = 1;
    params = init_simulation_params(rcs2_dB, pilot_ratios(r_idx), nIter, monteCarlo, is_genie, is_data_only);
    v2struct(params);

    Pd = 0;
    for m = 1:monteCarlo
        main_bistatic_isac;
        Pd = Pd + is_ref_detected;
    end
    Pd = Pd / monteCarlo;
    Pds(r_idx) = Pd;
end
toc

figure;
plot(pilot_ratios, Pds);

save("sim5v5.mat")



clear;      % sim1: data aided, sim2 genie (100% pilot), sim3 pilot only (only first stage)
close all;
clc;

addpath('modules');  % Optional folder for functions

pilot_ratios = (2:20)./100;
Pds = zeros(size(pilot_ratios));

tic
for r_idx = 1:length(pilot_ratios)
    r_idx
    clearvars -except r_idx pilot_ratios Pds;
    rcs2_dB = 2;
    nIter = 1;
    monteCarlo = 500;
    is_genie = 0;
    is_data_only = 1;
    params = init_simulation_params(rcs2_dB, pilot_ratios(r_idx), nIter, monteCarlo, is_genie, is_data_only);
    v2struct(params);

    Pd = 0;
    for m = 1:monteCarlo
        main_bistatic_isac;
        Pd = Pd + is_ref_detected;
    end
    Pd = Pd / monteCarlo;
    Pds(r_idx) = Pd;
end
toc

figure;
plot(pilot_ratios, Pds);

save("sim6v5.mat")