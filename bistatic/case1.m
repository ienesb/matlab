clear;      % sim1: data aided, sim2 genie (100% pilot), sim3 pilot only (only first stage)
close all;
clc;

addpath('modules');  % Optional folder for functions


rcs2_dB = -5:5;
rcs2 = db2pow(rcs2_dB);

nIter = 3;

Pds = zeros(length(rcs2), nIter);

tic
for r_idx = 1:length(rcs2)
    r_idx
    clearvars -except r_idx rcs2 rcs2_dB Pds nIter;
    pilot_ratio = 0.05;
    monteCarlo = 500;
    is_genie = 0;
    is_data_only = 0;
    params = init_simulation_params(rcs2_dB(r_idx), pilot_ratio, nIter, monteCarlo, is_genie, is_data_only);
    v2struct(params);

    Pd = zeros(1, nIter);
    for m = 1:monteCarlo
        main_bistatic_isac;
        Pd = Pd + is_ref_detected_array;
    end
    Pd = Pd ./ monteCarlo;
    Pds(r_idx, :) = Pd;
end
toc

figure;
plot(rcs2_dB, Pds);

save("sim7v5.mat")