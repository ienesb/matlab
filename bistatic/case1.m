clear;      % sim1: data aided, sim2 genie (100% pilot), sim3 pilot only (only first stage)
close all;
clc;

addpath('modules');  % Optional folder for functions


rcs2_dB = -10:10;
rcs2 = db2pow(rcs2_dB);
Pds = zeros(size(rcs2));

tic
for r_idx = 1:length(rcs2)
    r_idx
    clearvars -except r_idx rcs2 rcs2_dB Pds;
    params = init_simulation_params(rcs2(r_idx));
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

save("sim3.mat")