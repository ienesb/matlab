clear;      % sim4: data aided, sim5 genie (100% pilot), sim6 pilot only (only first stage)
close all;
clc;

addpath('modules');  % Optional folder for functions


pilot_ratios = (2:10) ./ 100;
Pds = zeros(size(pilot_ratios));

tic
for r_idx = 1:length(pilot_ratios)
    r_idx
    clearvars -except r_idx pilot_ratios Pds;
    params = init_simulation_params(pilot_ratios(r_idx));
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

save("sim6.mat")