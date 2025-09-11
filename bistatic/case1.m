clear;      % sim1: data aided, sim2 genie (100% pilot), sim3 pilot only (only first stage)
% close all;
clc;

addpath('modules');

rcs2_dB = -5:5;
rcs2 = db2pow(rcs2_dB);

nIter = 3;
pilot_ratio = 0.05;
monteCarlo = 100;
is_genie = 0;
is_data_only = 0;

Pds = zeros(length(rcs2), nIter);

tic
for r_idx = 1:length(rcs2)    
    params = init_simulation_params(rcs2_dB(r_idx), pilot_ratio, nIter, monteCarlo, is_genie, is_data_only);

    Pd = zeros(nIter, 1);
    parfor m = 1:monteCarlo
        % is_ref_detected = genie_sim(params);
        is_ref_detected = data_aided_sim(params);
        Pd = Pd + is_ref_detected;
    end

    Pd = Pd ./ monteCarlo;
    Pds(r_idx, :) = Pd;
end
toc

figure;
plot(rcs2_dB, Pds);

% save("sim.mat")