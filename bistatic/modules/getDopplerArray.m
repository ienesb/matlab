function [doppler_array, velocity_array] = getDopplerArray(Tsym, M, lambda)
    doppler_array = (-M/2:M/2-1)/ M / Tsym;
    velocity_array = doppler_array * lambda/2;
end