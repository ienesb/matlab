function [delay_array, range_array] = getDelayArray(delta_f, N)
    c = 3e8;
    delay_array = (0:N-1)/ N / delta_f;
    range_array = delay_array * c/2;
end