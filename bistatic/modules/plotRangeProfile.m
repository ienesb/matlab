function plotRangeProfile(HDD, params, is_indexed)
    v2struct(params);

    if nargin < 3
        is_indexed = 0;
    end    

    if ~is_indexed
        delays = getDelayArray(Tsym, N_fft);
    else
        delays = 1:size(HDD, 1);
    end

    figure;
    plot(delays, HDD(:, 129))
end