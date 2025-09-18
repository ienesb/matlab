function plotDDMap(DD, params, is_indexed)
    v2struct(params);
    [N, M] = size(DD);

    if nargin < 3
        is_indexed = 0;
    end

    if is_indexed
        delays = 0:(N-1);
        dopplers = 0:(M-1);
    else
        delays = delay_array;
        dopplers = doppler_array;
    end
    
    figure;
    imagesc(dopplers, delays, DD);
    colorbar;
    if ~is_indexed
        xlabel("Doppler (Hz)");
        ylabel("Delay (s)");
    end
end