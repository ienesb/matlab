function plotDDMap(DD, params, is_indexed, is_fftshifted)
    v2struct(params);
    [N, M] = size(DD);

    if nargin < 4
        is_fftshifted = 1;
        if nargin < 3
            is_indexed = 0;
        end
    end

    delays = getDelayArray(Tsym, N);
    dopplers = getDopplerArray(delta_f, M, is_fftshifted);

    if is_indexed
        delays = 0:(N-1);
        dopplers = 0:(M-1);
    end
    
    figure;
    imagesc(dopplers, delays, DD);
    colorbar;
    if ~is_indexed
        xlabel("Doppler (Hz)");
        ylabel("Delay (s)");
    end
end