function dopplers = getDopplerArray(delta_f, M, is_fftshifted)
    if nargin < 3
        is_fftshifted = 1;
    end
    dopplers = linspace(0, delta_f*(M-1)/M, M);
    if is_fftshifted
        dopplers = dopplers - dopplers(floor(length(dopplers)/2)+1);
    end
end