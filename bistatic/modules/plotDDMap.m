function plotDDMap(DD, N, M, delta_f)
    Tsym = 1 / delta_f;
    delays = linspace(0, Tsym*(N-1)/N, N) .* 1e6;
    dopplers = linspace(0, delta_f*(M-1)/M, M);

    figure;
    imagesc(dopplers, delays, DD);
    colorbar;
    xlabel("Doppler (Hz)");
    ylabel("Delay (\mus)");
end