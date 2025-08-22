function delays = getDelayArray(Tsym, N)
    delays = linspace(0, Tsym*(N-1)/N, N);
end