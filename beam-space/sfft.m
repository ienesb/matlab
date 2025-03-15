function Xdd = sfft(Xtf, M)
    Xdd = fft(ifft(Xtf).').' * M;
end