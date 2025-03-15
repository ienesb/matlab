function Xtf = isfft(Xdd, M)
    Xtf = fft(ifft(Xdd.').' / M);
end