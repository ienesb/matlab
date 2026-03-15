function Xtf = isfft(Xdd)
    [N, M] = size(Xdd);

    temp = 1/sqrt(N) * fft(Xdd, [], 1);
    Xtf = sqrt(M) * ifft(temp, [], 2);
    % temp = fft(Xdd, [], 1);
    % Xtf = ifft(temp, [], 2);
end