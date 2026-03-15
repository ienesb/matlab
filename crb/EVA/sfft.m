function Xdd = sfft(Xft) % freq -> delay, time -> doppler
    [N, M] = size(Xft);

    temp = sqrt(N) * ifft(Xft, [], 1);
    Xdd = 1/sqrt(M) * fft(temp, [], 2);
    % temp = ifft(Xft, [], 1);
    % Xdd = fft(temp, [], 2);
end
