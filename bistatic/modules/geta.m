function a = geta(theta, lambda, NT)
    d = lambda/2;
    a = exp(1i * 2 * pi / lambda * d * (0:NT-1).' * sind(theta));
end