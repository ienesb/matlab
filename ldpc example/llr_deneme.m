clear;
close all;

noiseVar = 1;

r = -1:0.0001:1;
llr = log(exp(-(r - 1/sqrt(2)).^2./noiseVar) ./ exp(-(r + 1/sqrt(2)).^2./noiseVar));

figure;
plot(r, abs(llr))