clear;
close all;
clc;

N = 10;
M = 5;

X = rand(N, M);

FN = dftmtx(N) / sqrt(N);
FM = dftmtx(M) / sqrt(M);

X1 = fft(X, N, 1) / sqrt(N);
X2 = ifft(X, N, 1) * sqrt(N);
X3 = fft(X, M, 2) / sqrt(M);
X4 = ifft(X, M, 2) * sqrt(M);


FNX = FN * X;
FNHX = FN' * X;
XFM = X * FM;
XFMH = X * FM';

sum(abs(X1(:) - FNX(:)))
sum(abs(X2(:) - FNHX(:)))
sum(abs(X3(:) - XFM(:)))
sum(abs(X4(:) - XFMH(:)))