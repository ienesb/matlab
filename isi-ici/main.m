clear;
close all;

M = 64;
N = 30;
df = 15e3;
T = 1/df;
fc = 5e9;

padLen = 0;
padType = 'NONE';

NT = 8;
NR = 8;





XDD = zeros(N, M);
FM = dftmtx(M)/sqrt(M);
FN = dftmtx(N)/sqrt(N);

Gtx = eye(N,N) ./ sqrt(T);

s = Gtx*XDD*FM';
s = s(:);