clear;

m  = 3e5;
n  = 1e4;

tic
piEst = computePi(m,n);
toc

pool = parpool("Threads");

tic
piEstParfor = computePiParfor(m,n);
toc

delete(gcp('nocreate'))
% timeSerial = timeit(@() computePi(m,n));