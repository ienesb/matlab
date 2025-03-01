tau = 1;
nu = 1;

T = 2;
dt = 1/T;

N = 5;
M = 6;

ns = (1:N)';
nsp = ns';

arg1 = (ns - nsp).*T - tau;

ms = (1:M)';
msp = ms';

arg2 = (ms - msp).*dt - nu;

arg = zeros(N,N,M,M);

