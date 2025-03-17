M = 16;
N = 32;
df = 120e3;
T = 1/df;
fc = 28e9;

tau = 0.1;
nu = 0.1;

bN = exp(-1j*2*pi*df*tau*(0:N-1)).';
bISI = exp(-1j*2*pi*df*tau*(0:M-1)/M).';

b_tau = kron(bN, bISI);

B = diag((b_tau));


cM = exp(1j*2*pi*T*nu*(0:M-1)).';
cICI = exp(1j*2*pi*T*nu*(0:N-1)/N).';

c_nu = kron(cM, cICI);

C = diag((c_nu));

