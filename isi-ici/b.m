M = 64;          % number of subcarriers
N = 30;          % number of subsymbols/frame
df = 15e3;       % make this the frequency bin spacing of LTE
T = 1/df;
fc = 5e9;        % carrier frequency in Hz
padLen = 0;     % make this larger than the channel delay spread channel in samples
padType = 'NONE';  % this example requires ZP for ISI mitigation

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

