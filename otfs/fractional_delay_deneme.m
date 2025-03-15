FD = 0.0003;
N = 10;
h = designFracDelayFIR(FD,N);

T = 1/44000;
f = 1000;

% Create an FIR filter object
fdfir = dsp.FIRFilter(h); 

% Generate some input
n = (1:100)';
t = n.*T;
x = cos(2*pi*f*T*n);

% Filter the input signal
y = fdfir(x);

figure; hold on;
plot(t, x);
plot(t, y);
legend('Original Sequence', 'Filter Output')
title('Raw Filter Output v.s. Input Sequence')

