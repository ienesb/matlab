fs = 48000; % 48 kHz sampling rate
T = 3; % 3 s total duration
f1 = 2000; % 1kHz
f2 = 2000; % 2kHz

n = 1:(fs*T) - 1; % sample indexes
n = n';
ts = n./fs; % time array

x1 = cos(2*pi*f1*n/fs);
x2 = frequencyOffset(x1, fs, -100);

X1 = fftshift(fft(x1));

X2 = fftshift(fft(x2));

freqs = n./T;
freqs = freqs - fs/2;  

figure; hold on;
plot(freqs, abs(X1));
plot(freqs, abs(X2));

figure; hold on;
plot(ts, x1);
plot(ts, x2);
 
