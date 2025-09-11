M = 16;       % Modulation order for 16QAM
nfft = 128;   % FFT length
cplen = 10;   % Cyclic prefix length
nSym = 1;     % Number of symbols per RE
nStreams = 1; % Number of transmit streams

dataIn = randi([0 M-1],nfft,nSym,nStreams);
qamSig = qammod(dataIn,M,'UnitAveragePower',true);
y1 = ofdmmod(qamSig,nfft,cplen);


temp = ifft(ifftshift((qamSig)));
temp = temp(:);
temp = [temp(end-cplen+1:end); temp];

isequal(y1, temp)

figure
stem(abs(y1));

figure
stem(abs(temp));