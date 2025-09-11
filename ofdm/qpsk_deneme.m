M = 4;

symLen = 100000;
bitLen = symLen * log2(4);

bits = randi([0, 1], bitLen, 1);

symbols = bits(1:2:end)*2 + bits(2:2:end);

qpskSymbols = qammod(symbols, M, UnitAveragePower=true);

snr_db = -10:10;
snr_lin = db2pow(snr_db);
EbN0 = snr_lin/2;
EbN0_db = pow2db(EbN0);

bers = zeros(size(snr_db));
sers = zeros(size(snr_db));

for idx = 1:length(snr_db)
    snr = db2pow(snr_db(idx));
    sigma2 = (1 / snr) / 2;
    noise = (randn(size(qpskSymbols)) + 1j*randn(size(qpskSymbols)))*sqrt(sigma2);
    receivedSignal = qpskSymbols + noise;
    receivedSymbols = qamdemod(receivedSignal, M);
    receivedBits = zeros(length(receivedSymbols)*2, 1);
    receivedBits(1:2:end) = floor(receivedSymbols/2);
    receivedBits(2:2:end) = mod(receivedSymbols, 2);
    
    bers(idx) = sum(receivedBits ~= bits)/bitLen;
    sers(idx) = sum(receivedSymbols ~= symbols)/symLen;
end

[berTheoretical, serTheoretical] = berawgn(EbN0_db, "psk", M, "nondiff");

figure;
semilogy(snr_db, bers); hold on;
semilogy(snr_db, berTheoretical);
xlabel("SNR [dB]");
ylabel("BER");
legend("Simulation", "Theoretical");
grid on;

figure;
semilogy(snr_db, sers); hold on;
semilogy(snr_db, serTheoretical);
xlabel("SNR [dB]");
ylabel("SER");
legend("Simulation", "Theoretical");
grid on;
