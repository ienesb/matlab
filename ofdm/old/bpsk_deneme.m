M = 2;

bitLen = 100000;

data = randi([0, M-1], bitLen, 1);

bpskSymbols = qammod(data, M);

snr_db = -10:1:10;

bers = zeros(size(snr_db));

for idx = 1:length(snr_db)
    snr = db2pow(snr_db(idx));
    sigma2 = (1 / snr) / 2;
    noise = randn(size(bpskSymbols))*sqrt(sigma2);
    receivedSignal = bpskSymbols + noise;
    % receivedSignal = awgn(bpskSymbols, snr);
    receivedSymbols = qamdemod(receivedSignal, M);
    bers(idx) = sum(receivedSymbols~= data)/bitLen;
end

berTheoretical = berawgn(snr_db, "psk", M, "nondiff");

figure;
semilogy(snr_db, bers); hold on;
semilogy(snr_db, berTheoretical);
xlabel("SNR [dB]");
ylabel("BER");
legend("Simulation", "Theoretical");
grid on;