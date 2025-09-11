M = 4;                 % Modulation alphabet
k = log2(M);           % Bits/symbol
numSC = 128;           % Number of OFDM subcarriers
cpLen = 32;            % OFDM cyclic prefix length
maxBitErrors = 100;    % Maximum number of bit errors
maxNumBits = 1e7;      % Maximum number of bits transmitted


ofdmMod = comm.OFDMModulator( ...
    FFTLength=numSC, ...
    CyclicPrefixLength=cpLen, ...
    Windowing=true, ...
    WindowLength=16);
ofdmDemod = comm.OFDMDemodulator( ...
    FFTLength=numSC, ...
    CyclicPrefixLength=cpLen);


errorRate = comm.ErrorRate(ResetInputPort=true);


ofdmDims = info(ofdmMod)


numDC = ofdmDims.DataInputSize(1)


frameSize = [k*numDC 1];


EbNoVec = (0:10)';
snrVec = EbNoVec + 10*log10(k) + 10*log10(numDC/numSC);


berVec = zeros(length(EbNoVec),3);
errorStats = zeros(1,3);


for m = 1:length(EbNoVec)
    snr = snrVec(m);
    
   while errorStats(2) <= maxBitErrors && errorStats(3) <= maxNumBits
       dataIn = randi([0,1],frameSize);                
       qpskTx = pskmod(dataIn,M,InputType="bit");      
       txSig = ofdmMod(qpskTx);                        
       powerDB = 10*log10(var(txSig));                 
       noiseVar = 10.^(0.1*(powerDB-snr));
       noise = sqrt(noiseVar/2)*complex(randn(size(txSig)), ...
           randn(size(txSig)));
       rxSig = txSig + noise;
       qpskRx = ofdmDemod(rxSig);                      
       dataOut = pskdemod(qpskRx,M,OutputType="bit");  
       errorStats = errorRate(dataIn,dataOut,0);       
   end
    
    berVec(m,:) = errorStats;                          
    errorStats = errorRate(dataIn,dataOut,1);          
end


berTheory = berawgn(EbNoVec,'psk',M,'nondiff');


figure
semilogy(EbNoVec,berVec(:,1),'*')
hold on
semilogy(EbNoVec,berTheory)
legend('Simulation','Theory','Location','Best')
xlabel('Eb/No (dB)')
ylabel('Bit Error Rate')
grid on
hold off


