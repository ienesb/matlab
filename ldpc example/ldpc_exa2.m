clear;
close all;
clc;

% load("C:\Users\İsmailEnesBülbül\Desktop\matlab\ldpc example\ldpc_H_90x100_dv9_dc10.mat")
load("C:\Users\İsmailEnesBülbül\Desktop\matlab\ldpc example\ldpc_H_90x100_systematic.mat")

H = sparse(logical(H));

pcmatrix = sparse(H);

cfgLDPCEnc = ldpcEncoderConfig(pcmatrix);
cfgLDPCDec = ldpcDecoderConfig(pcmatrix);

infoBits = randi([0, 1], 10, 1);

encodedBits = ldpcEncode(infoBits, cfgLDPCEnc);

M = 4;
maxnumiter = 10;
snr = [0.24 0.5 1 1.5 2];
numframes = 20000;

for ii = 1:length(snr)
    data = randi([0 1],cfgLDPCEnc.NumInformationBits,numframes,'int8');

    % Transmit and receive with LDPC coding
    encodedData = ldpcEncode(data,cfgLDPCEnc);
    modSignal = qammod(encodedData,M,InputType='bit');
    [rxsig, noisevar] = awgn(modSignal,snr(ii));
    llrOut = qamdemod(rxsig,M, ...
        OutputType='approxllr', ... % Or use 'llr' for exact but slower LLR calculation
        NoiseVariance=noisevar);
    rxbits = ldpcDecode(llrOut,cfgLDPCDec,maxnumiter);
    fprintf(['SNR = %2.1f\n   Coded: Error rate = %1.6f, ' ...
        'Number of errors = %d\n'], ...
        snr(ii),nnz(data~=rxbits)/numel(data),nnz(data~=rxbits));

    % Transmit and receive with no LDPC coding
    noCoding = qammod(data,M,InputType='bit');
    rxNoCoding = awgn(noCoding,snr(ii));
    rxBitsNoCoding = qamdemod(rxNoCoding,M,OutputType='bit');

    fprintf(['Noncoded: Error rate = %1.6f, ' ...
        'Number of errors = %d\n\n'], ...
        nnz(data~=rxBitsNoCoding)/numel(data),nnz(data~=rxBitsNoCoding))
end
