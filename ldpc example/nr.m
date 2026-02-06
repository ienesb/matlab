clear;
close all;
clc;

% Set up DL-SCH coding parameters
TBS = 13680;            % Transport block size, a positive integer
codeRate = 0.3;   % Target code rate, a real number between 0 and 1
rv = 0;                % Redundancy version, 0-3
modulation = 'QPSK';   % Modulation scheme, QPSK, 16QAM, 64QAM, 256QAM
nlayers = 1;           % Number of layers, 1-4 for a transport block
cbsInfo = nrDLSCHInfo(TBS,codeRate);
disp('DL-SCH coding parameters')
disp(cbsInfo)

switch modulation
    case 'QPSK'
        bitsPerSymbol = 2;
    case '16QAM'
        bitsPerSymbol = 4;
    case '64QAM'
        bitsPerSymbol = 6;
    case '256QAM'
        bitsPerSymbol = 8;
end

% Set up AWGN channel
EbNo = 1.25; % in dB
outlen = ceil(TBS/codeRate);
snrdB = convertSNR(EbNo,"ebno",...
    BitsPerSymbol=bitsPerSymbol,CodingRate=TBS/outlen);

% Random transport block data generation
in = randi([0 1],TBS,1,'int8');
% Transport block CRC attachment
tbIn = nrCRCEncode(in,cbsInfo.CRC);
% Code block segmentation and CRC attachment
cbsIn = nrCodeBlockSegmentLDPC(tbIn,cbsInfo.BGN); % 1: 8448, 2: 3840
% LDPC encoding
enc = nrLDPCEncode(cbsIn,cbsInfo.BGN);
% Rate matching and code block concatenation
chIn = nrRateMatchLDPC(enc,outlen,rv,modulation,nlayers);
% Symbol mapping
symOut = nrSymbolModulate(chIn,modulation);
symOutTrue = symOut;
% AWGN channel
[rxSig, noiseVar] = awgn(symOut,snrdB);
% Symbol demapping
rxllr = nrSymbolDemodulate(rxSig,modulation);
% Rate recovery
raterec = nrRateRecoverLDPC(rxllr,TBS,codeRate,rv,modulation,nlayers);
% LDPC decoding, with early termination and at most 12 iterations
decBits = nrLDPCDecode(raterec,cbsInfo.BGN,12);
% Code block desegmentation and CRC decoding
[blk,~] = nrCodeBlockDesegmentLDPC(decBits,cbsInfo.BGN,TBS+cbsInfo.L);
% Transport block CRC decoding
[out, err] = nrCRCDecode(blk,cbsInfo.CRC);
% Compare
blockError = any(out~=in)
sum(out~=in) / length(out)

% % Transport block CRC attachment
% tbIn = nrCRCEncode(out,cbsInfo.CRC);
% Code block segmentation and CRC attachment
cbsIn = nrCodeBlockSegmentLDPC(blk,cbsInfo.BGN);
% LDPC encoding
enc = nrLDPCEncode(cbsIn,cbsInfo.BGN);
% Rate matching and code block concatenation
chIn = nrRateMatchLDPC(enc,outlen,rv,modulation,nlayers);
% Symbol mapping
symOut = nrSymbolModulate(chIn,modulation);

sum(symOutTrue(:) ~= symOut) / length(symOut)


