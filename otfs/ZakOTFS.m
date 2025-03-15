clear;
close all;

M = 64;          % number of subcarriers
N = 30;          % number of subsymbols/frame
df = 15e3;       % make this the frequency bin spacing of LTE
fc = 5e9;        % carrier frequency in Hz
padLen = 10;     % make this larger than the channel delay spread channel in samples
padType = 'ZP';  % this example requires ZP for ISI mitigation

% Configure paths
chanParams.pathDelays      = [0  5   8  ]; % number of samples that path is delayed
chanParams.pathGains       = [1  0.7 0.5]; % complex path gain
chanParams.pathDopplers    = [0 -3/100   5/100  ]; % Doppler index as a multiple of fsamp/MN

fsamp = M*df;            % sampling frequency at the Nyquist rate
Meff = M + padLen;       % number of samples per OTFS subsymbol
numSamps = Meff * N;     % number of samples per OTFS symbol
T = ((M+padLen)/(M*df)); % symbol time (seconds)

% Calculate the actual Doppler frequencies from the Doppler indices
chanParams.pathDopplerFreqs = chanParams.pathDopplers * 1/(N*T); % Hz

berOTFSs = [];
berOFDMs = [];

SNRdBs = 0:1:10;
for SNRdB = SNRdBs

    % Pilot generation and grid population
    pilotBin = floor(N/2)+1;
    Pdd = zeros(M,N);
    Pdd(1,pilotBin) = exp(1i*pi/4); % populate just one bin to see the effect through the channel
    
    % OTFS modulation
    % txOut = helperOTFSmod(Pdd,padLen,padType);
    txOut = ZakOTFSmod(Pdd, M, N, padLen);
   
    % Send the OTFS modulated signal through the channel
    dopplerOut = dopplerChannel(txOut,fsamp,chanParams);
    
    % Add white Gaussian noise
    Es = mean(abs(pskmod(0:3,4,pi/4).^ 2));
    n0 = Es/(10^(SNRdB/10));
    chOut = awgn(dopplerOut,SNRdB,'measured');
    
    % Get a sample window
    rxIn = chOut(1:numSamps);
    
    % OTFS demodulation
    % Ydd = helperOTFSdemod(rxIn,M,padLen,0,padType);
    Ydd = ZakOTFSdemod(rxIn, M, N, padLen);
    
    % LMMSE channel estimate in the delay-Doppler domain
    Hdd = Ydd * conj(Pdd(1,pilotBin)) / (abs(Pdd(1,pilotBin))^2 + n0);
       
    [lp,vp] = find(abs(Hdd) >= 0.05);
    chanEst.pathGains = diag(Hdd(lp,vp));   % get path gains
    chanEst.pathDelays = lp - 1;            % get delay indices
    chanEst.pathDopplers = vp - pilotBin;   % get Doppler indices

    % Data generation
    Xgrid = zeros(M,N);
    Xdata = randi([0,1],2*M,N);
    Xgrid(1:M,:) = pskmod(Xdata,4,pi/4,InputType="bit");
    
    % OTFS modulation
    % txOut = helperOTFSmod(Xgrid,padLen,padType);
    txOut = ZakOTFSmod(Xgrid, M, N, padLen);
    
    % Add channel and noise
    dopplerOut = dopplerChannel(txOut,fsamp,chanParams);
    chOut = awgn(dopplerOut,SNRdB,'measured');
    
    % Form G matrix using channel estimates
    G = getG(M,N,chanEst,padLen,padType);
    
    rxWindow = chOut(1:numSamps);
    y_otfs = ((G'*G)+n0*eye(Meff*N)) \ (G'*rxWindow); % LMMSE
    
    % Xhat_otfs = helperOTFSdemod(y_otfs,M,padLen,0,padType); % OTFS demodulation
    Xhat_otfs = ZakOTFSdemod(y_otfs, M, N, padLen);
    
    XhatDataOTFS = pskdemod(Xhat_otfs,4,pi/4,OutputType="bit",OutputDataType="logical");
    [~,berOTFS] = biterr(Xdata,XhatDataOTFS);
    fprintf('OTFS BER with LMMSE equalization = %3.3e\n', berOTFS);
    berOTFSs(end+1) = berOTFS;
end

% figure; 
semilogy(SNRdBs, berOTFSs);
legend('OFDM', 'OTFS', 'ZakOTFS');