clear;
% close all;

M = 64;          % number of subcarriers
N = 30;          % number of subsymbols/frame
df = 15e3;       % make this the frequency bin spacing of LTE
fc = 5e9;        % carrier frequency in Hz
padLen = 10;     % make this larger than the channel delay spread channel in samples
padType = 'ZP';  % this example requires ZP for ISI mitigation
  
fsamp = M*df;            % sampling frequency at the Nyquist rate
Meff = M + padLen;       % number of samples per OTFS subsymbol
numSamps = Meff * N;     % number of samples per OTFS symbol
T = ((M+padLen)/(M*df)); % symbol time (seconds)

  
% Calculate the actual Doppler frequencies from the Doppler indices
% chanParams.pathDopplerFreqs = chanParams.pathDopplers * 1/(N*T); % Hz
% for k = 1:length(chanParams.pathDelays)
%     fprintf('Scatterer %d\n',k);
%     fprintf('\tDelay = %5.2f us\n', 1e6*chanParams.pathDelays(k)/(Meff*df));
%     fprintf('\tRelative Doppler shift = %5.0f Hz (%5.0f km/h)\n', ...
%         chanParams.pathDopplerFreqs(k), (physconst('LightSpeed')*chanParams.pathDopplerFreqs(k)/fc)*(3600/1000));
% end

% Pilot generation and grid population
pilotBin = floor(N/2)+1;
Pdd = zeros(M,N);
Pdd(1,pilotBin) = exp(1i*pi/4); % populate just one bin to see the effect through the channel

% figure;

for c=0:0.25:1

% Configure paths
    chanParams.pathDelays      = [0  5   8  ]; % number of samples that path is delayed
    chanParams.pathGains       = [1  0.7*c 0.5*c]; % complex path gain
    chanParams.pathDopplers    = [0 -3   5  ]; % Doppler index as a multiple of fsamp/MN
    chanParams.pathDopplerFreqs = chanParams.pathDopplers * 1/(N*T); % Hz
    berOTFSs = [];
    
    SNRdBs = 0:5:40;
    for SNRdB = SNRdBs
        ber = 0;
        for i=1:10
            % Channel Estimation
            % Transmit pilots over all subcarriers and symbols to sound the channel
            txOut = helperOTFSmod(Pdd,padLen,padType);       % transmit pilots over the entire grid
            dopplerOut = dopplerChannel(txOut,fsamp,chanParams);    % send through channel
            Es = mean(abs(pskmod(0:3,4,pi/4).^ 2));
            n0 = Es/(10^(SNRdB/10));
            chOut = awgn(dopplerOut,SNRdB,'measured');              % add noise
            rxIn = chOut(1:numSamps);
            Ydd = helperOTFSdemod(rxIn,M,padLen,0,padType);      % demodulate
            Hdd = Ydd * conj(Pdd(1,pilotBin)) / (abs(Pdd(1,pilotBin))^2 + n0); % LMMSE channel estimate
            [lp,vp] = find(abs(Hdd) >= 0.05);
            chanEst.pathGains = diag(Hdd(lp,vp));   % get path gains
            chanEst.pathDelays = lp - 1;            % get delay indices
            chanEst.pathDopplers = vp - pilotBin;   % get Doppler indices

            % Data generation
            Xgrid = zeros(M,N);
            Xdata = randi([0,1],2*M,N);
            Xgrid(1:M,:) = pskmod(Xdata,4,pi/4,InputType="bit");
            
            % Transmit data over the same channel and use channel estimates to equalize
            txOut = helperOTFSmod(Xgrid,padLen,padType);                        % transmit data grid
            
            dopplerOut = dopplerChannel(txOut,fsamp,chanParams);    % send through channel
            chOut = awgn(dopplerOut,SNRdB,'measured');              % add noise
            
            rxWindow = chOut(1:numSamps);
            G = getG(M,N,chanEst,padLen,padType);
            y_otfs = ((G'*G)+n0*eye(Meff*N)) \ (G'*rxWindow); % LMMSE
            Xhat_otfs = helperOTFSdemod(y_otfs,M,padLen,0,padType); % OTFS demodulation
            
            XhatDataOTFS = pskdemod(Xhat_otfs,4,pi/4,OutputType="bit",OutputDataType="logical");
            [~,berOTFS] = biterr(Xdata,XhatDataOTFS);
            ber = ber + berOTFS;
        end
        ber = ber/100;
        berOTFSs(end+1) = ber;
    end
     
    semilogy(SNRdBs, berOTFSs);
    hold on;
end

grid minor;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs. SNR (OTFS)');
