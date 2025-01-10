clear;
close all;

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


% Pilot generation and grid population
pilotBin = floor(N/2)+1;
Pdd = zeros(M,N);
Pdd(1,pilotBin) = exp(1i*pi/4); % populate just one bin to see the effect through the channel

figure;

for c=0:0.25:1

% Configure paths
    chanParams.pathDelays      = [0  5   8  ]; % number of samples that path is delayed
    chanParams.pathGains       = [1  0.7*c 0.5*c]; % complex path gain
    chanParams.pathDopplers    = [0 -3   5  ]; % Doppler index as a multiple of fsamp/MN
    chanParams.pathDopplerFreqs = chanParams.pathDopplers * 1/(N*T); % Hz
    berOFDMs = [];
    
    SNRdBs = 0:1:40;
    for SNRdB = SNRdBs
        ber = 0;
        for i=1:100
            % Channel Estimation
            % Transmit pilots over all subcarriers and symbols to sound the channel
            txOut = ofdmmod(exp(1i*pi/4)*ones(M,N),M,padLen);       % transmit pilots over the entire grid
            dopplerOut = dopplerChannel(txOut,fsamp,chanParams);    % send through channel
            Es = mean(abs(pskmod(0:3,4,pi/4).^ 2));
            n0 = Es/(10^(SNRdB/10));
            chOut = awgn(dopplerOut,SNRdB,'measured');              % add noise
            Yofdm = ofdmdemod(chOut(1:(M+padLen)*N),M,padLen);      % demodulate
            Hofdm = Yofdm * conj(Pdd(1,pilotBin)) / (abs(Pdd(1,pilotBin))^2 + n0); % LMMSE channel estimate
        
            % Data generation
            Xgrid = zeros(M,N);
            Xdata = randi([0,1],2*M,N);
            Xgrid(1:M,:) = pskmod(Xdata,4,pi/4,InputType="bit");
            
            % Transmit data over the same channel and use channel estimates to equalize
            txOut = ofdmmod(Xgrid,M,padLen);                        % transmit data grid
            
            dopplerOut = dopplerChannel(txOut,fsamp,chanParams);    % send through channel
            chOut = awgn(dopplerOut,SNRdB,'measured');              % add noise
            
            rxWindow = chOut(1:(M+padLen)*N);
            Yofdm = ofdmdemod(rxWindow,M,padLen);                   % demodulate
            
            Xhat_ofdm = conj(Hofdm) .* Yofdm ./ (abs(Hofdm).^2+n0); % equalize with LMMSE
            
            XhatDataOFDM = pskdemod(Xhat_ofdm,4,pi/4, ...
                OutputType="bit",OutputDataType="logical");         % decode
            [~,berOFDM] = biterr(Xdata,XhatDataOFDM);
            ber = ber + berOFDM;
        end
        ber = ber/100;
        berOFDMs(end+1) = ber;
    end
     
    semilogy(SNRdBs, berOFDMs);
    hold on;
end

grid minor;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs. SNR (OFDM)');
