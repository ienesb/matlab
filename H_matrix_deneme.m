clear;
close all;

M = 64;          % number of subcarriers
N = 30;          % number of subsymbols/frame
df = 15e3;       % make this the frequency bin spacing of LTE
fc = 5e9;        % carrier frequency in Hz
padLen = 0;     % make this larger than the channel delay spread channel in samples
padType = 'NONE';  % this example requires ZP for ISI mitigation

% Configure paths
chanParams.pathDelays      = [0  5   8  ]; % number of samples that path is delayed
chanParams.pathGains       = [1  0.7 0.5]; % complex path gain
chanParams.pathDopplers    = [0 -3   5  ]; % Doppler index as a multiple of fsamp/MN
    
fsamp = M*df;            % sampling frequency at the Nyquist rate
Meff = M + padLen;       % number of samples per OTFS subsymbol
numSamps = Meff * N;     % number of samples per OTFS symbol
T = 1/df;                   % symbol time (seconds)

% Calculate the actual Doppler frequencies from the Doppler indices
chanParams.pathDopplerFreqs = chanParams.pathDopplers * 1/(N*T); % Hz
chanParams.pathDelayTimes = chanParams.pathDelays * 1/fsamp; % second


pilotBin = floor(N/2)+1;
Pdd = zeros(M,N);
Pdd(1,pilotBin) = exp(1i*pi/4);

% Doppler Channel implementation

% OTFS modulation
txOut = helperOTFSmod(Pdd,padLen,padType);

% Send the OTFS modulated signal through the channel
dopplerOut = dopplerChannel(txOut,fsamp,chanParams);

% Get a sample window
rxIn = dopplerOut(1:numSamps);

% OTFS demodulation
Ydd = helperOTFSdemod(rxIn,M,padLen,0,padType);

figure;
surf(abs(Pdd));

figure;
surf(abs(Ydd));


% H matrix implementation
Ptf = isfft(Pdd, M);
Ytf = zeros(M, N);

for n=0:N-1
    for m=0:M-1
        for np=0:N-1
            for mp=0:M-1
                H = getH(n, np, m, mp, chanParams, T);
                Ytf(m+1, n+1) = Ytf(m+1, n+1) + H*Ptf(mp+1, np+1);
            end
        end
    end
end

Ydd = sfft(Ytf, M);

figure;
surf(abs(Ydd));
