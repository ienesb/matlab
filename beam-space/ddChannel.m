%% DD Channel, H Channel ve Doppler Channel implementasyonları kıyaslaması. DD ve H aynı sonucu veriyor.
clear;
close all;

M = 64;          % number of subcarriers
N = 30;          % number of subsymbols/frame
df = 15e3;       % make this the frequency bin spacing of LTE
fc = 5e9;        % carrier frequency in Hz
padLen = 0;     % make this larger than the channel delay spread channel in samples
padType = 'NONE';  % this example requires ZP for ISI mitigation

fsamp = M*df;            % sampling frequency at the Nyquist rate
Meff = M + padLen;       % number of samples per OTFS subsymbol
numSamps = Meff * N;     % number of samples per OTFS symbol
T = 1/df;                   % symbol time (seconds)

Tframe = N*T;   % frame duration
W = M*df;       % bandwidth

% Configure paths
chanParams.pathDelays      = [0  5   8  ]; % number of samples that path is delayed
chanParams.pathGains       = [1  0.7 0.5]; % complex path gain
chanParams.pathDopplers    = [0 -3   5  ]; % Doppler index as a multiple of fsamp/MN
    
% Calculate the actual Doppler frequencies from the Doppler indices
chanParams.pathDopplerFreqs = chanParams.pathDopplers * 1/(N*T); % Hz
chanParams.pathDelayTimes = chanParams.pathDelays * 1/fsamp; % second

chanParams.P = length(chanParams.pathDelays);

% pilotBin = floor(N/2)+1;
% Xdd = zeros(M,N);
% Xdd(1,pilotBin) = exp(1i*pi/4);

% Data generation
Xdd = zeros(M,N);
Xdata = randi([0,1],2*M,N);
Xdd(1:M,:) = pskmod(Xdata,4,pi/4,InputType="bit");

%% Doppler Channel implementation

% OTFS modulation
txOut = OTFSmod(Xdd,padLen,padType);

% Send the OTFS modulated signal through the channel
dopplerOut = dopplerChannel(txOut,fsamp,chanParams);

% Get a sample window
rxIn = dopplerOut(1:numSamps);

% OTFS demodulation
Yddv1 = OTFSdemod(rxIn,M,padLen,0,padType);

figure; 
surf(abs(Xdd));

figure;
surf(abs(Yddv1));


%% H matrix implementation

Xtf = isfft(Xdd, M);
Ytf = zeros(M, N);

for n=0:N-1
    for m=0:M-1
        for np=0:N-1
            for mp=0:M-1
                H = getH(n, np, m, mp, chanParams, T);
                Ytf(m+1, n+1) = Ytf(m+1, n+1) + H*Xtf(mp+1, np+1);
            end
        end
    end
end

Yddv2 = sfft(Ytf, M);

figure;
surf(abs(Yddv2));


%% dd implementation

Yddv3 = zeros(M,N);

Xdd = Xdd.';
Yddv3 = Yddv3.';
for p = 1:chanParams.P
    rhop = chanParams.pathGains(p);
    taup = chanParams.pathDelayTimes(p);
    nup = chanParams.pathDopplerFreqs(p);
    hp = rhop*exp(1j*2*pi*nup*taup);

    dd = dd_crosstalk_coefficientsv2(nup, taup, T, N, M);

    dd = permute(dd, [1,3,2,4]);

    temp = reshape(Xdd, [1, 1, N, M]);
    
    temp = dd.*temp;
    
    temp = sum(sum(temp, 4), 3);

    Yddv3 = Yddv3 + hp*temp;
end
Xdd = Xdd.';
Yddv3 = Yddv3.';

figure;
surf(abs(Yddv3));


