clear;
close all;

M = 16;          % number of subcarriers
N = 6;          % number of subsymbols/frame
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

ddv1 = zeros(N, N, M, M);

for p = 3:chanParams.P
    for k = 1:N
        for l = 1:M
            for kp = 1:N
                tic
                for lp = 1:M
                    taup = chanParams.pathDelayTimes(p);
                    nup = chanParams.pathDopplerFreqs(p);
                    ddv1(k, kp, l, lp) = dd_crosstalk_coefficients(k-1,kp-1,l-1,lp-1,nup,taup,T,N,M);
                end
                toc
            end
        end
    end
end

for p = 3:chanParams.P
    rhop = chanParams.pathGains(p);
    taup = chanParams.pathDelayTimes(p);
    nup = chanParams.pathDopplerFreqs(p);
    hp = rhop*exp(1j*2*pi*nup*taup);

    ddv2 = dd_crosstalk_coefficientsv2(nup, taup, T, N, M);
end


