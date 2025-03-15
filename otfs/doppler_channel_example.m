clear;
close all;
%% Doppler channel example
% Input and Output plotted in DD domain

M = 64;          % number of subcarriers
N = 30;          % number of subsymbols/frame
df = 15e3;       % make this the frequency bin spacing of LTE
fc = 5e9;        % carrier frequency in Hz
padLen = 0;     % make this larger than the channel delay spread channel in samples
padType = 'ZP';  % this example requires ZP for ISI mitigation
  
fsamp = M*df;            % sampling frequency at the Nyquist rate
Meff = M + padLen;       % number of samples per OTFS subsymbol
numSamps = Meff * N;     % number of samples per OTFS symbol
T = ((M+padLen)/(M*df)); % symbol time (seconds)

% Configure paths
    chanParams.pathDelays      = [0  5   8  ]; % number of samples that path is delayed
    chanParams.pathGains       = [1  0.7 0.5]; % complex path gain
    chanParams.pathDopplers    = [0 -3   5  ]; % Doppler index as a multiple of fsamp/MN
    chanParams.pathDopplerFreqs = chanParams.pathDopplers * 1/(N*T); % Hz  

% Calculate the actual Doppler frequencies from the Doppler indices
chanParams.pathDopplerFreqs = chanParams.pathDopplers * 1/(N*T); % Hz
for k = 1:length(chanParams.pathDelays)
    fprintf('Scatterer %d\n',k);
    fprintf('\tDelay = %5.2f us\n', 1e6*chanParams.pathDelays(k)/(Meff*df));
    fprintf('\tRelative Doppler shift = %5.0f Hz (%5.0f km/h)\n', ...
        chanParams.pathDopplerFreqs(k), (physconst('LightSpeed')*chanParams.pathDopplerFreqs(k)/fc)*(3600/1000));
end

% Pilot generation and grid population
pilotBin = floor(N/2)+1;
Pdd = zeros(M,N);
Pdd(1,pilotBin) = exp(1i*pi/4); 

txOut = izak(Pdd);

dopplerOut = dopplerChannel(txOut,fsamp,chanParams); 

rxIn = dopplerOut(1:numSamps);
Ydd = zak(rxIn, Meff);

figure;
xa = 0:1:N-1;
ya = 0:1:M-1;
mesh(xa,ya,abs(Pdd));
view([-9.441 62.412]);
title('Delay-Doppler Channel Input in DD domain');
xlabel('Normalized Doppler');
ylabel('Normalized Delay');
zlabel('Magnitude');

figure;
xa = 0:1:N-1;
ya = 0:1:M-1;
mesh(xa,ya,abs(Ydd));
view([-9.441 62.412]);
title('Delay-Doppler Channel Output in DD domain');
xlabel('Normalized Doppler');
ylabel('Normalized Delay');
zlabel('Magnitude');
