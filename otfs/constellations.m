clear;
close all;
%% generate constellations figure for OFDM and SFFT-OTFS

M = 64;          % number of subcarriers
N = 30;          % number of subsymbols/frame
df = 15e3;       % make this the frequency bin spacing of LTE
fc = 5e9;        % carrier frequency in Hz
padLen = 10;     % make this larger than the channel delay spread channel in samples
padType = 'ZP';  % this example requires ZP for ISI mitigation

% Configure paths
chanParams.pathDelays      = [0  5   8  ]; % number of samples that path is delayed
chanParams.pathGains       = [1  0.7 0.5]; % complex path gain
chanParams.pathDopplers    = [0 -3   5  ]; % Doppler index as a multiple of fsamp/MN
    
fsamp = M*df;            % sampling frequency at the Nyquist rate
Meff = M + padLen;       % number of samples per OTFS subsymbol
numSamps = Meff * N;     % number of samples per OTFS symbol
T = ((M+padLen)/(M*df)); % symbol time (seconds)

% Calculate the actual Doppler frequencies from the Doppler indices
chanParams.pathDopplerFreqs = chanParams.pathDopplers * 1/(N*T); % Hz

SNRdB = 40;

% Pilot generation and grid population
pilotBin = floor(N/2)+1;
Pdd = zeros(M,N);
Pdd(1,pilotBin) = exp(1i*pi/4); % populate just one bin to see the effect through the channel

% OTFS modulation
txOut = helperOTFSmod(Pdd,padLen,padType);

% Send the OTFS modulated signal through the channel
dopplerOut = dopplerChannel(txOut,fsamp,chanParams);

% Add white Gaussian noise
Es = mean(abs(pskmod(0:3,4,pi/4).^ 2));
n0 = Es/(10^(SNRdB/10));
chOut = awgn(dopplerOut,SNRdB,'measured');

for k = 1:length(chanParams.pathDelays)
    fprintf('Scatterer %d\n',k);
    fprintf('\tDelay = %5.2f us\n', 1e6*chanParams.pathDelays(k)/(Meff*df));
    fprintf('\tRelative Doppler shift = %5.0f Hz (%5.0f km/h)\n', ...
        chanParams.pathDopplerFreqs(k), (physconst('LightSpeed')*chanParams.pathDopplerFreqs(k)/fc)*(3600/1000));
end

% Get a sample window
rxIn = chOut(1:numSamps);

% OTFS demodulation
Ydd = helperOTFSdemod(rxIn,M,padLen,0,padType);

% LMMSE channel estimate in the delay-Doppler domain
Hdd = Ydd * conj(Pdd(1,pilotBin)) / (abs(Pdd(1,pilotBin))^2 + n0);

% figure;
% xa = 0:1:N-1;
% ya = 0:1:M-1;
% mesh(xa,ya,abs(Hdd));
% view([-9.441 62.412]);
% title('Delay-Doppler Channel Response H_{dd} from Channel Sounding');
% xlabel('Normalized Doppler');
% ylabel('Normalized Delay');
% zlabel('Magnitude');

[lp,vp] = find(abs(Hdd) >= 0.05);
chanEst.pathGains = diag(Hdd(lp,vp));   % get path gains
chanEst.pathDelays = lp - 1;            % get delay indices
chanEst.pathDopplers = vp - pilotBin;   % get Doppler indices


% Data generation
Xgrid = zeros(M,N);
Xdata = randi([0,1],2*M,N);
Xgrid(1:M,:) = pskmod(Xdata,4,pi/4,InputType="bit");

% Transmit pilots over all subcarriers and symbols to sound the channel
txOut = ofdmmod(exp(1i*pi/4)*ones(M,N),M,padLen);       % transmit pilots over the entire grid
dopplerOut = dopplerChannel(txOut,fsamp,chanParams);    % send through channel
chOut = awgn(dopplerOut,SNRdB,'measured');              % add noise
Yofdm = ofdmdemod(chOut(1:(M+padLen)*N),M,padLen);      % demodulate
Hofdm = Yofdm * conj(Pdd(1,pilotBin)) / (abs(Pdd(1,pilotBin))^2 + n0); % LMMSE channel estimate

% Transmit data over the same channel and use channel estimates to equalize
txOut = ofdmmod(Xgrid,M,padLen);                        % transmit data grid

dopplerOut = dopplerChannel(txOut,fsamp,chanParams);    % send through channel
chOut = awgn(dopplerOut,SNRdB,'measured');              % add noise

rxWindow = chOut(1:(M+padLen)*N);
Yofdm = ofdmdemod(rxWindow,M,padLen);                   % demodulate

Xhat_ofdm = conj(Hofdm) .* Yofdm ./ (abs(Hofdm).^2+n0); % equalize with LMMSE

constDiagOFDM = comm.ConstellationDiagram( ...
    'ReferenceConstellation',pskmod(0:3,4,pi/4), ...
    'XLimits',[-2 2], ...
    'YLimits',[-2 2], ...
    'Title','OFDM with Single-Tap FDE');
constDiagOFDM(Xhat_ofdm(:));

XhatDataOFDM = pskdemod(Xhat_ofdm,4,pi/4, ...
    OutputType="bit",OutputDataType="logical");         % decode
[~,berOFDM] = biterr(Xdata,XhatDataOFDM);
fprintf('OFDM BER with single-tap equalizer = %3.3e\n', berOFDM);

% OTFS modulation
txOut = helperOTFSmod(Xgrid,padLen,padType);

% Add channel and noise
dopplerOut = dopplerChannel(txOut,fsamp,chanParams);
chOut = awgn(dopplerOut,SNRdB,'measured');

% Form G matrix using channel estimates
G = getG(M,N,chanEst,padLen,padType);

rxWindow = chOut(1:numSamps);
y_otfs = ((G'*G)+n0*eye(Meff*N)) \ (G'*rxWindow); % LMMSE

Xhat_otfs = helperOTFSdemod(y_otfs,M,padLen,0,padType); % OTFS demodulation

constDiagOTFS = comm.ConstellationDiagram( ...
    'ReferenceConstellation',pskmod(0:3,4,pi/4), ...
    'XLimits',[-2 2], ...
    'YLimits',[-2 2], ...
    'Title','OTFS with Time-Domain LMMSE Equalization');
constDiagOTFS(Xhat_otfs(:));

XhatDataOTFS = pskdemod(Xhat_otfs,4,pi/4,OutputType="bit",OutputDataType="logical");
[~,berOTFS] = biterr(Xdata,XhatDataOTFS);
fprintf('OTFS BER with LMMSE equalization = %3.3e\n', berOTFS);
berOFDMs(end+1) = berOFDM;
berOTFSs(end+1) = berOTFS;
