clear;
close all;

M = 64;          % number of subcarriers
N = 30;          % number of subsymbols/frame
df = 15e3;       % make this the frequency bin spacing of LTE
fc = 5e9;        % carrier frequency in Hz
padLen = 0;     % make this larger than the channel delay spread channel in samples
padType = 'NONE';  % this example requires ZP for ISI mitigation

Na = 64;
Nrf = 1;
Ns = 1;

B = 1;

theta_min = -pi/4;
theta_max = pi/4;
fov = theta_max - theta_min;

fsamp = M*df;            % sampling frequency at the Nyquist rate
Meff = M + padLen;       % number of samples per OTFS subsymbol
numSamps = Meff * N;     % number of samples per OTFS symbol
T = 1/df;                   % symbol time (seconds)

Tframe = N*T;   % frame duration
W = M*df;       % bandwidth

% Configure paths
chanParams.pathDelays      = [0 5   8  ]; % number of samples that path is delayed
chanParams.pathGains       = [1 0.7 0.5]; % complex path gain
chanParams.pathDopplers    = [0 -3   5  ]; % Doppler index as a multiple of fsamp/MN
    
% Calculate the actual Doppler frequencies from the Doppler indices
chanParams.pathDopplerFreqs = chanParams.pathDopplers * 1/(N*T); % Hz
chanParams.pathDelayTimes = chanParams.pathDelays * 1/fsamp; % second

chanParams.P = length(chanParams.pathDelays);
    
% pilotBin = floor(N/2)+1;
% Xdd = zeros(M, N);
% Xdd(1,pilotBin) = exp(1i*pi/4);

% Data generation
Xdd = zeros(M,N);
Xdata = randi([0,1],2*M,N);
Xdd(1:M,:) = pskmod(Xdata,4,pi/4,InputType="bit");

x = Xdd(:);
y = zeros(N*M, 1);
     
G = 0;
for p = 1:chanParams.P
    rhop = chanParams.pathGains(p);
    taup = chanParams.pathDelayTimes(p);
    nup = chanParams.pathDopplerFreqs(p);
    hp = rhop*exp(1j*2*pi*nup*taup);

    dd = dd_crosstalk_coefficientsv2(nup, taup, T, N, M);
    ddm = dd2m(dd, M, N);
    
    % G matrix
    G = G + hp * ddm;
end
    
y = G*x;

ks = -(N-1):N-1;
ls = 0:M-1;

doppler_grids = ks ./ (N*T);
delay_grids = ls ./ (M*df);
  
Ss = zeros(length(doppler_grids), length(delay_grids));

idx_doppler = 1;
idx_delay = 1;


for nu = doppler_grids
    tic
    for tau = delay_grids
        dd = dd_crosstalk_coefficientsv2(nu, tau, T, N, M);
        ddm = dd2m(dd, M, N);
        G = ddm;
        temp = G*x;
        denum = norm(temp).^2;
        temp = y'*temp;        
        num = abs(temp).^2;
        S = num / denum;
        Ss(idx_doppler, idx_delay) = S;
        idx_delay = idx_delay + 1;
    end
    toc
    idx_delay = 1;
    idx_doppler = idx_doppler + 1;
end
