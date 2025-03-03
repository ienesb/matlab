clear;
close all;

M = 64;          % number of subcarriers
N = 30;          % number of subsymbols/frame
df = 15e3;       % make this the frequency bin spacing of LTE
fc = 5e9;        % carrier frequency in Hz
padLen = 0;     % make this larger than the channel delay spread channel in samples
padType = 'NONE';  % this example requires ZP for ISI mitigation

Na = 64;
Nrf = 4;
Ns = 1;

B = 6;

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
chanParams.pathDelays      = [0  5   8  ]; % number of samples that path is delayed
chanParams.pathGains       = [1  0.7 0.5]; % complex path gain
chanParams.pathDopplers    = [0 -3   5  ]; % Doppler index as a multiple of fsamp/MN
chanParams.pathAoAs        = [0 -20 35].*(pi/180); % Doppler index as a multiple of fsamp/MN
    
% Calculate the actual Doppler frequencies from the Doppler indices
chanParams.pathDopplerFreqs = chanParams.pathDopplers * 1/(N*T); % Hz
chanParams.pathDelayTimes = chanParams.pathDelays * 1/fsamp; % second

chanParams.P = length(chanParams.pathDelays);

b = 1;

U = getU(Na, Nrf, B, b);
F = getF(Na, Ns);



pilotBin = floor(N/2)+1;
Xdd = zeros(M, N);
Xdd(1,pilotBin) = exp(1i*pi/4);
Xdd = Xdd.';

% Data generation
% Xdd = zeros(M,N);
% Xdata = randi([0,1],2*M,N);
% Xdd(1:M,:) = pskmod(Xdata,4,pi/4,InputType="bit");

Ydd = zeros(N, M, Nrf);

G = 0;
G1 = 0;

for p = 1:chanParams.P
    rhop = chanParams.pathGains(p);
    taup = chanParams.pathDelayTimes(p);
    nup = chanParams.pathDopplerFreqs(p);
    hp = rhop*exp(1j*2*pi*nup*taup);
    phip = chanParams.pathAoAs(p);

    a = array_response(phip, Na);

    dd = dd_crosstalk_coefficientsv2(nup, taup, T, N, M);
    ddm = dd2m(dd, M, N);
    dd = permute(dd, [1,3,2,4]);

    temp = reshape(Xdd, [1, 1, N, M]);
    temp = dd.*temp;
    temp = sum(sum(temp, 4), 3);

    C = U' * (a * a') * F;
    % G = G + kron(C, ddm) * hp;

    G = G + getGb(nup, taup, phip, b) * hp;

    for n = 1:Nrf
        Ydd(:, :, n) = Ydd(:, :, n) + temp *  C(n)*hp;
    end
end


% for n = 1:Nrf
%     figure;
%     surf(abs(Ydd(:, :, n).'));
% end

Xdd = Xdd.';
x = Xdd(:);

% Ydd = Ydd';
Ydd = permute(Ydd, [2, 1, 3]);

yv1 = Ydd(:);
yv2 = G*x;

diff = yv1 - yv2;

figure;
plot(abs(diff));

% Y = reshape(yv2, M, N);
% 
% figure;
% surf(abs(Y))


ks = 0:N-1;
ls = 0:M-1;

doppler_grids = ks ./ (N*T);
delays_grids = ls ./ (M*df);

bs = 1:B;
ns = (1:Nrf)';

angle_grids = theta_min + (bs-1)*fov/B + (2*ns-1)/2*fov/(B*Nrf);
angle_grids = angle_grids(:);

