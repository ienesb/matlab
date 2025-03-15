clear;
% close all;

M = 64;          % number of subcarriers
N = 30;          % number of subsymbols/frame
df = 15e3;       % make this the frequency bin spacing of LTE
fc = 5e9;        % carrier frequency in Hz
padLen = 0;     % make this larger than the channel delay spread channel in samples
padType = 'NONE';  % this example requires ZP for ISI mitigation

Na = 64;
Nrf = 4;
Ns = 1;


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
chanParams.pathDelays      = [5   8  ]; % number of samples that path is delayed
chanParams.pathGains       = [0.7 0.5]; % complex path gain
chanParams.pathDopplers    = [-3   5  ]; % Doppler index as a multiple of fsamp/MN
chanParams.pathAoAs        = [-20 35].*(pi/180);
    
% Calculate the actual Doppler frequencies from the Doppler indices
chanParams.pathDopplerFreqs = chanParams.pathDopplers * 1/(N*T); % Hz
chanParams.pathDelayTimes = chanParams.pathDelays * 1/fsamp; % second

chanParams.P = length(chanParams.pathDelays);
    
pilotBin = floor(N/2)+1;
Xdd = zeros(M, N);
Xdd(1,pilotBin) = exp(1i*pi/4);
x = Xdd(:);

results = zeros(3, 11);
Bidx = 1;
for B = [6 10 18]
B
ys = zeros(N*M*Nrf,B);
rs = zeros(N*M*Nrf,B);

F = getF(Na, Ns);

for b = 1:B
    U = getU(Na, Nrf, B, b);   
    
    G = 0;
    for p = 1:chanParams.P
        rhop = chanParams.pathGains(p);
        taup = chanParams.pathDelayTimes(p);
        nup = chanParams.pathDopplerFreqs(p);
        hp = rhop*exp(1j*2*pi*nup*taup);
        phip = chanParams.pathAoAs(p);
    
        a = array_response(phip, Na);
    
        dd = dd_crosstalk_coefficientsv2(nup, taup, T, N, M);
        ddm = dd2m(dd, M, N);
    
        C = hp * U' * (a * a') * F;
        G = G + kron(C, ddm);
    end

    % sigma_w = 1;
    % w = normrnd(0, sigma_w, [N*M*Nrf,1]);
    
    r = G*x;
    rs(:, b) = r; 
end

ks = 0:N-1;
ls = 0:M-1;

doppler_grids = ks ./ (N*T);
delays_grids = ls ./ (M*df);

bs = 1:B;
ns = (1:Nrf)';

angle_grids = theta_min + (bs-1)*fov/B + (2*ns-1)/2*fov/(B*Nrf);
angle_grids = angle_grids(:)';



test_indexes = [28;
                6;
                7];

gt = [1 1 0];
% result = zeros(1,11);
i = 1;

for sigma_w = 0:5:50
    sigma_w
    for epoch = 1:20
        w = normrnd(0, sigma_w, [N*M*Nrf, B]);
        ys = rs + w;
    
        for idx = test_indexes
            nu = doppler_grids(idx(1));
            tau = delays_grids(idx(2));
            phi = angle_grids(idx(3));
        
            test_grid = generateTestGrid(nu, tau, phi, N, M, T, B, fov, Nrf, 5);
            th = 0;
            for test = test_grid
                nu_test = test(1);
                tau_test = test(2);
                phi_test = test(3);
                th = th + getS(nu_test, tau_test, phi_test, ys, x, B);
            end
            th = th/100;
            if getS(nu, tau, phi, ys, x, B) > th
                % result(i) = result(i) + 1;
                results(Bidx, i) = results(Bidx, i) + 1;
            end
    
        end
    end
    i = i + 1;
end
Bidx = Bidx + 1;
end