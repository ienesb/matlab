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
chanParams.pathDelays      = [0 5   8  ]; % number of samples that path is delayed
chanParams.pathGains       = [1 0.7 0.5]; % complex path gain
chanParams.pathDopplers    = [0 -3   5  ]; % Doppler index as a multiple of fsamp/MN
chanParams.pathAoAs        = [0 -20 35].*(pi/180);
    
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
ys = zeros(Nrf*N*M, B);

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
        
        % G matrix
        C = U' * (a * a') * F;
        G = G + hp * kron(C, ddm);
    end
    
    ys(:, b) = G*x;
end

ks = 0:N-1;
ls = 0:M-1;

doppler_grids = ks ./ (N*T);
delay_grids = ls ./ (M*df);

bs = 1:B;
ns = (1:Nrf)';

angle_grids = theta_min + (bs-1)*fov/B + (2*ns-1)/2*fov/(B*Nrf);
angle_grids = angle_grids(:)';

a = array_response(angle_grids, Na);

a = permute(a, [1,3,2]);
aH = conj(permute(a, [2,1,3]));

aaH = pagemtimes(a, aH);
Cs = zeros(Nrf, Ns, length(angle_grids), B);
for b = 1:6
    U = getU(Na, Nrf, B, b); 
    Cs(:, :, :, b) = pagemtimes(U', pagemtimes(aaH, F));
end
  
% Ss = zeros(length(angle_grids), length(doppler_grids), length(delay_grids));
% 
% idx_doppler = 1;
% idx_delay = 1;
% tic
% for nu = doppler_grids
%     for tau = delay_grids
%         num = zeros(length(angle_grids), 1);
%         denum = zeros(length(angle_grids), 1);
%         for b = 1:6
%             C = Cs(:,:,:,b);
%             dd = dd_crosstalk_coefficientsv2(nu, tau, T, N, M);
%             ddm = dd2m(dd, M, N);
%             Gs = my_kron_nd_fast(C, ddm);
%             temp = pagemtimes(Gs, x);
%             d = pagenorm(temp).^2;
%             d = d(:);
%             denum = denum + d;
%             temp = pagemtimes(ys(:, b)', temp);
%             num = num + temp(:);
%         end
%         num = abs(num).^2;
%         S = num ./ denum;
%         Ss(:, idx_doppler, idx_delay) = S;
%         idx_delay = idx_delay + 1;
%     end
%     idx_delay = 1;
%     idx_doppler = idx_doppler + 1;
% end
% toc

dtau = 1/(M*df);
dnu = 1/(N*T);
dphi = fov/(B*Nrf);

tau = delay_grids(12);
nu = doppler_grids(12);
phi = angle_grids(12);

taus = tau + [-3*dtau -2*dtau 2*dtau 3*dtau];
nus = nu + [-3*dnu -2*dnu 2*dnu 3*dnu];
phis = phi + [-3*dphi -2*dphi 2*dphi 3*dphi];