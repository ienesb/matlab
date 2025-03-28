clear;
close all;

c = physconst('LightSpeed');

M = 16;
N = 32;
df = 120e3;
T = 1/df;
fc = 28e9;
lambda = c/fc;

NT = 8;
NR = 8;

sigma = 1;

targets.velocities = [20, 20, 20, 20, 20];
targets.ranges = [60 150 150 360 460];
targets.dopplers = 2*targets.velocities/lambda;
targets.delays = 2*targets.ranges/c;
targets.angles = [10, 10, 15, 10, 10] .*(pi/180);
targets.SNRs = [20, 15, 5, 25, 10];
targets.K = length(targets.ranges);
targets.gains = sqrt(db2pow(targets.SNRs) .* sigma^2);
% targets.gains = sqrt(db2pow(targets.SNRs) .* sigma^2) .* exp(1j*rand(1, targets.K)*360);

com.JR = 0; % dB
com.SNR = 10; % dB
com.angles = ones(1,11) .* (-30 * (pi/180));
com.K = length(com.angles);

SNR = db2pow(com.SNR);
JR = db2pow(com.JR);
a_multipath = randn(1, com.K-1);
P_multipath = sum(a_multipath.^2);

a0 = sqrt(JR * P_multipath);

total_power = a0^2 + P_multipath;
scale_factor = sqrt((SNR * sigma^2) / total_power);

a0 = a0 * scale_factor;
a_multipath = a_multipath * scale_factor;

com.gains = [a0, a_multipath];

% rho = 0.4;


Drad = 0;
for k = 1:targets.K
    alpha = targets.gains(k); % gain
    theta = targets.angles(k); % AoA
    atx = getatx(theta, NT);
    Drad = Drad + abs(alpha)^2 * (atx * atx');
end
Drad = Drad .* (1/sigma^2);

Dcom = 0;
for k = 1:com.K
    alpha = com.gains(k); % gain
    theta = com.angles(k); % AoA
    atx = getatx(theta, NT);
    Dcom = Dcom + abs(alpha)^2 * (atx * atx');
end
Dcom = Dcom .* (1/sigma^2);

% D = rho * Drad + (1-rho) * Dcom;
% [V,~] = eig(D);
% b_opt = V(:, end);
% 
% p = ones(N*M, 1);
% Ws = p * b_opt.';
% Ws = permute(Ws, [1,3,2]);
% Ws = reshape(Ws, [N, M, NT]);
% Ws = Ws.*sqrt(NT);

angles = (-80:80) .* (pi/180);
atxs = getatx(angles, NT);
rhos = [0.1 0.4 0.9];

figure;
hold on;

for idx = 1:length(rhos)
    rho = rhos(idx);
    D = rho * Drad + (1-rho) * Dcom;
    [V,e] = eig(conj(D));
    b_opt = V(:, end);
    beampattern = abs(atxs.' * b_opt).^2;
    beampattern = pow2db(beampattern);
    plot(angles .* (180/pi), beampattern);
end

grid minor;
xlabel("Angle");
ylabel("|a_{tx}^T \beta(\rho)|")
legend("\rho = 0.1", "\rho = 0.4", "\rho = 0.9")
