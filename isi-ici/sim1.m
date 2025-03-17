clear;
close all;

c = physconst('LightSpeed');

M = 16;
N = 32;
df = 480e3;
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
targets.gains = sqrt(10.^(targets.SNRs/10) .* sigma^2) .* exp(1j*rand(1, targets.K)*360);

Ws = zeros(N, M, NT);
W = randi([1, NT], N, M);

for idx = 1:NT
    temp = W;
    temp(temp~=idx) = 0;
    temp(temp==idx) = 1;
    Ws(:, :, idx) = temp;
end

Xdata = randi([0, 63], N*M, 1);
XDD = qammod(Xdata, 64, UnitAveragePower=true);
XDD = reshape(XDD, N, M);
XDDs = XDD .* Ws;

FM = dftmtx(M)/sqrt(M);
FN = dftmtx(N)/sqrt(N);
F = dftmtx(N*M)/sqrt(N*M);

Gtx = eye(N,N);
% Gtx = eye(N,N) ./ sqrt(T);

ss = pagemtimes(Gtx, pagemtimes(XDDs, FM'));
S = reshape(ss, N*M, NT);

Y = 0;
for k = 1:targets.K
    tau = targets.delays(k); % delay
    nu = targets.dopplers(k); % doppler
    alpha = targets.gains(k); % gain
    theta = targets.angles(k); % AoA

    B = getB(tau, df, N, M);
    C = getC(nu, T, N, M);

    atx = getatx(theta, NT);
    arx = getarx(theta, NR);

    Y = Y + alpha * C * F' * B * F * S * atx * arx';
end

noise = (randn(N*M, NR) + 1i*randn(N*M, NR)) * (sigma/sqrt(2));
Y_with_noise = Y + noise;

% figure 5
ranges = 0:500;
delays = 2*ranges/c;
velocity = 20;
nu = 2*velocity/lambda;

range_profile = zeros(size(ranges));

for idx = 1:length(delays)
    tau = delays(idx);
    B = getB(tau, df, N, M);
    C = getC(nu, T, N, M);
    range_profile(idx) = norm(S' * F' * B' * F * C' * Y_with_noise, "fro")^2;
end
range_profile = mag2db(range_profile);

Nc = 100;
cfar = phased.CFARDetector;
cfar.ThresholdFactor = "Custom";
cfar.CustomThresholdFactor = 1.05;
cfar.ThresholdOutputPort = true;
cfar.NumTrainingCells = Nc;

[detected,th] = cfar(range_profile.',1:length(range_profile));


figure; hold on;
plot(ranges, range_profile);
plot(ranges, th);
xline(60, '--', 'LineWidth', 1.5);
xline(150, '--', 'LineWidth', 1.5);
xline(360, '--', 'LineWidth', 1.5);
xline(460, '--', 'LineWidth', 1.5);
xlabel("Range (m)");
ylabel("Magnitude (dB)");
legend("range profile", "cfar threshold", "ground truth");
grid minor;

% figure 6
% ranges = 0:500;
% delays = 2*ranges/c;
% velocities = -100:100;
% dopplers = 2*velocities/lambda;
% 
% profiles = zeros(length(ranges), length(velocities));
% 
% for idx1 = 1:length(delays)
%     for idx2 = 1:length(dopplers)
%         tau = delays(idx1);
%         nu = dopplers(idx2);
% 
%         B = getB(tau, df, N, M);
%         C = getC(nu, T, N, M);
%         profiles(idx1, idx2) = norm(S' * F' * B' * F * C' * Y, "fro")^2;
%     end
% end